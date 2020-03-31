micSim = function (initPop, immigrPop = NULL, transitionMatrix, absStates = NULL, 
          initStates = c(), initStatesProb = c(), maxAge = 99, simHorizon, 
          fertTr = c(), dateSchoolEnrol = "09/01") 
{
  if (is.null(initPop)) 
    stop("No starting population has been defined.")
  if (!is.null(initPop)) {
    if (paste(colnames(initPop), collapse = "/") != "ID/birthDate/initState") 
      stop("Matrix specifying the starting population has not been defined properly.")
  }
  if (!is.null(immigrPop)) {
    if (paste(colnames(immigrPop), collapse = "/") != "ID/immigrDate/birthDate/immigrInitState") 
      stop("Matrix specifying immigrants has not been defined properly.")
  }
  if (is.null(transitionMatrix)) 
    stop("Matrix defining transition pattern und functions has not been defined properly.")
  if (maxAge <= 0) 
    stop("The maximal age until which individual life courses are simulated should exceed zero.")
  if (length(simHorizon) != 2) 
    stop("The simulation horizon has not been defined properly.")
  if (class(simHorizon)[1] != "dates" & class(simHorizon)[2] != 
      "dates") 
    stop("The simulation horizon has not been defined properly.")
  if (is.null(absStates)) 
    absStates <- setdiff(colnames(transitionMatrix), rownames(transitionMatrix))
  if (length(fertTr) > 0) {
    if (is.null(initStates) | is.null(initStatesProb)) 
      stop("For children potentially born during simulation no inital state(s) and/or corresponding occurrence probabilities have been defined.")
  }
  if (length(dateSchoolEnrol) == 0) {
    dateSchoolEnrol <- "09/01"
  }
  queue <- matrix(NA, ncol = 6, nrow = 0)
  colnames(queue) <- c("ID", "currTime", "currState", "currAge", 
                       "nextState", "timeToNextState")
  t.clock <- as.numeric(simHorizon[1])
  transitions <- matrix(NA, ncol = 5, nrow = 0)
  colnames(transitions) <- c("ID", "From", "To", "transitionTime", 
                             "transitionAge")
  isLeapYear <- function(year) {
    if (((year%%4 == 0) & (year%%100 != 0)) || (year%%400 == 
                                                0)) 
      return(TRUE)
    return(FALSE)
  }
  getAgeInYears <- function(days) {
    date <- chron(days)
    c.y <- as.numeric(as.character(years(date)))
    completeYears <- c.y - 1970
    daysInCY <- ifelse(isLeapYear(c.y), 366, 365)
    fracCY <- as.numeric(date - chron(paste(1, "/", 1, "/", 
                                            c.y), format = c(dates = "d/m/y")) + 1)/daysInCY
    return(completeYears + fracCY)
  }
  isBirthEvent <- function(currState, destState) {
    if (length(fertTr) == 0) 
      return(FALSE)
    fert <- matrix(unlist(strsplit(fertTr, split = "->")), 
                   ncol = 2)
    cS <- unlist(strsplit(currState, "/"))
    if (!("f" %in% cS)) 
      return(FALSE)
    dS <- unlist(strsplit(destState, "/"))
    for (i in 1:nrow(fert)) {
      ff <- fert[i, ]
      oS <- unlist(strsplit(ff[1], "/"))
      bS <- unlist(strsplit(ff[2], "/"))
      if (!(F %in% (oS %in% cS)) & !(F %in% (bS %in% dS))) {
        return(TRUE)
      }
    }
    return(FALSE)
  }
  addNewNewborn <- function(birthTime) {
    birthState <- sample(apply(initStates, 1, paste, collapse = "/"), 
                         size = 1, replace = T, prob = initStatesProb)
    if (is.null(immigrPop)) {
      id <- as.numeric(max(as.numeric(initPop[, "ID"]))) + 
        1
    }
    else {
      id <- as.numeric(max(c(as.numeric(immigrPop[, "ID"]), 
                             as.numeric(initPop[, "ID"])))) + 1
    }
    birthDate <- dates(chron(birthTime, format = c(dates = "d/m/Y", 
                                                   times = "h:m:s"), out.format = c(dates = "d/m/year", 
                                                                                    times = "h:m:s")))
    newInd <- c(id, as.character(birthDate), birthState)
    initPop <<- rbind(initPop, newInd)
    nE <- getNextStep(c(id, birthState, 0, birthTime))
  }
  isSchoolEnrolment <- function(currState, destState) {
    enrol <- F
    if (T %in% ("no" %in% unlist(strsplit(currState, "/"))) & 
        T %in% ("low" %in% unlist(strsplit(destState, "/")))) 
      enrol <- T
    return(enrol)
  }
  getNextStep <- function(inp, isIMInitEvent = F) {
    id <- as.numeric(unlist(inp[1]))
    currState <- as.character(unlist(inp[2]))
    currAge <- as.numeric(unlist(inp[3]))
    calTime <- as.numeric(unlist(inp[4]))
    lagToWaitingTime <- ifelse(isIMInitEvent, (as.numeric(calTime) - 
                                                 as.numeric(simHorizon[1]))/365.25, 0)
    ageInYears <- getAgeInYears(currAge)
    possTr <- transitionMatrix[which(rownames(transitionMatrix) %in% 
                                       currState), ]
    possTr <- possTr[which(possTr != 0)]
    nextEventMatrix <- matrix(0, ncol = 2, nrow = length(possTr))
    ranMaxAge <- maxAge - ageInYears
    ranMaxYear <- (as.numeric(simHorizon[2]) - calTime)/365.25
    ran <- min(ranMaxYear, ranMaxAge)
    ranAge <- c(ageInYears, ageInYears + ranMaxAge)
    ranYear <- chron(c(calTime, as.numeric(calTime) + ran * 
                         365.25))
    historiesInd <- transitions[as.numeric(transitions[, 
                                                       "ID"]) %in% id & as.numeric(transitions[, "transitionTime"]) <= 
                                  calTime, , drop = F]
    initPopInd <- initPop[as.numeric(initPop[, "ID"]) %in% 
                            id, ]
    birthTime <- initPopInd["birthDate"]
    initState <- as.character(unlist(initPopInd["initState"]))
    if (as.numeric(birthTime) < as.numeric(simHorizon[1]) | 
        id %in% as.numeric(immigrPop[, "ID"])) {
      dur <- rbind(c(initState, NA), cbind(historiesInd[, 
                                                        "To"], historiesInd[, "transitionTime"]))
      dur <- cbind(dur, c(diff(as.numeric(dur[, 2])), 
                          0))
      colnames(dur) <- c("TransitionTo", "AtTime", "durUntil")
      dur[which(is.na(dur[, "AtTime"])), "durUntil"] <- NA
    }
    else {
      birthTime <- initPop[as.numeric(initPop[, "ID"]) %in% 
                             id, "birthDate"]
      dur <- rbind(c(initState, birthTime), cbind(historiesInd[, 
                                                               "To"], historiesInd[, "transitionTime"]))
      dur <- cbind(dur, c(diff(as.numeric(dur[, 2])), 
                          0))
      colnames(dur) <- c("TransitionTo", "AtTime", "durUntil")
    }
    for (i in 1:length(possTr)) {
      tr <- possTr[i]
      destState <- names(tr)
      cS <- unlist(strsplit(currState, "/"))
      dS <- unlist(strsplit(destState, "/"))
      covToCh <- which((cS == dS) == F)
      durSinceLastCovCh <- Inf
      if (length(covToCh) == 1) {
        covHist <- do.call(rbind, sapply(dur[, "TransitionTo"], 
                                         strsplit, split = "/"))[, covToCh]
        idd <- which(covHist == cS[covToCh])
        if (length(idd) > 1) {
          if (F %in% (diff(idd) == 1)) {
            y <- rev(idd)[c(-1, diff(rev(idd))) == -1]
            idd <- rev(y)[c(diff(rev(y)), 1) == 1]
          }
        }
        durSinceLastCovCh <- sum(as.numeric(dur[idd, 
                                                "durUntil"]))
        if (is.na(durSinceLastCovCh)) 
          durSinceLastCovCh <- 0
      }
      if (length(covToCh) > 1 & (!destState %in% absStates)) {
        cat("Recognized a possible transition implying a change of two or more covariates.", 
            "Concerning the derivation of the time being elapsed since the last transition this feature is not yet implemented.", 
            "Current State: ", currState, " -> Possible transition to ", 
            destState, "\n")
      }
      indRateFctDET <- function(x) {
        res <- eval(do.call(tr, args = list(age = trunc(ageInYears) + 
                                              x, calTime = trunc(1970.001 + calTime/365.25) + 
                                              x, duration = trunc(durSinceLastCovCh/365.25) + 
                                              x)))
        return(res)
      }
      ranAccuracyInDays <- (0:(trunc(ran * 365.25) + 1))/365.25
      detE <- indRateFctDET(ranAccuracyInDays)
      daysToTrInYears <- (which(detE == Inf)[1] - 1)/365.25
      if (Inf %in% detE) {
        timeToNext <- daysToTrInYears
      }
      else {
        u <- -log(1 - runif(1))
        indRateFct <- function(x) {
          ageIn <- ageInYears + x
          calIn <- 1970.001 + calTime/365.25 + x
          durIn <- durSinceLastCovCh/365.25 + x
          res <- eval(do.call(tr, args = list(age = ageIn, 
                                              calTime = calIn, duration = durIn)))
          if (TRUE %in% (res < 0)) 
            stop("I have found negative rate value/s for transition: ", 
                 tr, "\n\n                 This is implausible. Please check this. Simulation has been stopped.\n")
          return(res)
        }
        if (sum(indRateFct(0:ran)) == 0) {
          intHaz <- 0
        }
        else {
          intHaz <- try(integrate(Vectorize(indRateFct), lower = 0, 
                                  upper = ran)$value, silent = TRUE)
          if (inherits(intHaz, "try-error")) {
            intHaz <- integrate(Vectorize(indRateFct), lower = 0, 
                                upper = ran, stop.on.error = FALSE, rel.tol = 0.01)$value
          }
        }
        if (u <= intHaz) {
          invHazFct <- function(x) {
            try.res <- try(integrate(Vectorize(indRateFct), lower = 0, 
                                     upper = x)$value - u, silent = TRUE)
            if (inherits(try.res, "try-error")) {
              try.res <- integrate(Vectorize(indRateFct), lower = 0, 
                                   upper = x, stop.on.error = FALSE, rel.tol = 0.01)$value - 
                u
            }
            return(try.res)
          }
          timeToNext <- uniroot(invHazFct, interval = c(0, 
                                                        ran))$root
        }
        else {
          timeToNext <- Inf
        }
      }
      nextEventMatrix[i, 1] <- destState
      nextEventMatrix[i, 2] <- (timeToNext + lagToWaitingTime) * 
        365.25
    }
    nE <- nextEventMatrix[which(nextEventMatrix[, 2] == 
                                  min(as.numeric(nextEventMatrix[, 2]))), , drop = F]
    if (dim(nE)[1] > 1) 
      nE <- nE[1, , drop = F]
    if (nE[1, 2] != Inf) {
      tt <- chron(as.numeric(calTime) + as.numeric(nE[1, 
                                                      2]) - ageInYears%%1, out.format = c(dates = "d/m/year", 
                                                                                          times = "h/m/s"))
      if (isSchoolEnrolment(currState, nE[1, 1])) {
        enYear <- years(tt)
        if (as.numeric(months(tt)) <= 9) {
          enDate <- chron(paste(enYear, dateSchoolEnrol, 
                                sep = "/"), format = c(dates = "y/m/d"), 
                          out.format = c(dates = "d/m/year"))
        }
        else {
          enYear <- as.numeric(as.character(enYear))
          enDate <- chron(paste(enYear, dateSchoolEnrol, 
                                sep = "/"), format = c(dates = "y/m/d"), 
                          out.format = c(dates = "d/m/year"))
        }
        diffToEn <- as.numeric(enDate - tt)
        nE[1, 2] <- as.numeric(nE[1, 2]) + diffToEn
      }
      queue <<- rbind(queue, c(id, t.clock, currState, 
                               currAge - lagToWaitingTime * 365.25, nE[1, 1], 
                               nE[1, 2]))
    }
    return(nE)
  }
  cat("Initialization ... \n")
  IN <- data.frame(ID = initPop[, "ID"], currState = initPop[, 
                                                             "initState"], age = (simHorizon[1] - initPop[, "birthDate"]), 
                   calTime = rep(as.numeric(simHorizon[1]), dim(initPop)[1]), 
                   stringsAsFactors = FALSE)
  init <- apply(IN, 1, getNextStep)
  if (!is.null(immigrPop)) {
    IM <- data.frame(ID = immigrPop[, "ID"], currState = immigrPop[, 
                                                                   "immigrInitState"], age = immigrPop[, "immigrDate"] - 
                       immigrPop[, "birthDate"], calTime = as.numeric(immigrPop[, 
                                                                                "immigrDate"]), stringsAsFactors = FALSE)
    immigrInitPop <- immigrPop[, c("ID", "birthDate", "immigrInitState")]
    colnames(immigrInitPop)[3] <- "initState"
    initPop <- rbind(initPop, immigrInitPop)
    imit <- apply(IM, 1, getNextStep, isIMInitEvent = T)
  }
  cat("Simulation is running ... \n")
  currYear <- as.numeric(as.character(years(simHorizon[1])))
  while (dim(queue)[1] > 0 & t.clock <= as.numeric(simHorizon[2])) {
    queue <- queue[order(as.numeric(queue[, "currTime"]) + 
                           as.numeric(queue[, "timeToNextState"])), , drop = F]
    indS <- queue[1, ]
    queue <- queue[-1, , drop = F]
    t.clock <- as.numeric(indS["currTime"]) + as.numeric(indS["timeToNextState"])
    cY <- as.numeric(as.character(years(t.clock)))
    if (t.clock > as.numeric(simHorizon[2])) 
      break
    if (cY > currYear) {
      cat("Year: ", cY, "\n")
      currYear <- cY
    }
    age <- as.numeric(indS["currAge"]) + as.numeric(indS["timeToNextState"])
    transitions <- rbind(transitions, c(indS[c("ID", "currState", 
                                               "nextState")], t.clock, age))
    if (!indS["nextState"] %in% absStates) {
      if (isBirthEvent(indS["currState"], indS["nextState"])) {
        addNewNewborn(t.clock)
      }
      res <- getNextStep(c(indS[c("ID", "nextState")], 
                           age, t.clock))
    }
  }
  transitions <- transitions[order(as.numeric(transitions[, 
                                                          1])), , drop = F]
  if (nrow(transitions) == 0) {
    transitionsOut <- data.frame(ID = initPop[, "ID"], From = rep(NA, 
                                                                  nrow(initPop)), To = rep(NA, nrow(initPop)), transitionTime = rep(NA, 
                                                                                                                                    nrow(initPop)), transitionAge = rep(NA, nrow(initPop)), 
                                 stringsAsFactors = FALSE)
    cat("Simulation has finished.\n")
    cat("Beware that along the simulation horizon the individual/s considered do/es not experience any transition/s.\n")
    cat("------------------\n")
  }
  else {
    cat("Simulation has finished.\n------------------\n")
    transitionsOut <- data.frame(ID = transitions[, "ID"], 
                                 From = transitions[, "From"], To = transitions[, 
                                                                                "To"], transitionTime = dates(chron(as.numeric(transitions[, 
                                                                                                                                           "transitionTime"]), out.format = c(dates = "d/m/year", 
                                                                                                                                                                              times = "h:m:s"))), transitionAge = round(getAgeInYears(as.numeric(transitions[, 
                                                                                                                                                                                                                                                             "transitionAge"])), 2), stringsAsFactors = FALSE)
  }
  pop <- merge(initPop, transitionsOut, all = T, by = "ID")
  pop <- pop[order(as.numeric(pop[, 1])), ]
  return(pop)
}