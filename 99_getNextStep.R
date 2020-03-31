# getNextStep.R
# uses Vectorize in integrate
# March 2020

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
