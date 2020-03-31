# functions.micsim.R
# function for initial population
# March 2020

# retain this function even if not using a random pop; it's required for a randomly generated immigrant population
getRandInitState <- function(birthDate, starting.probs = c(0.90,0.05,0.05)){ 
  age <- trunc(as.numeric(simHorizon[1] - birthDate)/365.25)
  s1 <- sample(c('f','m'), 1) # 50:50 gender split
  s2 <- '0' # ifelse(age<=18, fert[1], sample(fert, 1))
  s3 <- sample(c("S","E",'I1'), size=1, prob = starting.probs) # 
  initState <- paste(c(s1, s3), collapse="/") # removed s2, no fertility
  return(initState)
}

