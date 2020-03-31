# 1_microsim_covid.R
# run the COVID model using microsimulation based on a model using ODEs
# March 2020
#library(MicSim) # using my slightly modified versions
library(chron)
library(snowfall)
library(snow)
library(rlecuyer)
#library(stringr)
library(dplyr)
source('99_functions.micsim.R') # functions to start micsim
source('99_functions.micsim.rates.R') # functions for transition rates
source('99_micSim.R') # adapted from MicSim by me (Vectorize)
source('99_getNextStep.R') # adapted from MicSim by me (Vectorize)
source('https://github.com/cran/MicSim/raw/master/R/auxFctMicSim.r') # from github repo of MicSim

# save results here
if(exists('outfile')==FALSE){outfile = 'data/simresults.RData'}

# key parameters
IncubPeriod = 5 # Incubation period, days
DurMildInf = 6   # Duration of mild infections
FracSevere = 15 / 100 # Fraction of infections that are severe
FracCritical = 5  / 100 # Fraction of infections that are critical
FracMild = 1 - FracSevere - FracCritical # Fraction of infections that are mild
ProbDeath = 0.4  # Probability of dying given critical infection
CFR = ProbDeath * FracCritical # Case fatality rate (fraction of infections resulting in death)
TimeICUDeath = 8 # Time from ICU admission to death, days
DurHosp = 6 # Duration of hospitalization (severe infections), days
N.start = 1000 # population size
starting.probs = c(0.98,0.01,0.01) # Starting probabilities to seed outbreak in states Susceptible, Exposed and Mild
b1 = 0.33 # Transmission rate (mild infections)
b21 = 0.1 # Transmission rate (severe infections, relative to mild)
b31 = 0.1 # Transmission rate (critical infections, relative to mild)
max.day = 200 # maximum day to run simulation to

# quick check of death rate
u = (1/TimeICUDeath)*(CFR/FracCritical)

# Set simulation period
maxAge <- 105 # maximum age
sex <- c("m",'f') # 
first.date = "01/01/2020"
# set initial simulation horizon
days = 0
start.date = format(as.Date(first.date, format='%d/%m/%Y') + days, '%d/%m/%Y') # add days, then convert back to character
end.date = format(as.Date(first.date, format='%d/%m/%Y') + days + 1, '%d/%m/%Y') # add days, then convert back to character
simHorizon <- setSimHorizon(startDate=start.date, endDate=end.date)

# Define the health states through which individuals transition
health <- c("S", "E", "I1", "I2", "I3", "R", "D")
stateSpace <- expand.grid(sex=sex, health=health) # 
absStates <- c("dead") # absorbing states; no migration

## Definition of an initial population 
# Use a randomly-generated population
# a) random uniform
initBirthDatesRange <- chron(dates=c("31/12/1930", "31/12/2007"), format=c(dates="d/m/Y"), 
                             out.format=c(dates="d/m/year"))
birthDates <- dates(initBirthDatesRange[1] + runif(N.start, min=0, max=diff(initBirthDatesRange)))
# b) matching the Australian adult population
load('data/AusPop.RData')
sampled.ages = base::sample(ages$age, replace=TRUE, size=N.start, prob = ages$p)
birthDates <- dates(as.numeric(as.Date(start.date,'%d/%m/%Y')) - 365.25*sampled.ages)

# create initial population
initPop <- data.frame(ID=1:N.start, birthDate=birthDates, initState=sapply(birthDates, getRandInitState), stringsAsFactors = FALSE)
# convert character format of birthDate to chron
initPop <- transform(initPop, birthDate = dates(chron(birthDate, format = c(dates = "d/month/Y"), out.format = c(dates="d/m/year"))))

## Transition pattern and assignment of functions specifying transition rates
sicktrMatrix <- cbind(c("m/S->m/E", "m/E->m/I1", "m/I1->m/I2", "m/I2->m/I3", "m/I1->m/R", "m/I2->m/R", "m/I3->m/R", "m/I3->m/D",
                        "f/S->f/E", "f/E->f/I1", "f/I1->f/I2", "f/I2->f/I3", "f/I1->f/R", "f/I2->f/R", "f/I3->f/R", "f/I3->f/D"), 
                      c("S.2.E",    "E.2.I1",    "I1.2.I2",    "I2.2.I3", "I1.2.R", "I2.2.R", "I3.2.R","I3.2.D",
                        "S.2.E",    "E.2.I1",    "I1.2.I2",    "I2.2.I3", "I1.2.R", "I2.2.R", "I3.2.R","I3.2.D"))
allTransitions <- rbind(sicktrMatrix)

# absorbing transitions (long-term deaths)
absTransitions <- c("dead", "mortRates") # state and function
# overall transition matrix:
transitionMatrix <- buildTransitionMatrix(allTransitions=allTransitions, absTransitions=absTransitions, stateSpace=stateSpace)

# loop through days
transitions = NULL
for (days in 0:max.day){
  if(days%%10==0){cat('Day =', days, '.\n', sep='')} # 
  
  start.date = format(as.Date(first.date, format='%d/%m/%Y') + days, '%d/%m/%Y') # add days, then convert back to character
  end.date = format(as.Date(first.date, format='%d/%m/%Y') + days + 1, '%d/%m/%Y') # add days, then convert back to character
  simHorizon <- setSimHorizon(startDate=start.date, endDate=end.date)
  
  # dynamically update initial state
  if(days > 0){
    initPop = group_by(pop, ID) %>%
      arrange(ID) %>%
      slice(n()) %>% # latest transition if two per day
      mutate(
        update = ifelse(is.na(To)==TRUE, initState, To)) %>% # leave as before if not moved
        ungroup() %>%
        select(-transitionTime, -transitionAge, -From, -To, -initState) %>%
      rename('initState' = 'update') # rename
      #filter(!initState %in% c('dead')) # remove deaths 
      #      filter(!initState %in% c('dead','f/D','f/R','m/D','m/R')) # remove deaths and recovered (no longer need to estimate transitions)
    initPop = data.frame(initPop) # needs to be a data frame
  }
  
  # break out of loop if population is depleted
  N = nrow(initPop) # update population size (account for losses to death and recovery)
  cat('N=', N,'\n')
  if(N==1){break}
  
  # update daily I1 to I3 numbers
  I1 = length(grep('I1',initPop$initState))
  I2 = length(grep('I2',initPop$initState))
  I3 = length(grep('I3',initPop$initState))
  I.numbers = c(I1, I2, I3)

# Run microsimulation 
  pop <- micSim(initPop=initPop, immigrPop=NULL, transitionMatrix=transitionMatrix, 
                      absStates=absStates, initStates=NULL, initStatesProb=NULL, 
                      maxAge=maxAge, simHorizon=simHorizon)
  # concatenate transitions
  to.add = filter(pop, !is.na(From)) %>%
    mutate(day = days)
  transitions = bind_rows(transitions, to.add)
}

# save the results and meta-data on the parameters
meta = list("IncubPeriod"=IncubPeriod,"DurMildInf"=DurMildInf,"FracSevere"=FracSevere,
                  "FracCritical"=FracCritical,"ProbDeath"=ProbDeath,"DurHosp"=DurHosp,
                  "TimeICUDeath"=TimeICUDeath,"N.start"=N.start, 'first.date'=first.date, 
                  'b1'=b1,'b21'=b21,'b31'=b31,
                  'max.day'=max.day, starting.probs=starting.probs)
save(transitions, meta, file=outfile)
