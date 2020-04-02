# 1_microsim_covid_vary_hospital.R
# run the COVID model using microsimulation based on a model using ODEs
# version with varying parameters
# version that looks at the impact on the hospital
# lyra version
# March 2020
#library(MicSim) # using my slightly modified versions
library(chron)
library(snowfall)
library(snow)
library(rlecuyer)
#library(stringr)
#library(dplyr) # not for lyra
source('99_functions.micsim.R') # functions to start micsim
source('99_functions.micsim.rates.R') # functions for transition rates
source('99_micSim.R') # adapted from MicSim by me (Vectorize)
source('99_getNextStep.R') # adapted from MicSim by me (Vectorize)
source('https://github.com/cran/MicSim/raw/master/R/auxFctMicSim.r') # from github repo of MicSim

# save results here
if(exists('outfile')==FALSE){outfile = 'data/simresults.RData'}

# key parameters
IncubPeriod = 5 # Incubation period, days
IncubPeriod = runif(n=1, min=0.9*IncubPeriod, max=1.1*IncubPeriod)
DurMildInf = 6   # Duration of mild infections
DurMildInf = runif(n=1, min=0.9*DurMildInf, max=1.1*DurMildInf)
FracSevere = 15 / 100 # Fraction of infections that are severe
FracSevere = runif(n=1, min=0.9*FracSevere, max=1.1*FracSevere)
FracCritical = 5  / 100 # Fraction of infections that are critical
FracCritical = runif(n=1, min=0.9*FracCritical, max=1.1*FracCritical)
FracMild = 1 - FracSevere - FracCritical # Fraction of infections that are mild
ProbDeath = 0.4  # Probability of dying given critical infection
ProbDeath = runif(n=1, min=0.9*ProbDeath, max=1.1*ProbDeath)
CFR = ProbDeath * FracCritical # Case fatality rate (fraction of infections resulting in death)
TimeICUDeath = 8 # Time from ICU admission to death, days
TimeICUDeath = runif(n=1, min=0.9*TimeICUDeath, max=1.1*TimeICUDeath)
DurHosp = 6 # Duration of hospitalization (severe infections), days
DurHosp = runif(n=1, min=0.9*DurHosp, max=1.1*DurHosp)
N.start = 1000 # population size
starting.probs = c(0.99,0.01,0.00) # Starting probabilities to seed outbreak in states Susceptible, Exposed and Mild
b1 = 0.33 # Transmission rate (mild infections)
b1 = runif(n=1, min=0.9*b1, max=1.1*b1)
b21 = 0.1 # Transmission rate (severe infections, relative to mild)
b21 = runif(n=1, min=0.9*b21, max=1.1*b21)
b31 = 0.1 # Transmission rate (critical infections, relative to mild)
b31 = runif(n=1, min=0.9*b31, max=1.1*b31)
max.day = 200 # maximum day to run simulation to
# hospital parameters:
qld.pop = 5115451
presentHOther = 180830

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
health <- c("S", "HS", "E", "HE", "I1", "H1", "I2", "I3", "R", "HR", "D")
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
initPop$initState[1:5] = 'm/E' # temporary, first 5 exposed
# convert character format of birthDate to chron
initPop <- transform(initPop, birthDate = dates(chron(birthDate, format = c(dates = "d/month/Y"), out.format = c(dates="d/m/year"))))

## Transition pattern and assignment of functions specifying transition rates
sicktrMatrix <- cbind(c("m/S->m/E", "m/S->m/HS", "m/HS->m/S", "m/R->m/HR", "m/HR->m/R", "m/E->m/HE", "m/HE->m/E", "m/I1->m/H1", "m/H1->m/I1", "m/E->m/I1", "m/I1->m/I2", "m/I2->m/I3", "m/I1->m/R", "m/I2->m/R", "m/I3->m/R", "m/I3->m/D",
                        "f/S->f/E", "f/s->f/HS", "f/HS->f/S", "f/R->f/HR", "f/HR->f/R", "f/E->f/HE", "f/HE->f/E", "f/I1->f/H1", "f/H1->f/I1", "f/E->f/I1", "f/I1->f/I2", "f/I2->f/I3", "f/I1->f/R", "f/I2->f/R", "f/I3->f/R", "f/I3->f/D"), 
                      c("S.2.E",    "S.2.HS",    "HS.2.S",    "R.2.HR",    "HR.2.R",    "E.2.HE",    "HE.2.E",    "I1.2.H1",    "H1.2.I1",    "E.2.I1",    "I1.2.I2",    "I2.2.I3", "I1.2.R", "I2.2.R", "I3.2.R","I3.2.D",
                        "S.2.E",    "S.2.HS",    "HS.2.S",    "R.2.HR",    "HR.2.R",    "E.2.HE",    "HE.2.E",    "I1.2.H1",    "H1.2.I1",    "E.2.I1",    "I1.2.I2",    "I2.2.I3", "I1.2.R", "I2.2.R", "I3.2.R","I3.2.D"))
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
  
  # dynamically update initial state (different in lyra)
  if(days > 0){
    #
    initPopx = by(pop, pop$ID, tail, n=1) # latest transition if two per day 
    initPop = NULL # convert list to frame
    for (j in 1:length(initPopx)){
      initPop = rbind(initPop, initPopx[[j]])
    }
    initPop$update = initPop$initState # start with initial state
    initPop$update[is.na(initPop$To)==FALSE] = initPop$To[is.na(initPop$To)==FALSE] # update if changed
    initPop = initPop[, c('ID','birthDate','update')]
    names(initPop)[3] = 'initState'
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
  # concatenate transitions (different in lyra)
  index = is.na(pop$From) == FALSE
  if(sum(index)>0){
    to.add = pop[index,]
    to.add$day = days
    transitions = rbind(transitions, to.add)
  }
}

# save the results and meta-data on the parameters
meta = list("IncubPeriod"=IncubPeriod,"DurMildInf"=DurMildInf,"FracSevere"=FracSevere,
                  "FracCritical"=FracCritical,"ProbDeath"=ProbDeath,"DurHosp"=DurHosp,
                  "TimeICUDeath"=TimeICUDeath,"N.start"=N.start, 'first.date'=first.date, 
            'b1'=b1,'b21'=b21,'b31'=b31,
                  'max.day'=max.day, 'starting.probs'=starting.probs)
save(transitions, meta, file=outfile)
