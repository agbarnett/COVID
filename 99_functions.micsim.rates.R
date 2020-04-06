# functions.micsim.rates.R
# rate functions for micsim
# April 2020

# Definition of (possible) transition rates (that vary along age, calendar time, 
# and duration). BEWARE: Each function that describes transition rates has to feature
# at least age (in years) as input parameter. Additionally, calTime (in years) and duration (in years) 
# might enter the function.

# susceptible to exposed
S.2.E <- function(age, calTime, duration){
  seas = 1 + seas.amp*cos(2*pi*(calTime - seas.phase)/365.25)
  I1 = I.numbers[1]
  I2 = I.numbers[2]
  I3 = I.numbers[3]
  B1 = (b1/N.pop) * 365 # Transmission rate (mild infections); use starting population size
  b2 = b21 *B1  # Transmission rate (severe infections, relative to mild)
  b3 = b31 *B1  # Transmission rate (critical infections, relative to mild)
  rate = (B1*I1 + b2*I2 + b3*I3) * seas
  return(rate)
}

# exposed to I1 (mild)
E.2.I1 <- function(age, calTime, duration){
  a = 1/ IncubPeriod;
  rate = a * 365 # 'a' in Alison's model
  return(rate)
}

# I1 to I2 (mild to severe)
I1.2.I2 <- function(age, calTime, duration){
  g1 = (1/DurMildInf)*FracMild
  p1 = (1/DurMildInf) - g1
  rate = p1 * 365 # 'p1' in Alison's model
  return(rate)
}

# I2 to I3 (sever to critical)
I2.2.I3 <- function(age, calTime, duration){
  age.curve = 1 # off for now
  p2 = (1/DurHosp)*(FracCritical/(FracSevere+FracCritical))
  rate = p2 * 365 * age.curve # 'p2' in Alison's model
  return(rate)
}

# I1 to R (mild to recovered)
I1.2.R <- function(age, calTime, duration){
  age.curve = 1 + ((age - 50)/10)*age.slope # higher rate of recovery in younger people
  g1 = (1/DurMildInf)*FracMild
  rate = g1 * 365 * age.curve # 'g1' in Alison's model
  return(rate)
}

# I1 to H1 (mild presenting to hospital)
I1.2.H1 <- function(age, calTime, duration){
  rate = -log(1-daily.prob.mild) # daily probability to rate
  return(rate)
}

# H1 to I1 (mild presenting at hospital back to I1)
H1.2.I1 <- function(age, calTime, duration){
  rate = ifelse(duration > 0.5/365.25, Inf, 0) # instant return after 1 day
  return(rate)
}

# I2 to R (severe to recovered)
I2.2.R <- function(age, calTime, duration){
  p2 = (1/DurHosp)*(FracCritical/(FracSevere+FracCritical))
  g2 = (1/DurHosp)-p2
  rate = g2 * 365 # 'g2' in Alison's model (gamma2)
  return(rate)
}

# I3 to R  (critical to recovered)
I3.2.R <- function(age, calTime, duration){
  age.curve = 1 + ((age - 50)/10)*age.slope # higher rate of recovery in younger people
  if(FracCritical==0){
    u=0
  }else{
    u=(1/TimeICUDeath)*(CFR/FracCritical)
  }
  g3 = (1/TimeICUDeath) - u
  rate = g3 * 365 * age.curve # 'g3' in Alison's model
  return(rate)
}

# I3 to D
I3.2.D <- function(age, calTime, duration){
  age.curve = 1 # keep at 1 for now
  if(FracCritical==0){
    u=0
    }else{
    u = (1/TimeICUDeath)*(CFR/FracCritical)
    }
  rate <- u * 365 * age.curve # u (mu) in Alison's model
  return(rate)
}

## exposed to hospital and back

# E to HE
E.2.HE <- function(age, calTime, duration){
  daily = -log(1-(presentHOther/qld.pop))/31 # based on Q health hospital presentation data
  rate = daily * 365
  return(rate)
}
# HE to E
HE.2.E <- function(age, calTime, duration){
  rate = ifelse(duration > 0.5/365.25, Inf, 0) # instant return after 1 day
  return(rate)
}

## susceptible to hospital and back

# S to HS
S.2.HS <- function(age, calTime, duration){
  daily = -log(1-(presentHOther/qld.pop))/31 # based on Q health hospital presentation data
  rate = daily * 365
  return(rate)
}

# HS to S
HS.2.S <- function(age, calTime, duration){
  rate = ifelse(duration > 0.5/365.25, Inf, 0) # instant return after 1 day
  return(rate)
}

## recovered to hospital and back

# R to HR
R.2.HR <- function(age, calTime, duration){
  daily = -log(1-(presentHOther/qld.pop))/31 # based on Q health hospital presentation data
  rate = daily * 365
  return(rate)
}
# HR to R
HR.2.R <- function(age, calTime, duration){
  rate = ifelse(duration > 0.5/365.25, Inf, 0) # instant return after 1 day
  return(rate)
}

## long-term annual mortality rates (function of age, same for men and women)
mortRates <- function(age, calTime, duration){
  a <- 0.00003
  b <- 0.097
  rate <- a*exp(b*age)
  return(rate)
}
