# functions.micsim.rates.R
# rate functions for micsim
# March 2020

# Definition of (possible) transition rates (that vary along age, calendar time, 
# and duration). BEWARE: Each function that describes transition rates has to feature
# at least age (in years) as input parameter. Additionally, calTime (in years) and duration (in years) 
# might enter the function.

# susceptible to exposed
S.2.E <- function(age, calTime, duration){
  I1 = I.numbers[1]
  I2 = I.numbers[2]
  I3 = I.numbers[3]
  B1 = (b1/N) * 365 # Transmission rate (mild infections)
  b2 = b21 *B1  # Transmission rate (severe infections, relative to mild)
  b3 = b31 *B1  # Transmission rate (critical infections, relative to mild)
  rate = B1*I1 + b2*I2 + b3*I3
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
  p2 = (1/DurHosp)*(FracCritical/(FracSevere+FracCritical))
  rate = p2 * 365 # 'p2' in Alison's model
  return(rate)
}

# I1 to R (mild to recovered)
I1.2.R <- function(age, calTime, duration){
  g1 = (1/DurMildInf)*FracMild
  rate = g1 * 365 # 'g1' in Alison's model
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
  if(FracCritical==0){
    u=0
  }else{
    u=(1/TimeICUDeath)*(CFR/FracCritical)
  }
  g3 = (1/TimeICUDeath) - u
  rate = g3 * 365 # 'g3' in Alison's model
  return(rate)
}

# I3 to D
I3.2.D <- function(age, calTime, duration=0){
  if(FracCritical==0){
    u=0
    }else{
    u = (1/TimeICUDeath)*(CFR/FracCritical)
    }
  rate <- u * 365 # u (mu) in Alison's model
  return(rate)
}

## long-term annual mortality rates (function of age, same for men and women)
mortRates <- function(age, calTime, duration=0){
  a <- 0.00003
  b <- 0.097
  rate <- a*exp(b*age)
  rate = 0 # for comparing with Alison
  return(rate)
}
