#TEST6 with same parameters
library(ibmcraftr)
library(deSolve)


#parameters in the global environment
init.pop <- c(S=300, E=0, I_S=1, R_T=0, S_I=0, I_UA=0, I_DA=0, D = 300) 
maxtime <- 1000
mu_i = 1/(55*364) #birth rate of human
mu_S = 1/(55*364) #normal death rate of human

pmu <- rate2prob(mu_i) #birth/death probability of human, RUNS FASTER if defined here

omega_0 = 1/(5*364) #total loss of immunity (5 years), #then the model will need to run over 5 years
mu_E = 1 #same death rate as susceptibles
f = 1/10 #duration of E is approximately 10 days, Collins
mu_IS = 3 # increased death rate due to untreated malaria, 3 times higher chance
gamma = 1/(20) #recovery rate for both natural recovery and treatment
rep = .70 # treatment+report coverage
mu_RT = 1.0001 # increased death rate from treatment, 1.0001 times
sigma = 1/(6*30) # loss of prophylactic effect of drugs/treatment
mu_SI = 1 # same death rate as mu_S
xi1 = .3 #30% of reinfected individuals become sympatomatic
xi2 = .6 #60% of reinfected & Asymptomatic individuals become detectable with some Dx tests
omega_IUA = 1/365 # loss of all parasites at one year after becoming subpatent
omega_IDA = 1/100 # loss of detectability of parasites
mu_IUA = 1 #same death rate as mu_S
mu_IDA = 1 #same death rate as mu_S
timecounter = 0

#mosquitos
#init
S_M=500
I_M=200
#parameters
beta_ = 1/10 #birth rate of mosquitos
kappa_SM = 1/10 #death rate of mosquitos
kappa_IM = 1/10 #in search of equilibrium :) 1/8, #death rate of infected mosquitos

a = 3 #human feeding rate
b = .006 #.01, #probability of disease transmission per bite for humans
c_ = .02 #.03 #probability a mosquito becomes infected after biting an infected human

#seasonality
seas_switch <- 1 #logical switch for seasonality
amp <- .2 #.7 or 1 Lisa's advice
phi <- 0
magnitude <- 1

#parameterizing for ODE
times <- seq(0, maxtime, by = 1)
param_ode <- c(lambda = lambda)

#transient values in ODE
transient_ode <- c("", "", "")


#initialization for IBM
pop <- syn_pop(init.pop)

#transient values in IBM
#1.parameter values
transient_ibm_b4 <- c(
  "state_sum <- colSums(pop)",
  "N <- sum(state_sum[-length(state_sum)])",
  "M <- S_M + I_M",
  "m <- M/N",
  "x <- (state_sum[3]+state_sum[6]+state_sum[7])/N",
  "z <- I_M/M",
  "seas <- amp*cos(2*pi*(timecounter-phi)/365)+magnitude",
  "lam_M <- a*c_*x*seas",
  "timecounter <- timecounter + 1"
)
lam_Ht <- "lam_H <-  rate2prob(m*a*b*z)"
lam_Ht_dynamic <-  "lam_H <- param[[1]][[3]][[1]] <- rate2prob(m*a*b*z)"
#2.ODE functions inside IBM
transient_ibm_within <- c("S_M <- S_M + M*beta_ - S_M*kappa_SM - S_M*lam_M",
                          "I_M <- I_M + S_M*lam_M - I_M*kappa_IM")

#evaluation of the transient parameter values
eval(parse(text = c(transient_ibm_b4, lam_Ht)))
transient_all <- c(transient_ibm_b4, lam_Ht_dynamic, transient_ibm_within)

#model structures
#ODE
rateofchange <- c("dS <- -lambda*S",
                  "dI <- +lambda*S", "")

#IBM
param_ibm <- list(
  list(1, c(2,8), c(lam_H,pmu)),
  list(2, c(3,8), c(rate2prob(f),pmu*mu_E)),
  list(3, c(4,7,8), c(rep*rate2prob(gamma),(1-rep)*rate2prob(gamma),pmu*mu_IS)),
  list(4, c(5,8),c(rate2prob(sigma),pmu*mu_RT)),
  list(5, c(1,3,6,7,8), c(rate2prob(omega_0),lam_H*xi1,lam_H*(1-xi2)*(1-xi1),lam_H*xi2*(1-xi1),pmu*mu_SI)),
  list(6, c(3,5,7,8), c(lam_H*xi1,rate2prob(omega_IUA),lam_H*(1-xi1),pmu*mu_IUA)),
  list(7, c(3,6,8), c(lam_H*xi1,rate2prob(omega_IDA),pmu*mu_IDA)),
  list(8,1,pmu)
)


result <-  run_state_trans(timesteps = maxtime, param=param_ibm, transient = transient_all, pop)

#ODE function
testfun <- function(t, state, parameters)
{
  with(as.list(c(state, parameters)),
       {
         # define variables
         eval(parse(text = transient_ode))
         
         # rate of change
         eval(parse(text = rateofchange))
         
         # return the rate of change
         list(c(dS, dI))
       })
}

out <-  ode(y = init.pop, times = times, func = testfun,parms = param_ode)

#plot HERE
plot(result[, 1], type = 'l', col = 'blue', main = paste("IBM VS ODE"))
lines(result[, 2], type = 'l', col = 'red')
lines(out[, 2], type = 'l', col = 'blue')
lines(out[, 3], type = 'l', col = 'red')
