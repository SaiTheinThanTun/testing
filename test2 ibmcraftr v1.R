#test2: transition between S and E
#parameters and initial values for ibmcraftr

library(ibmcraftr)

rm(list = ls())
#humans
#init
init.pop <- c(S=1000, E=0, I_S=1, R_T=0, S_I=0, I_UA=0, I_DA=0, D = 1000)
pop <- syn_pop(init.pop)
#parameters

mu_i = 1/(55*364) #birth rate of human
mu_S = 1/(55*364) #normal death rate of human

pmu <- rate2prob(mu_i) #birth/death probability of human

#transients
state_sumt <- "state_sum <- colSums(pop)"
Nt <- "N <- sum(state_sum[-length(state_sum)])" #exclude Deaths
Mt <- "M <- S_M + I_M"
mt <- "m <- M/N" #ratio of mosquito to human
xt <- "x <- (state_sum[3]+state_sum[6]+state_sum[7])/N" #3,6,7 -> I_S,I_UA,I_DA
zt <- "z <- I_M/M"
lam_Mt <- "lam_M <- a*c*x"
lam_Ht <- "lam_H <- m*a*b*z"

#transient ODE
S_Mt <- "S_M <- "

#mosquitos
#init
S_M=500
I_M=200
#parameters
beta = 1/10 #birth rate of mosquitos
kappa_SM = 1/10 #death rate of mosquitos
kappa_IM = 1/10 #in search of equilibrium :) 1/8, #death rate of infected mosquitos

a = 3 #human feeding rate
b = .006 #.01, #probability of disease transmission per bite for humans
c = .02 #.03 #probability a mosquito becomes infected after biting an infected human

#?run_state_trans parameters

timesteps <- 10000
param <- list(
  list(1,8,pmu),
  list(8,1,pmu)
)
transient <- c(Mt)

result <- run_state_trans(timesteps,param,transient=transient,pop)


tail(result) #summary results at the end of the timesteps