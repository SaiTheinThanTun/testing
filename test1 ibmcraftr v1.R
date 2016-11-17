#test1: transition between birth and death
#parameters and initial values for ibmcraftr

library(ibmcraftr)

init.pop <- c(S=1000, E=0, I_S=1, R_T=0, S_I=0, I_UA=0, I_DA=0, D = 1000)
pop <- syn_pop(init.pop)

mu_i = 1/(55*364) #birth rate of human
mu_S = 1/(55*364) #normal death rate of human

pmu <- rate2prob(mu_i)
StoD <- cumprob(pmu)
#?run_state_trans

timesteps <- 10000
param <- list(
  list(1,8,pmu),
  list(8,1,pmu)
)

result <- run_state_trans(timesteps,param,pop)


tail(result) #summary results at the end of the timesteps