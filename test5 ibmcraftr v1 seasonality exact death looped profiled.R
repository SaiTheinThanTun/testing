#test5: test5 seasonality , add exact deaths in all compartments, looped profiled
#parameters and initial values for ibmcraftr
library(ibmcraftr)
#library(profvis)

system.time({
no_sims <- 100
lci <- .05
hci <- .95
#####################################
#######################CAUTION!!!#### It will take a large amount of time!!!!!!!!!!!!!!!
##################################### You can skip to "loading RData" section w/o running the loop

#rm(list = ls()) #####remove this when looping to run many times####



#creating a list of simulations
sims <- list()
for(i in 1:no_sims){
#humans
#init
init.pop <- c(S=1000, E=0, I_S=1, R_T=0, S_I=0, I_UA=0, I_DA=0, D = 1000)
pop <- syn_pop(init.pop)
#parameters

mu_i = 1/(55*364) #birth rate of human
mu_S = 1/(55*364) #normal death rate of human

pmu <- rate2prob(mu_i) #birth/death probability of human

omega_0 = 1/(5*364) #total loss of immunity (5 years), #then the model will need to run over 5 years
mu_E = 1 #same death rate as susceptibles
f = 1/10 #duration of E is approximately 10 days, Collins
mu_IS = 3 # increased death rate due to untreated malaria, 3 times higher chance
#gamma = 1/(10), # recovery with treatment
#gamma = 1/(20), # natural recovery, approximately from Collins
gamma = 1/(20) #recovery rate for both natural recovery and treatment
rep = .70 # treatment+report coverage
mu_RT = 1.0001 # increased death rate from treatment, 1.0001 times
sigma = 1/(6*30) # loss of prophylactic effect of drugs/treatment
mu_SI = 1 # same death rate as mu_S
xi1 = .3 #30% of reinfected individuals become sympatomatic
xi2 = .6 #60% of reinfected & Asymptomatic individuals become detectable with some Dx tests
#lam_SI = .60, #FOI to become truely sympatomatic infection is reduced by 40%
#lam_IUA = .85, #FOI of Undectable Asymptomatic infection becoming symptomatic again is decreased by 15% 
#lam_IDA = .75, #FOI of Dectable Asymptomatic infection becoming symptomatic again is decreased by 25%
omega_IUA = 1/365 # loss of all parasites at one year after becoming subpatent
omega_IDA = 1/100 # loss of detectability of parasites
mu_IUA = 1 #same death rate as mu_S
mu_IDA = 1 #same death rate as mu_S
timecounter = 0

#seasonality
seas_switch <- 1 #logical switch for seasonality
amp <- .2
phi <- 0
magnitude <- 1

#transients #this should be done outside of the function once
state_sumt <- "state_sum <- colSums(pop)"
Nt <- "N <- sum(state_sum[-length(state_sum)])" #exclude Deaths
Mt <- "M <- S_M + I_M"
mt <- "m <- M/N" #ratio of mosquito to human
xt <- "x <- (state_sum[3]+state_sum[6]+state_sum[7])/N" #3,6,7 -> I_S,I_UA,I_DA
zt <- "z <- I_M/M"
seast <- "seas <- amp*cos(2*pi*(timecounter-phi)/365)+magnitude"
lam_Mt <- "lam_M <- a*c_*x*seas"
timecountert <- "timecounter <- timecounter + 1"
lam_Ht <- "lam_H <-  rate2prob(m*a*b*z)"
lam_Ht_dyanmic <- "lam_H <- param[[1]][[3]] <- rate2prob(m*a*b*z)"

transient_para <- c(state_sumt,Nt,Mt,mt,xt,zt,seast,lam_Mt,timecountert)

#transient ODE
S_Mt <- "S_M <- S_M + M*beta_ - S_M*kappa_SM - S_M*lam_M"
I_Mt <- "I_M <- I_M + S_M*lam_M - I_M*kappa_IM"

transient_ode <- c(S_Mt,I_Mt)

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


timesteps <- 10000
#transients

eval(parse(text = c(transient_para,lam_Ht)))

param <- list(
  list(1, c(2,8), c(lam_H,pmu)),
  list(2, c(3,8), c(rate2prob(f),pmu*mu_E)),
  list(3, c(4,7,8), c(rep*rate2prob(gamma),(1-rep)*rate2prob(gamma)),pmu*mu_IS),
  list(4, c(5,8),c(rate2prob(sigma),pmu*mu_RT)),
  list(5, c(1,3,6,7,8), c(rate2prob(omega_0),lam_H*xi1,lam_H*(1-xi2)*(1-xi1),lam_H*xi2*(1-xi1)),pmu*mu_SI),
  list(6, c(3,5,7,8), c(lam_H*xi1,rate2prob(omega_IUA),lam_H*(1-xi1)),pmu*mu_IUA),
  list(7, c(3,6,8), c(lam_H*xi1,rate2prob(omega_IDA)),pmu*mu_IDA),
  list(8,1,pmu)
)


transient_all <- c(transient_para, lam_Ht_dyanmic, transient_ode)
result <- run_state_trans(timesteps,param,transient=transient_all,pop)
    sims[[i]] <- run_state_trans(timesteps,param,transient=transient_all,pop)
  }
})

#for single run
# plot(result[,1]+result[,5], ylim=c(0, N), type='l', col='blue')
# lines(result[,3], type='l', col='red')
# lines(result[,6]+result[,7], type='l', col='purple')

##############################
############end###############
##############################

#saving the dataset
save(sims,file=paste("100iterations_",Sys.Date(),".RData",sep=""))

#plotting

tmp_avg <- rep(NA,no_sims)
avg_sims <- matrix(NA,nrow(sims[[1]]),ncol(sims[[1]])) #initializing a blank dataset of summary table
for(i in 1:ncol(sims[[1]])){ #outer loop for the columns
  for(j in 1:nrow(sims[[1]])){ #inner loop for the rows
    for(k in 1:no_sims){#innermost loop for no. of simulations(3rd dimension)
      tmp_avg[k] <- sims[[k]][j,i]
    }
    avg_sims[j,i] <- mean(tmp_avg)
  }
}

#lower CI (LCI)
tmp_lci <- rep(NA,no_sims)
lci_sims <- matrix(NA,nrow(sims[[1]]),ncol(sims[[1]])) #initializing a blank dataset of summary table
for(i in 1:ncol(sims[[1]])){ #outer loop for the columns
  for(j in 1:nrow(sims[[1]])){ #inner loop for the rows
    for(k in 1:no_sims){#innermost loop for no. of simulations(3rd dimension)
      tmp_lci[k] <- sims[[k]][j,i]
    }
    lci_sims[j,i] <- quantile(tmp_lci, probs= lci, na.rm=TRUE)
  }
}

#high CI (HCI)
tmp_hci <- rep(NA,no_sims)
hci_sims <- matrix(NA,nrow(sims[[1]]),ncol(sims[[1]])) #initializing a blank dataset of summary table
for(i in 1:ncol(sims[[1]])){ #outer loop for the columns
  for(j in 1:nrow(sims[[1]])){ #inner loop for the rows
    for(k in 1:no_sims){#innermost loop for no. of simulations(3rd dimension)
      tmp_hci[k] <- sims[[k]][j,i]
    }
    hci_sims[j,i] <- quantile(tmp_hci, probs = hci, na.rm=TRUE)
  }
}

colnames(hci_sims) <- colnames(lci_sims) <- colnames(avg_sims) <- c('S','E','I_S', 'R_T','S_I','I_UA','I_DA','D') #column names for the summary table
ts <- 1:length(avg_sims[,1])
avg_sims <- cbind(avg_sims,ts)
hci_sims <- cbind(hci_sims,ts)
lci_sims <- cbind(lci_sims,ts)

#saving averages
save(avg_sims,hci_sims,lci_sims, file=paste("100it_avg_hci_lci_",Sys.Date(),".RData",sep=""))



# ###loading Data####
# load("D:\\Dropbox\\IBM project_Sai\\testingRData\\100it_avg_hci_lci.RData")
# load("D:\\Dropbox\\IBM project_Sai\\testingRData\\100iterations.RData")

####plot S vs I ####
par(mar=c(5,4,4,4))
plot(avg_sims[,9],avg_sims[,1]+avg_sims[,5], type="l", col="blue", axes=FALSE, xlab="", ylab="", main=paste(no_sims," iterations: S vs I"))# lambda",lam_h,"and CI",lci,'-',hci))
polygon(c(avg_sims[,9], rev(avg_sims[,9])),c(hci_sims[,1]+hci_sims[,5], rev(lci_sims[,1]+lci_sims[,5])),col=rgb(0,0,100,50,maxColorValue=255), border=NA)
axis(2, ylim=c(0,17),col="blue") 
mtext("Susceptible humans",side=2,line=2.5) 

box()
par(new=TRUE)
plot(avg_sims[,9],avg_sims[,3], type="l", col="red", axes=FALSE, xlab="", ylab="")
polygon(c(avg_sims[,9], rev(avg_sims[,9])), c(hci_sims[,3], rev(lci_sims[,3])),col=rgb(100,0,0,50,maxColorValue = 255), border=NA)
axis(4, ylim=c(0,17),col="red") 
mtext("Infected humans",side=4, line=2.5)

axis(1,pretty(range(avg_sims[,9]),10))
mtext("Time (days)",side=1,col="black",line=2.5)

legend("top",legend=c("Susceptibles","Infected"),
       text.col=c("blue","red"),pch= "__", col=c("blue","red"))


####plot S vs asympt ####
par(mar=c(5,4,4,4))
plot(avg_sims[,9],avg_sims[,1]+avg_sims[,5], type="l", col="blue", axes=FALSE, xlab="", ylab="", main=paste(no_sims," iterations: S vs Asymptomatics"))# lambda",lam_h,"and CI",lci,'-',hci))
polygon(c(avg_sims[,9], rev(avg_sims[,9])),c(hci_sims[,1]+hci_sims[,5], rev(lci_sims[,1]+lci_sims[,5])),col=rgb(0,0,100,50,maxColorValue=255), border=NA)
axis(2, ylim=c(0,17),col="blue") 
mtext("Susceptible humans",side=2,line=2.5) 
box()
par(new=TRUE)
plot(avg_sims[,9],avg_sims[,6]+avg_sims[,7], type="l", col="purple", axes=FALSE, xlab="", ylab="")
polygon(c(avg_sims[,9], rev(avg_sims[,9])), c(hci_sims[,6]+hci_sims[,7], rev(lci_sims[,6]+lci_sims[,7])),col=rgb(100,0,0,50,maxColorValue = 255), border=NA)
axis(4, ylim=c(0,17),col="purple") 
mtext("Asymptomatic Infected humans",side=4, line=2.5)

axis(1,pretty(range(avg_sims[,9]),10))
mtext("Time (days)",side=1,col="black",line=2.5)

legend("top",legend=c("Susceptibles","Asymptomatic Infected"),
       text.col=c("blue","purple"),pch= "__", col=c("blue","purple"))


#combination plot
plot(avg_sims[,9],avg_sims[,1]+avg_sims[,5], type="l", col="blue", ylim=c(0,max(avg_sims[,1])), xlab="Days", ylab="No. individuals", main=paste(no_sims," iterations"))# lambda",lam_h,"and CI",lci,'-',hci))
polygon(c(avg_sims[,9], rev(avg_sims[,9])),c(hci_sims[,1]+hci_sims[,5], rev(lci_sims[,1]+lci_sims[,5])),col=rgb(0,0,100,50,maxColorValue=255), border=NA)


lines(avg_sims[,9],avg_sims[,3], type="l", col="red")
polygon(c(avg_sims[,9], rev(avg_sims[,9])), c(hci_sims[,3], rev(lci_sims[,3])),col=rgb(100,0,0,50,maxColorValue = 255), border=NA)


lines(avg_sims[,9],avg_sims[,6]+avg_sims[,7], type="l", col="purple")
polygon(c(avg_sims[,9], rev(avg_sims[,9])), c(hci_sims[,6]+hci_sims[,7], rev(lci_sims[,6]+lci_sims[,7])),col=rgb(50,0,0,50,maxColorValue = 255), border=NA)

legend("top",legend=c("Susceptibles","Infected","Asymptomatic Infected"),
       text.col=c("blue","red","purple"),pch= "__", col=c("blue", "red","purple"))

#combination plot compared with ODE

library(deSolve)
library(maemod)
maxtime <- 10000

out <- maemodrun("D:\\Dropbox\\IBM project_Sai\\ODE\\SIRSI.txt", timegrid = seq(0,maxtime,1)) #scenario2

plot(avg_sims[,9],avg_sims[,1]+avg_sims[,5], type="l", col="blue", ylim=c(0,max(avg_sims[,1])), xlab="Days", ylab="No. individuals", main=paste(no_sims," iterations"))# lambda",lam_h,"and CI",lci,'-',hci))
polygon(c(avg_sims[,9], rev(avg_sims[,9])),c(hci_sims[,1]+hci_sims[,5], rev(lci_sims[,1]+lci_sims[,5])),col=rgb(0,0,100,50,maxColorValue=255), border=NA)


lines(avg_sims[,9],avg_sims[,3], type="l", col="red")
polygon(c(avg_sims[,9], rev(avg_sims[,9])), c(hci_sims[,3], rev(lci_sims[,3])),col=rgb(100,0,0,50,maxColorValue = 255), border=NA)


lines(avg_sims[,9],avg_sims[,6]+avg_sims[,7], type="l", col="purple")
polygon(c(avg_sims[,9], rev(avg_sims[,9])), c(hci_sims[,6]+hci_sims[,7], rev(lci_sims[,6]+lci_sims[,7])),col=rgb(50,0,0,50,maxColorValue = 255), border=NA)


#humans
lines(out[,2]+out[,6], type='l', ylim=c(0,max(c(out[,2],out[,4]))), col="blue", main="Humans")
lines(out[,4], col="red")
lines(out[,7]+out[,8], col="purple")

legend("top",legend=c("Susceptibles","Infected","Asymptomatic Infected"),
       text.col=c("blue","red","purple"),pch= "__", col=c("blue", "red","purple"))