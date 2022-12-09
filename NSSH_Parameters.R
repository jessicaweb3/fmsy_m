
##########################------------   NSSH Parameters -------------#########################

Blim<- 2.5e6 ### Blim = 2.5m tonnes (official advice sheet)   

recage <- 2
amax <- 15     ## 14 age classes 2-15  
age <- recage:amax
ac<-length(age)

fbar<-c(3,10)    # in reality age 5-12 

mbar<-c(3,10) 

averageF <- 0.142 #ICES assessment 1988-2019

tmax<- 250   # simulation time

### Natural mortality ###
Nmort<-c(.9,rep(.15,(ac-1)))
Mvariation <- 0.015


##### Ricker #####

start.rec <-  16338375000   #billions

reccap <- 89781000000 #max(nsshass$Rec.age2)*1.5

alphar <- 9.296661

betar <- -2.963954e-07

# recruitment 5.1 - stochasticy 
ricker.sd <- 0.8294368
                      ## normally distributed error term (on log-scale)

# recruitment 5.2 - stochasticy + autocorrelation
AR1par <- 0.2618611          # autocorrelation bt years 


##### BH ####
alphabh <- 0.0001636856
betabh <- 7.495156e-11

bh.sd <- 9.075588e-11

##### Maturity #####
a50 <- 4.419048
env <- 0.4814544

maturity <- round(1/(1+exp(-((age-a50)/env))),2)



##### Weights #####
k <- 0.3482399
b <- 3.054745 
Winf <- 0.4025885 


weights <- Winf * (1-exp(-k * age))^b 


##### F selectivity #####
s50 <- 4.419048  #6.402239
ss <- 0.6673977   #1.610838 


Fsel <- 1/(1+exp(-((age-s50)/ss)))


#Fsel <- c(0, 0.5, rep(1,(ac-2))) 

