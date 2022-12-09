

#######################------------   BEAKED REDFISH Parameters -------------####################

Blim <- 227000  #

recage <- 2
amax <- 40    
age <- recage:amax
ac<-length(age)

fbar <- c(17,38)  # in assessment it is Ages 19+ ratio 
mbar <- c(17, 38)



averageF <- 0.022 #from assessment 

tmax<- 250   # simulation time


### Natural mortality ###
Nmort<-rep(.05,ac)


####### Ricker ######

start.rec <- 218296300  #millions

reccap <- 796500000 #max(redfishass2$Rec.age2thousand)*1000 *1.5

alphar <- 5.712472       #asessment 2018
betar <- -6.318309e-07    #2018


ricker.sd<- 0.8874451 #0.8004844 

AR1par<- 0.6168852   #0.8196807 


### BH ###
#alphabh <- 0.001629795  #2019 assessment
#betabh <- 7.481104e-09  #2019 assessment

alphabh <- 0.0005289017 #2018 assessment
betabh <- 9.166759e-09 #2018 assessment 

### maturity ###
a50 <- 12.31206
env <- 2.180238

maturity <- round(1/(1+exp(-((age-a50)/env))),2) 



### weights ###

k <- 0.05147513 
b <- 2.018299 
Winf <- 1.937196        


weights <- Winf * (1-exp(-k * age))^b 



### selectivity ####
s50 <- 10.09438
ss <- 1.268491


Fsel <- 1/(1+exp(-((age-s50)/ss)))






