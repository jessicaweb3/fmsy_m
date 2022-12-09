

########################------------   NEA COD Parameters -------------#######################

#ICES asssessment based data 

Blim<- 220000 ### Blim tonnes (official advice sheet)

recage <- 3
amax <- 15     ## 13 age classes 3-15  
age <- recage:amax
ac<-length(age)
fbar<-c(2,7)    # in reality age 5-10
mbar <- c(2,7)


averageF <- 0.6 # ICES assessment 

tmax<- 250   # simulation time

### Natural mortality ###
Nmort<-rep(.2,ac)

####### Ricker ######


start.rec <-  748681600  # 

reccap <-  3866920500  # max(codass$Rec.age3)*1.5

alphar <- 7.953965 

betar <- -1.176552e-06

ricker.sd <- 0.6447752

AR1par <- 0.5190351


###### B-H ######

alphabh <- 0.0004234719  #0.4234719 

betabh <- 8.498407e-10 #8.498407e-07

#bh.sd <- 646770.4

#AR1par <- 0.5080555

### maturity ###

a50 <- 7.919142
env <- 1.184689

maturity <- round(1/(1+exp(-((age-a50)/env))),2)


### weights ###
k <- 0.1513596
b <- 5.127882
Winf <- 28.04999 


weights <- Winf * (1-exp(-k * age))^b 


### f sel ###
s50 <- 5.811967
ss <- 1.307185


Fsel <- 1/(1+exp(-((age-s50)/ss)))

#Fsel <- c(0.5, 0.5, 0.5, rep(1, ac-3))  

