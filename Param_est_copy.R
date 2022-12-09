library(tidyverse)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(reshape2)
library(openxlsx)
library(TMB)
library(here)
library(png)

 
#####------------     NOTES     -------------------#####

# - nll for parameter estimation is log-likelihood 
# - SSB is real numbers & SSB1 is theoretical numbers  


#------------------------- DATASETS -------------------------------#

nsshass <-read.xlsx(here('Scripts_R/NSSH assessment.xlsx'),sheet=1)

nsshmaturity<-read.xlsx(here('Scripts/NSSH assessment.xlsx'),sheet=2)
nsshweight <- read.xlsx(here('Scripts/NSSH assessment.xlsx'), sheet=3)
nsshfmort<-read.xlsx(here('Scripts/NSSH assessment.xlsx'),sheet=4)


codass <-read.xlsx(here('Scripts_R/Cod_assessment.xlsx'),sheet=1)

#for biased plotting     codass <- codass[-74,]     
# colnames(codass)[5] <- "SSB"  
# colnames(codass)[6] <- "SSBhigh"      
# colnames(codass)[7] <- "SSBlow"

codmaturity<-read.xlsx(here('Scripts/Cod_assessment.xlsx'),sheet=2)
codweight <- read.xlsx(here('Scripts/Cod_assessment.xlsx'), sheet=3)
codfmort<-read.xlsx(here('Scripts/Cod_assessment.xlsx'),sheet=4)

#for biased plotting 
    #redfishass <-read.xlsx(here('Scripts/S.mentella_assessment.xlsx'),sheet=1) # Advice2020
    # colnames(redfishass)[3] <- "SSB"
#colnames(redfishass)[4] <- "SSBhigh"      
#colnames(redfishass)[5] <- "SSBlow"


redfishass <-read.xlsx(here('Scripts_R/S.mentella_assessment.xlsx'),sheet=2) # Advice2018

redfishass2 <-read.xlsx(here('Scripts_R/S.mentella_assessment.xlsx'),sheet=3)  # AFWG2019

redfishmaturity<-read.xlsx(here('Scripts_R/S.mentella_assessment.xlsx'),sheet=4)
redfishweight <- read.xlsx(here('Scripts_R/S.mentella_assessment.xlsx'), sheet=5)
redfishfmort<-read.xlsx(here('Scripts_R/S.mentella_assessment.xlsx'),sheet=6)


#######################################################################
######---------------      NSSH           --------------- #############
#######################################################################


#### Setting up general parameters 

age <- c(2:15)

recage<-2

start.rec <- mean(nsshass$Rec.age2)*1000

#### Recruitment: important to correct for the rec.age offset ####
tsl<-length(nsshass$Year)
Rec <- nsshass$Rec.age2[(recage+1):tsl]*1000 #assessment recruitment
SSB <- nsshass$SSB[1:(tsl-(recage))]         # assessment ssb

Rec <- Rec/1000000000  # for plotting
SSB <- SSB/1000000     #for plotting

## 5. Ricker 
data<-list()
data$SSB<-SSB #c(0,SSB)
data$Rec<-Rec #c(0,rep(mean(Rec),length(Rec)))

#nlsricker<- nls(Rec~alpha*SSB*exp(beta*SSB),data=data,start=list(alpha=1e7,beta=1e-7))
nlsrssb <- nls(log(Rec/SSB) ~ alpha+(beta*SSB), data=data,start=list(alpha=5,beta=-1))

alphar <- summary(nlsrssb)$coefficients[1]
betar <- summary(nlsrssb)$coefficients[2]

#Recruitment5 <- exp(alpha)*data$SSB*exp(beta*data$SSB)
RSSB5 <- alphar +betar *data$SSB
ricker.sd<-sd(residuals(nlsrssb))
ricker.mean<-mean(residuals(nlsrssb))
autocorrelation1<-acf(residuals(nlsrssb))
AR1par<-autocorrelation1$acf[2]

Rvar.std <- mean(exp(rnorm(1e6,0,ricker.sd))) 


SSB1 <- seq(0,9e6,1e5) # SSB1 theoretical numbers to test alpha & beta estimates 

SSB1 <- SSB1/1000000  #for plotting

nsshRecruitment5 <- exp(alphar)*SSB1*exp(betar *SSB1)  #- Deterministic stock-recruitment version

Rvariation<-rnorm(length(SSB1),0,ricker.sd)        #- Norm. dist error term (on log-scale)
Recruitment5.1 <- exp(alphar)*SSB1*exp(betar *SSB1) * exp(Rvariation)/Rvar.std #- Stochasticy term

Rvariationacf<-c(Rvariation[1],AR1par*Rvariation[1:(length(SSB1)-1)]+Rvariation[2:length(SSB1)]) # autocorrelation
Recruitment5.2 <- exp(alphar)*SSB1 *exp(betar *SSB1) *exp(Rvariationacf)/Rvar.std #- stochasticity term


plot(nsshRecruitment5~SSB1,type="l",col="red") #xlim=c(0, 9), ylim=c(0, 60), xlab="SSB1 million t",  ylab="Rec billions")   #,ylim=c(0,7e7))
points(Recruitment5.1~SSB1,col="blue",pch=1)
points(Recruitment5.2~SSB1,col="purple",pch=2)
points(Rec~SSB,pch=16)


##---- Compare with assessed SSB --- SSB real data 
Recruitment5 <- exp(alphar)*SSB*exp(betar*SSB)  #- Deterministic stock-recruitment version
Rvariation<-rnorm(length(SSB),0,ricker.sd)      #- Norm. dist. error term (on log-scale)
Recruitment5.1 <-  exp(alphar)*SSB *exp(betar *SSB) *exp(Rvariation)/Rvar.std #- Stochasticity 

Rvariationacf<-c(Rvariation[1],AR1par*Rvariation[1:(length(SSB)-1)]+Rvariation[2:length(SSB)]) #- AR1 autocorr. error term
Recruitment5.2 <-  exp(alphar)*SSB *exp(betar *SSB) *exp(Rvariationacf)/Rvar.std #- Stochasticity

plot(Recruitment5~SSB,type="l",col="red")
points(Recruitment5.1~SSB,col="blue",pch=1)
points(Recruitment5.2~SSB,col="purple",pch=2)
points(Rec~SSB,pch=16)

c(sum((Rec-Recruitment5)^2),sum((Rec-Recruitment5.1)^2),sum((Rec-Recruitment5.2)^2))

###### BH (trying again 2/2 2022)
# 1/R = beta + alpha * 1/SSB algebraic transformation of BH

tsl<-length(nsshass$Year)
Rec <- nsshass$Rec.age2[(recage+1):tsl]*1000 #assessment recruitment
SSB <- nsshass$SSB[1:(tsl-(recage))]         # assessment ssb

Rec <- Rec/1000000000  # for plotting
SSB <- SSB/1000000     #for plotting

data<-list()
data$SSB<-SSB #c(0,SSB)
data$Rec<-Rec #c(0,rep(mean(Rec),length(Rec)))

## BH algebraic transformation, estimating base of ICES data  
bhnlsrssb <- nls(1/Rec ~ beta + alpha * (1/SSB), data=data, start = list(alpha=5, beta=-1)) 

alphabh <- summary(bhnlsrssb)$coefficients[1]
betabh <- summary(bhnlsrssb)$coefficients[2]

bh.sd<-sd(residuals(bhnlsrssb))
bhvar.std <- mean(exp(rnorm(1e6,0, bh.sd))) 
bhvariation<- rnorm(length(SSB1),0,bh.sd) 

SSB1 <- seq(0, 9e6, 1e5) # SSB1 theoretical numbers to test alpha & beta estimates 
SSB1 <- SSB1/1000000  #for plotting

nsshRecruitment6 <- 1/(betabh + alphabh * 1/SSB1) # testing the alpha and beta, works!
Recruitment6.1 <- 1/(betabh + alphabh * 1/SSB1) * (1/(bhvariation))/bh.sd #- Stochasticy term 

plot(nsshRecruitment6 ~ SSB1,type="l",col="red")

#######   dont use any of this 
## BH - alpha and beta estimated in TMB - 2 options ##
#Recruitment6 <- exp(alphabh +log(SSB)-log(exp(betabh)*SSB))

#Recruitment6 <- (alphabh*SSB)/(1+betabh*SSB)   ### dont use this
#plot(Recruitment6 ~ SSB,type="l",col="red")


#----- Plotting all ---- #
rec.vector <- c(nsshRecruitment5, nsshRecruitment6)
type <-rep(c("Ricker5", "BH"),each=length(nsshRecruitment5))
rec.df <- data.frame(Type=type, SSB=rep(SSB1,2),Rec=rec.vector)

rec.df %>% ggplot(aes(x=SSB,y=Rec, color=Type)) + geom_line(size=1.5) +
  scale_color_brewer(palette="Accent") +
    theme_bw() + theme(panel.background = element_blank(), 
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"), axis.text.y = element_text(size=15), 
    axis.text.x = element_text(size=15), axis.title.x = element_text(size=20), 
    axis.title.y = element_text(size=20))


##### Maturity #####
  
maturity2 <-melt(nsshmaturity,id.vars="Year")
maturity2$Age <- rep(0:15,each=dim(nsshmaturity)[1])

maturity3 <- maturity2 %>% filter(Year>1987 & Age>1) 

maturity3 %>% 
  ggplot(aes(x=Age,y=value,color=as.factor(Year))) + geom_line()

compile('Scripts/maturity.cpp')
dyn.load(dynlib('Scripts/maturity'))

data<-list()
data$age<-maturity3$Age       # changed it from maturity2 to maturity3 to use correct filter
data$mprop<-maturity3$value

param <- list()
param$a50 <- 4
param$env <- .2
param$logsigma <-0

obj <- MakeADFun(data, param,DLL="maturity")   # MakeADFun - automatic differentiation function 
opt <- nlminb(obj$par, obj$fn, obj$gr)

## parameter estimates
a50 <- opt$par[1]
env <- opt$par[2]

maturity <- round(1/(1+exp(-((age-a50)/env))),2)   ## used in simulation 


maturity3$mprop.est <- 1/(1+exp(-((maturity3$Age-a50)/env)))

nsshmat <- maturity3 %>% 
  ggplot(aes(x=Age,y=value)) + geom_point(size= 2.5) +
  geom_line(inherit.aes=F,aes(x=Age,y=mprop.est), size= 0.6) + labs(y= "Maturity", x= "Age")+ theme_bw() + theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.y = element_text(size=20), axis.text.x = element_text(size=20), axis.title.x = element_text(size=30), axis.title.y = element_text(size=30)) + labs(tag="a") 




##### Weights #####

weights2 <-melt(nsshweight,id.vars="Year")
weights2$Age <- rep(0:15,each=dim(nsshweight)[1])

weights3 <- weights2 %>% filter(Year>1987 & Age>1) 

weights3 %>% 
  ggplot(aes(x=Age,y=value,color=as.factor(Year))) + geom_line()

compile('Scripts/weight.cpp')
dyn.load(dynlib('Scripts/weight'))

data<-list()
data$age<-weights3$Age         # changed it from weights2 to weights3 to use correct filter
data$wprop<-weights3$value

param <- list()
param$k <- 0.4
param$b <- 3
param$Winf <- 0.4
param$logsigma <- 0

obj <- MakeADFun(data, param,DLL="weight")
opt <- nlminb(obj$par, obj$fn, obj$gr)        

## parameter estimates
k <- opt$par[1]
b <- opt$par[2]
Winf <- opt$par[3]       


weights <- Winf * (1-exp(-k * age))^b     ### used in simulation 


weights3$wprop.est <- Winf * (1-exp(-k * weights3$Age))^b


nsshwei <- weights3 %>% 
  ggplot(aes(x=Age,y=value)) + geom_point(size= 2.5) +
  geom_line(inherit.aes=F,aes(x=Age,y=wprop.est), size= 0.6) + labs(y= "Weight (kg)", x= "Age") + theme_bw() + theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.y = element_text(size=20), axis.text.x = element_text(size=20), axis.title.x = element_text(size=30), axis.title.y = element_text(size=30))+ labs(tag="a")



 ##### Selectivity #####
fmort2<-melt(nsshfmort,id.vars="Year")
fmort2$Age <- rep(2:12,each=dim(nsshfmort)[1])

fmort2 %>% 
  ggplot(aes(x=Age,y=value, color=as.factor(Year))) + geom_point()

### ad hoc fix: everything above age 5 is fully selected
fmort2 <- fmort2 %>% mutate(sel=ifelse(Age>5,1,value/max(value[Age<5])),sel=ifelse(sel>1,1,sel))

fmort2 %>% 
  ggplot(aes(x=Age,y=sel, color=as.factor(Year))) + geom_point()

### using aggregated data
fmort3 <- fmort2 %>% group_by(Age) %>% summarise(sel=mean(sel))

### ad hoc fix: everything below age 3 is not fished
fmort3 <- fmort3 %>% mutate(sel=case_when(Age<3  ~ 0,TRUE ~ sel))


compile('Scripts/sel.cpp')
dyn.load(dynlib('Scripts/sel'))

data<-list()
data$age<-fmort3$Age
data$sel<-fmort3$sel

param <- list()
param$s50 <- 7
param$ss <- 1
param$logsigma <-0

obj <- MakeADFun(data, param,DLL="sel")
opt <- nlminb(obj$par, obj$fn, obj$gr)

## parameter estimates
s50 <- opt$par[1]
ss <- opt$par[2]


Fsel <- 1/(1+exp(-((age-s50)/ss)))  ### used in simulation 

fmort3$sel.est <- 1/(1+exp(-((fmort3$Age-s50)/ss)))

nsshfsel <- fmort3 %>% 
  ggplot(aes(x=Age,y=sel)) + geom_point(size= 2.5) +
  geom_line(inherit.aes=F,aes(x=Age,y=sel.est), size= 0.6)+ labs(y= "Fishing selectivity", x= "Age") + theme_bw() + theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.y = element_text(size=20), axis.text.x = element_text(size=20), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20))+ labs(tag="a") + scale_x_continuous(breaks=c(2,6, 10))




#######################################################################
######---------------      NEA COD           --------------- #############
#######################################################################


age <- c(3:15)     ## maximum age 3-15 
recage<- 3

start.rec <- mean(codass$Rec.age3)*1000

#### Recruitment: important to correct for the rec.age offset ####

tsl<-length(codass$Year)
Rec <- codass$Rec.age3[(recage+1):tsl]*1000      
SSB <- codass$SSBtonnes[1:(tsl-(recage))]

Rec <- Rec/1000000000  # for plotting
SSB <- SSB/1000000     #for plotting


## 5. Ricker 
data<-list()
data$SSB<-SSB #c(0,SSB)
data$Rec<-Rec #c(0,rep(mean(Rec),length(Rec)))


#nlsricker<- nls(Rec~alpha*SSB*exp(beta*SSB),data=data,start=list(alpha=1e7,beta=1e-7))
nlsrssb <- nls(log(Rec/SSB) ~ alpha+(beta*SSB),data=data,start=list(alpha=5,beta=-1), na.action = na.omit)

alphar <- summary(nlsrssb)$coefficients[1]
betar <- summary(nlsrssb)$coefficients[2]


#Recruitment5 <- exp(alpha)*data$SSB*exp(beta*data$SSB)
RSSB5 <- alphar +betar *data$SSB
ricker.sd<-sd(residuals(nlsrssb))
ricker.mean<-mean(residuals(nlsrssb))
autocorrelation<-acf(residuals(nlsrssb))
AR1par <-autocorrelation$acf[2]


# Theoretical SSB1
SSB1<-seq(0, 6e6, 1e3)

SSB1<- SSB1/1000000 # for plotting

Recruitment5 <- exp(alphar)*SSB1*exp(betar *SSB1)  #Deterministic stock-recruitment version

Rvariation<-rnorm(length(SSB1),0,ricker.sd)        #Norm. dist. error term (on log-scale)
Recruitment5.1 <- exp(alphar)*SSB1*exp(betar *SSB1) * exp(Rvariation) #- Stochasticity  

Rvariationacf<-c(Rvariation[1],AR1par*Rvariation[1:(length(SSB1)-1)]+Rvariation[2:length(SSB1)]) #-  autocorrelation
Recruitment5.2 <- exp(alphar)*SSB1 *exp(betar *SSB1) *exp(Rvariationacf) #- stochasticity 

plot(Recruitment5~SSB1,type="l",col="red", xlim=c(0,6), ylim=c(0,3), xlab="SSB1 million t",  ylab="Rec billions")
points(Recruitment5.1~SSB1,col="blue",pch=1)
points(Recruitment5.2~SSB1,col="purple",pch=2)
points(Rec~SSB,pch=16)




##---- Compare with assessed SSB --- SSB real data 
Recruitment5 <- exp(alphar)*SSB*exp(betar*SSB)     #- Deterministic stock-recruitment version

Rvariation<-rnorm(length(SSB),0,ricker.sd)         #- Norm. dist. error term (on log-scale)
Recruitment5.1 <- exp(alphar)*SSB*exp(betar *SSB) * exp(Rvariation)   #- stochasticity 

Rvariationacf<-c(Rvariation[1],AR1par*Rvariation[1:(length(SSB)-1)]+Rvariation[2:length(SSB)]) #- AR1 autocorr to error term
Recruitment5.2 <- exp(alphar)*SSB *exp(betar *SSB) *exp(Rvariationacf)  #- stochasticity


plot(Recruitment5~SSB,col="red", type="l", xlab="SSB1 million t",  ylab="Rec billions")
points(Recruitment5.1~SSB,col="blue",pch=1)
points(Recruitment5.2~SSB,col="purple",pch=2)
points(Rec~SSB,pch=16)

c(sum((Rec-Recruitment5)^2),sum((Rec-Recruitment5.1)^2),sum((Rec-Recruitment5.2)^2))


###### BH (trying again 2/2 2022)
# 1/R = beta + alpha * 1/SSB algebraic transformation of BH

#run first lines of NEA rec to load data

data<-list()
data$SSB<-SSB #c(0,SSB)
data$Rec<-Rec #c(0,rep(mean(Rec),length(Rec)))

## BH algebraic transformation, estimating base of ICES data  
bhalge <- nls(1/Rec ~ beta + alpha * (1/SSB), data=data, start = list(alpha=5, beta=-1)) 

alphabh <- summary(bhalge)$coefficients[1]
betabh <- summary(bhalge)$coefficients[2]

# SSB1 theoretical numbers to test alpha & beta estimates 
SSB1 <- seq(0,9e6,1e5)
#SSB1 <- SSB1/1000000  #for plotting

codRec6 <- betabh + alphabh * (1/SSB) # testing the alpha and beta, works!

plot(1/codRec6 ~ SSB,type="l",col="red")

#######   dont use any of this 
## BH - alpha and beta estimated in TMB - 2 options ##
#Recruitment6 <- exp(alphabh +log(SSB)-log(exp(betabh)*SSB))

#Recruitment6 <- (alphabh*SSB)/(1+betabh*SSB)   ### dont use this
#plot(Recruitment6 ~ SSB,type="l",col="red")




## plot
rec.vector<-c(Rec,Recruitment5,Recruitment5.1, Recruitment5.2)
type<-rep(c("Real","Ricker","Ricker5.1", "Ricker5.2"),each=length(Rec))
rec.df<-data.frame(Type=type,SSB=rep(SSB,4),Rec=rec.vector)

rec.df %>% ggplot(aes(x=SSB,y=Rec,color=Type)) + geom_point(size=3) +
  scale_color_brewer(palette="Accent") + theme_bw()


##### BEVERTON-HOLT

# ## 6. BH  - NB: fit depends heavily on starting values and sucks (tends to become constant)
# compile('Scripts/bh.cpp')
# dyn.load(dynlib('Scripts/bh'))
# 
# data<-list()
# data$ssb<-SSB
# data$logR<-log(Rec)
# 
# param <- list()
# param$loga <- 1
# param$logb <- 1
# param$logsigma <-0
# 
# obj <- MakeADFun(data, param,DLL="bh")
# optbh <- nlminb(obj$par, obj$fn, obj$gr)
# 
# Recruitment6 <- (optbh$par[1]*SSB)/(1+optbh$par[2]*SSB) #exp(optbh$par[1]+log(SSB)-log(exp(optbh$par[2])*SSB))
# 
# alphabh <- optbh$par[1]
# betabh <- optbh$par[2]
# 
# 
# plot(Recruitment6~SSB, col="red", pch=3)
# 
# ## plot
# rec.vector<-c(Rec, Recruitment5, Recruitment6)
# type<-rep(c("Real","Ricker","BH"),each=length(Rec))
# rec.df<-data.frame(Type=type,SSB=rep(SSB,3),Rec=rec.vector)
# 
# rec.df %>% ggplot(aes(x=SSB,y=Rec,color=Type)) + geom_point(size=3) +
#   scale_color_brewer(palette="Accent") + theme_bw()

                

   ##### BEVERTON-HOLT continued (another way)

data<-list()
data$SSB<-SSB #c(0,SSB)
data$Rec<-Rec #c(0,rep(mean(Rec),length(Rec)))


nlsrssb <- nls(1/Rec ~ beta + alpha / SSB, data=data,start=list(alpha=5,beta=-1))

alphabh <- summary(nlsrssb)$coefficients[1]
betabh <- summary(nlsrssb)$coefficients[2]


#Recruitment6 <- exp(alpha)*data$SSB*exp(beta*data$SSB)
RSSB6 <- SSB/(alphabh +betabh *SSB)
bh.sd <- 1/sd(residuals(nlsrssb))
bh.mean <- mean(residuals(nlsrssb))
autocorrelation <- acf(residuals(nlsrssb))
AR1par <- autocorrelation$acf[2]


# Theoretical SSB1
SSB1<-seq(0, 4e6, 1e5)

Recruitment6 <- SSB1/(alphabh +betabh *SSB1)  #Deterministic stock-recruitment version


Rvariation<-rnorm(length(SSB1),0,bh.sd)         #Norm. dist. errorterm (on log-scale)
Recruitment6.1 <- SSB1/(alphabh +betabh *SSB1) + Rvariation   #- stochasticity 


Rvariationacf<-c(Rvariation[1],AR1par*Rvariation[1:(length(SSB1)-1)]+Rvariation[2:length(SSB1)]) #- AR1 autocorr to error term
Recruitment6.2 <- SSB1/(alphabh +betabh *SSB1) + Rvariationacf  #- stochasticity


plot(Recruitment6~SSB1,col="red", type="l",ylim=c(0,max(Rec)*1.2))
points(Recruitment6.1~SSB1,col="blue",pch=1)
points(Recruitment6.2~SSB1,col="purple",pch=2)
points(Rec~SSB,pch=16)


##---- Compare with assessed SSB --- SSB real data 
Recruitment6 <- SSB/(alphabh +betabh *SSB)    #Deterministic stock-recruitment version


Rvariation<-rnorm(length(SSB),0,bh.sd)         #Norm. dist. errorterm (on log-scale)


Recruitment6.1 <- SSB/(alphabh +betabh *SSB) + Rvariation   #- stochasticity 


Rvariationacf<-c(Rvariation[1],AR1par*Rvariation[1:(length(SSB)-1)]+Rvariation[2:length(SSB)]) #- AR1 autocorr to error term
Recruitment6.2 <- SSB/(alphabh +betabh *SSB) + Rvariationacf  #- stochasticity



plot(Recruitment6~SSB,col="red", pch=3,ylim=c(0,max(Rec)*1.2))
points(Recruitment6.1~SSB,col="blue",pch=1)
points(Recruitment6.2~SSB,col="purple",pch=2)
points(Rec~SSB,pch=16)

c(sum((Rec-Recruitment5)^2),sum((Rec-Recruitment5.1)^2),sum((Rec-Recruitment5.2)^2))


## plot
rec.vector<-c(Rec, Recruitment5, Recruitment6)
type<-rep(c("Real","Ricker","BH"),each=length(Rec))
rec.df<-data.frame(Type=type,SSB=rep(SSB,3),Rec=rec.vector)

rec.df %>% ggplot(aes(x=SSB,y=Rec,color=Type)) + geom_point(size=3) +
scale_color_brewer(palette="Accent") + theme_bw()




##### Maturity #####

maturity2 <-melt(codmaturity,id.vars="Year")
maturity2$Age <- rep(3:15,each=dim(codmaturity)[1])

maturity2 %>% 
  ggplot(aes(x=Age,y=value,color=as.factor(Year))) + geom_line()

compile('Scripts/maturity.cpp')
dyn.load(dynlib('Scripts/maturity'))

data<-list()
data$age<-maturity2$Age
data$mprop<-maturity2$value

param <- list()
param$a50 <- 4
param$env <- .2
param$logsigma <-0

obj <- MakeADFun(data, param,DLL="maturity")
opt <- nlminb(obj$par, obj$fn, obj$gr)

## parameter estimates
a50 <- opt$par[1]
env <- opt$par[2]

maturity <- round(1/(1+exp(-((age-a50)/env))),2)   ## used in simulation 


maturity2$mprop.est <- 1/(1+exp(-((maturity2$Age-a50)/env)))

codmat<-maturity2 %>% 
  ggplot(aes(x=Age,y=value)) + geom_point(size= 2.5) +
  geom_line(inherit.aes=F,aes(x=Age,y=mprop.est), size= 0.6) + labs(y= "Maturity", x= "Age") + theme_bw() + theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.y = element_text(size=20), axis.text.x = element_text(size=20), axis.title.x = element_text(size=30), axis.title.y = element_text(size=30))+ labs(tag="b")



##### Weights #####
weights2 <-melt(codweight,id.vars="Year")
weights2$Age <- rep(3:15,each=dim(codweight)[1])


weights2 %>% 
  ggplot(aes(x=Age,y=value,color=as.factor(Year))) + geom_line()

compile('Scripts/weight.cpp')
dyn.load(dynlib('Scripts/weight'))

data<-list()
data$age<-weights2$Age
data$wprop<-weights2$value

param <- list()
param$k <- 0.4
param$b <- 3
param$Winf <- 0.4
param$logsigma <- 0

obj <- MakeADFun(data, param,DLL="weight")
opt <- nlminb(obj$par, obj$fn, obj$gr)        

## parameter estimates
k <- opt$par[1]
b <- opt$par[2]
Winf <- opt$par[3]       


weights <- Winf * (1-exp(-k * age))^b     ### used in simulation 


weights2$wprop.est <- Winf * (1-exp(-k * weights2$Age))^b


codwei<-weights2 %>% 
  ggplot(aes(x=Age,y=value)) + geom_point(size= 2.5) +
  geom_line(inherit.aes=F,aes(x=Age,y=wprop.est), size= 0.6) + labs(y= "Weight (kg)", x= "Age") + theme_bw() + theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.y = element_text(size=20), axis.text.x = element_text(size=20), axis.title.x = element_text(size=30), axis.title.y = element_text(size=30))+ labs(tag="b")




##### Selectivity #####

fmort2<-melt(codfmort,id.vars="Year")
fmort2$Age <- rep(3:15,each=dim(codfmort)[1])

fmort2 %>% 
  ggplot(aes(x=Age,y=value,color=as.factor(Year))) + geom_point()


### using aggregated data
fmort3 <- fmort2 %>% group_by(Age) %>% summarise(fmean=mean(value))

fmort3 <- fmort3 %>% mutate(sel=fmean/max(fmean))


compile('Scripts/sel.cpp')
dyn.load(dynlib('Scripts/sel'))

data<-list()
data$age<-fmort3$Age
data$sel<-fmort3$sel

param <- list()
param$s50 <- 7
param$ss <- 1
param$logsigma <-0

obj <- MakeADFun(data, param,DLL="sel")
opt <- nlminb(obj$par, obj$fn, obj$gr)

## parameter estimates
s50 <- opt$par[1]
ss <- opt$par[2]


Fsel <- 1/(1+exp(-((age-s50)/ss)))  ### used in simulation 

fmort3$sel.est <- 1/(1+exp(-((fmort3$Age-s50)/ss)))

codfsel<- fmort3 %>% 
  ggplot(aes(x=Age,y=sel)) + geom_point(size= 2.5) +
  geom_line(inherit.aes=F,aes(x=Age,y=sel.est), size= 0.6)+ labs(y= "Fishing selectivity", x= "Age") + theme_bw() + theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.y = element_text(size=20), axis.text.x = element_text(size=20), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20))+ labs(tag="b")




#######################################################################
######---------------      BEAKED REDFISH           --------------- #############
#######################################################################

recage <- 2
amax <- 40    
age <- recage:amax

start.rec <- mean(redfishass2$Rec.age2thousand)*1000

#### Recruitment ####
tsl <- length(redfishass2$Year)  #old assesment 2018
Rec <- redfishass2$Rec.age2thousand[(recage+1):tsl]*1000  #old assesment 2018
SSB <- redfishass2$StockBiomass.t[1:(tsl-(recage))]    #old assesment 2018


#tsl <- length(redfishass$Year)  #old assesment 2019
#Rec <- redfishass$Rec.age2.1000[(recage+1):tsl]*1000  #old assesment 2019
#SSB <- redfishass$SSB.tonnes[1:(tsl-(recage))] 


Rec <- Rec/1000000000  # for plotting million
SSB <- SSB/1000000     #for plotting thousands


## 5. Ricker 
data<-list()
data$SSB<-SSB #c(0,SSB)
data$Rec<-Rec #c(0,rep(mean(Rec),length(Rec)))

#nlsricker<- nls(Rec~alpha*SSB*exp(beta*SSB),data=data,start=list(alpha=1e7,beta=1e-7))
nlsrssb <- nls(log(Rec/SSB) ~ alpha+(beta*SSB),data=data,start=list(alpha=5,beta=-1))

alphar <- summary(nlsrssb)$coefficients[1]
betar <- summary(nlsrssb)$coefficients[2]


Recruitment6 <- exp(alphar)*data$SSB*exp(betar*data$SSB)
RSSB5 <- alphar +betar *data$SSB
ricker.sd <-sd(residuals(nlsrssb))
ricker.mean<-mean(residuals(nlsrssb))
autocorrelation<-acf(residuals(nlsrssb))
AR1par <-autocorrelation$acf[2]


# Theoretical SSB1
SSB1<-seq(0,9e6,1e5)

SSB1 <- SSB1/1000000 #for plotting

redRecruitment5 <- exp(alphar)*SSB1*exp(betar *SSB1)  #- Deterministic stock-recruitment version
Rvariation<-rnorm(length(SSB1),0,ricker.sd)        #- Norm dist. error term (on log-scale)
Recruitment5.1 <- exp(alphar)*SSB1*exp(betar *SSB1) * exp(Rvariation) #- stochasticity

Rvariationacf<-c(Rvariation[1],AR1par*Rvariation[1:(length(SSB1)-1)]+Rvariation[2:length(SSB1)]) #-  autocorrelation
Recruitment5.2 <- exp(alphar)*SSB1 *exp(betar *SSB1) *exp(Rvariationacf) #- stochasticity 

plot(redRecruitment5_2018~SSB1,type="l",col="blue") #xlim=c(0, 3000000), ylim=c(0, 10e+08),  xlab="SSB1 million t",  ylab="Rec billions")
points(Recruitment5.1~SSB1,col="blue",pch=1)
points(Recruitment5.2~SSB1,col="purple",pch=2)
points(Rec~SSB,pch=16)

points(y=Rec, x=SSB)
points(y=1/redRec6_2018, x=SSB1, col="blue")
points(y=1/redRec6_2019, x=SSB1, col="red")
points(y=redRecruitment5_2019, x=SSB1, col="red", pch=3)


##---- Compare with assessed SSB --- SSB real data 
Recruitment5 <- exp(alphar)*SSB*exp(betar*SSB)     #- Deterministic stock-recruitment version
Rvariation<-rnorm(length(SSB),0,ricker.sd)         #- Norm. dist. error term (on log-scale)
Recruitment5.1 <- exp(alphar)*SSB*exp(betar *SSB) * exp(Rvariation)   #- stochasticity 

Rvariationacf<-c(Rvariation[1],AR1par*Rvariation[1:(length(SSB)-1)]+Rvariation[2:length(SSB)]) #- AR1 autocorr to error term
Recruitment5.2 <- exp(alphar)*SSB *exp(betar *SSB) *exp(Rvariationacf)  #- stochasticity

plot(redRecruitment5~SSB1,type="l",col="red") # xlim= c(0, 10e+05), ylim=c(0, 2e+09))
points(Recruitment5~SSB,col="blue",pch=1)
points(Recruitment5.2~SSB,col="purple",pch=2)
points(Rec~SSB,pch=16)

c(sum((Rec-Recruitment5)^2),sum((Rec-Recruitment5.1)^2),sum((Rec-Recruitment5.2)^2))


###### BH (trying again 2/2 2022)
# 1/R = beta + alpha * 1/SSB algebraic transformation of BH

tsl <- length(redfishass2$Year)  #old assesment 2018
Rec <- redfishass2$Rec.age2thousand[(recage+1):tsl]*1000  #old assesment 2018
SSB <- redfishass2$StockBiomass.t[1:(tsl-(recage))]    #old assesment 2018

#tsl <- length(redfishass$Year)  #old assesment 2019
#Rec <- redfishass$Rec.age2.1000[(recage+1):tsl]*1000  #old assesment 2019
#SSB <- redfishass$SSB.tonnes[1:(tsl-(recage))] 

Rec <- Rec/1000000000  # for plotting million
SSB <- SSB/1000000     #for plotting thousands

data<-list()
data$SSB<-SSB #c(0,SSB)
data$Rec<-Rec #c(0,rep(mean(Rec),length(Rec)))

## BH algebraic transformation, estimating base of ICES data  
bhalge <- nls(1/Rec ~ beta + alpha * (1/SSB), data=data, start = list(alpha=5, beta=-1)) 

alphabh <- summary(bhalge)$coefficients[1]
betabh <- summary(bhalge)$coefficients[2]

# SSB1 theoretical numbers to test alpha & beta estimates 
SSB1 <- seq(0,9e6,1e5)
SSB1 <- SSB1/1000000  #for plotting

redRecruitment6 <- betabh + alphabh * (1/SSB1) # testing the alpha and beta, works!

plot(1/redRecruitment6 ~ SSB1,type="l",col="red")
points(x=SSB, y=Rec)



#######   dont use any of this 
## BH - alpha and beta estimated in TMB - 2 options ##
#Recruitment6 <- exp(alphabh +log(SSB)-log(exp(betabh)*SSB))

#Recruitment6 <- (alphabh*SSB)/(1+betabh*SSB)   ### dont use this
#plot(Recruitment6 ~ SSB,type="l",col="red")



#---- Plotting all ----#
rec.vector <- c(redRecruitment5, 1/redRecruitment6)
type <-rep(c("Ricker5", "BH"),each=length(redRecruitment5))
rec.df <- data.frame(Type=type, SSB=rep(SSB1,2),Rec=rec.vector)

rec.df %>% ggplot(aes(x=SSB,y=Rec, color=Type)) + geom_line(size=1.5) +
  scale_color_brewer(palette="Accent") +
  theme_bw() + theme(panel.background = element_blank(), 
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"), axis.text.y = element_text(size=15), 
                     axis.text.x = element_text(size=15), axis.title.x = element_text(size=20), 
                     axis.title.y = element_text(size=20))


##### Maturity #####

maturity <-melt(redfishmaturity,id.vars="Year")
maturity$Age <- rep(6:19,each=dim(redfishmaturity)[1])

maturity %>% 
  ggplot(aes(x=Age,y=value,color=as.factor(Year))) + geom_line()

compile('Scripts_R/maturity.cpp')
dyn.load(dynlib('Scripts_R/maturity'))

data<-list()
data$age<-maturity$Age
data$mprop<-maturity$value

param <- list()
param$a50 <- 4
param$env <- .2
param$logsigma <-0

obj <- MakeADFun(data, param,DLL="maturity")
opt <- nlminb(obj$par, obj$fn, obj$gr)

## parameter estimates
a50 <- opt$par[1]
env <- opt$par[2]

maturity <- round(1/(1+exp(-((unique(age)-a50)/env))),2)   ## used in simulation 

maturity$mprop.est <- 1/(1+exp(-((maturity$Age-a50)/env)))

redmat<- maturity %>% 
  ggplot(aes(x=Age,y=value)) + geom_point(size= 2.5) +
  geom_line(inherit.aes=F,aes(x=Age,y=mprop.est), size= 0.6) +
  labs(y= "Maturity", x= "Age") + 
  theme_bw() + theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.y = element_text(size=20), axis.text.x = element_text(size=20), axis.title.x = element_text(size=30), axis.title.y = element_text(size=30)) + labs(tag="c")


maturityforplot <- list() 
maturityforplot$Age <- c(6:40)
maturityforplot$mprop.est <- 1/(1+exp(-((c(6:40)-a50)/env)))
maturityforplot <- as.data.frame(maturityforplot)


redmat <- ggplot(data=maturityforplot, inherit.aes=F, aes(x=Age,y=mprop.est)) + 
  geom_line(color="red",size= 0.6) +
  labs(y= "Maturity", x= "Age") + 
  theme_bw() + theme(panel.background = element_blank(), 
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"), axis.text.y = element_text(size=20), 
                     axis.text.x = element_text(size=20), axis.title.x = element_text(size=30),
                     axis.title.y = element_text(size=30)) + labs(tag="c")

redmat + geom_point(data= maturity,aes(x=Age,y=value), size= 2.5)




##### Weights #####
weights2 <-melt(redfishweight,id.vars="Year")
weights2$Age <- rep(6:19,each=dim(redfishweight)[1])

weights2 %>% 
  ggplot(aes(x=Age,y=value,color=as.factor(Year))) + geom_line()

compile('Scripts_R/weight.cpp')
dyn.load(dynlib('Scripts_R/weight'))

data<-list()
data$age<-weights2$Age
data$wprop<-weights2$value

param <- list()
param$k <- 0.4
param$b <- 3
param$Winf <- 0.4
param$logsigma <- 0

obj <- MakeADFun(data, param,DLL="weight")
opt <- nlminb(obj$par, obj$fn, obj$gr)       

## parameter estimates
k <- opt$par[1]
b <- opt$par[2]
Winf <- opt$par[3]       

weights <- Winf * (1-exp(-k * age))^b     ### used in simulation 

weights2$wprop.est <- Winf * (1-exp(-k * weights2$Age))^b


redwei<- weights2 %>% 
  ggplot(aes(x=Age,y=value)) + geom_point(size= 3) +
  geom_line(inherit.aes=F,aes(x=Age,y=wprop.est), size= 0.6) + 
  labs(y= "Weight (kg)", x= "Age") + theme_bw() + theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.y = element_text(size=20), axis.text.x = element_text(size=20), axis.title.x = element_text(size=30), axis.title.y = element_text(size=30)) + labs(tag="c")


weightforplot <- list() 
weightforplot$Age <- c(6:40)
weightforplot$wprop.est <- Winf * (1-exp(-k * (c(6:40))))^b 
weightforplot <- as.data.frame(weightforplot)


redwei <- ggplot(data=weightforplot, inherit.aes=F, aes(x=Age,y=wprop.est)) + 
  geom_line(color="red",size= 0.6)  + 
  theme_bw() + theme(panel.background = element_blank(), 
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"), axis.text.y = element_text(size=20), 
                     axis.text.x = element_text(size=20), axis.title.x = element_blank(),
                     axis.title.y = element_blank()) 

redwei + geom_point(data= weights2, aes(x=Age,y=value), size= 2.5)


##### Selectivity #####

fmort2<-melt(redfishfmort,id.vars="Year")
fmort2$Age <- rep(2:19,each=dim(redfishfmort)[1])

fmort2 %>% 
  ggplot(aes(x=Age,y=value,color=as.factor(Year))) + geom_point()


### using aggregated data
fmort3 <- fmort2 %>% group_by(Age) %>% summarise(fmean=mean(value))

fmort3 <- fmort3 %>% mutate(sel=fmean/max(fmean))


compile('Scripts/sel.cpp')
dyn.load(dynlib('Scripts/sel'))
data<-list()
data$age<-fmort3$Age
data$sel<-fmort3$sel

param <- list()
param$s50 <- 7
param$ss <- 1
param$logsigma <-0

obj <- MakeADFun(data, param,DLL="sel")    # MakeADFun - automatic differentiation function 
opt <- nlminb(obj$par, obj$fn, obj$gr)

## parameter estimates
s50 <- opt$par[1]
ss <- opt$par[2]

Fsel <- 1/(1+exp(-((age-s50)/ss)))  ### used in simulation 

fmort3$sel.est <- 1/(1+exp(-((fmort3$Age-s50)/ss)))

redfsel <- fmort3 %>% 
  ggplot(aes(x=Age,y=sel)) + geom_point(size= 2.5) +
  geom_line(inherit.aes=F,aes(x=Age,y=sel.est), size= 0.6) + 
  labs(y= "Fishing selectivity", x= "Age") + 
  theme_bw() + theme(panel.background = element_blank(), 
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"), 
                     axis.text.y = element_text(size=20), 
                     axis.text.x = element_text(size=20), 
                     axis.title.x = element_text(size=20), 
                     axis.title.y = element_text(size=20)) + labs(tag="c")



fselforplot <- list() 
fselforplot$Age <- c(2:40)
fselforplot$fselprop.est <- 1/(1+exp(-(((c(2:40))-s50)/ss)))
fselforplot <- as.data.frame(fselforplot)

redfsel <- ggplot(data=fselforplot, inherit.aes=F, aes(x=Age,y=fselprop.est)) + 
  geom_line(color="red",size= 0.6)  + 
  theme_bw() + theme(panel.background = element_blank(), 
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"), axis.text.y = element_text(size=20), 
                     axis.text.x = element_text(size=20), axis.title.x = element_blank(),
                     axis.title.y = element_blank()) 

redfsel + geom_point(data= fmort3, aes(x=Age,y=sel), size= 2.5)



grid.arrange(nsshmat, codmat, redmat, ncol=3)
grid.arrange(redwei, nsshwei, codwei)
grid.arrange(nsshfsel, codfsel, redfsel)



ggsave(filename="fsel_nssh.pdf",
       plot=last_plot(),
       width = 100, 
       height = 80, 
       units = "mm")


###### PLOTTINg all ############
maturitynssh$species <- "NSSH"
maturitycod$species <- "cod"
maturityred$species <- "redfish"


maturitynssh <- as.data.frame(maturitynssh)
maturitycod <- as.data.frame(maturitycod)
maturityred <- as.data.frame(maturityred)

cbind(maturitycod, maturitynssh, maturityred, capematurity)


map2(maturitynssh, maturitycod, maturityred, capmat, left_join)


matall <- merge(maturitynssh, maturitycod, maturityred, by.x = "mprop.est", by.y = "Age", by.z = "species", all.x = T, all.y = T, all.z= T)

params <- fmort3 %>% 
  ggplot(aes(x=Age,y=sel)) + geom_point() +
  geom_line(inherit.aes=F,aes(x=Age,y=sel.est))+ labs(y= "Fishing selectivity", x= "age") + theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 

params <- params + geom_point(aes(y= )) 

grid.arrange(redmat, nsshmat, codmat)
grid.arrange(redwei, nsshwei, codwei)
grid.arrange(redfsel, nsshfsel, codfsel)







