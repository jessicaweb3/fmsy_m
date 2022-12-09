
# F. Zimmermann & J. Tengvall 
#last edit 13/10 2021

#### Simulating M populations   ####

library(tidyverse)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(reshape2)
library(openxlsx)
library(TMB)
library(here) #setting working directory - here for Rproj
library(furrr)
library(data.table) 

#library(purrr)   #mac library for map_df

 
# Function for simulation "am= age model" - [Year= rows, age= col] One N or Ca matrix per M (and F), all F runs per M.

am <- function(r, sp, start.rec, Fsel, multiplier, tmax, maturity, weights, boot, fseq, fmax) #Mvariation
{
  ac<-length(age) # number of age classes
  N<-matrix(nrow=tmax+1,ncol=ac)     ## Fished pop. matrix
  Ca<-matrix(nrow=tmax+1,ncol=ac)    ## catch matrix per age and year     
  Ct<- NA     ## Catch per year  kg
  SSB <- NA                            ## SSB in t
  SSB1 <- NA       #use to test SSB > SSB[year1]/1000
  Recruitment <- NA
  Rvariation <- NA
  
  Caamean <- matrix(nrow=tmax+1, ncol=ac)
  
  Finput<- seq(0, fmax, length=fseq)  ## increments of fishing mortality
  
  # Overall means and references 
  Cmean<-vector(length=length(Finput))
  SSBmean<-vector(length=length(Finput))
  Recmean <- vector(length=length(Finput))
  Collapse<-vector(length=length(Finput))
  Fbar <- vector(length=length(Finput))
  Mbar <- vector(length=length(Finput))
  SSBcv <- vector(length=length(Finput))
  Reccv <- vector(length=length(Finput))
  Blim.any <- vector(length=length(Finput))
  Blim.mean <- vector(length=length(Finput))
  Blim.prop <- vector(length=length(Finput))
  
  outputlogMsd.df <- data.table()  # data.table instead of df saves space  
  outputTimeseries.df <- data.frame()
  
  catchatage.df <- data.table()
  
  # M vector 
  Mfinal <- Nmort*multiplier  # Option 1: * multiplier ## age-dependent M vector * multiplier
  

  #### Replicates of fishing rates  ####

  for(b in 1:boot) {
    for(f in seq_along(Finput)) { 
      
      Finput1 <- Finput[f] # creates Fbar per age
      
      ## Year 1 #no fishing
      if(sp %in% c("nssh","red")){N[1,] <-  start.rec * exp(-Mfinal*1:ac)      
      N[2,] <- start.rec * exp(-Mfinal*1:ac)
      
      SSB[1:2] <- sum(N[1,] * maturity * weights)/1000       # SSB ton
      Rvariation[1:2] <- c(0,0)}
      
      else if(sp=="cod"){N[1,] <-  start.rec* exp(-Mfinal*1:ac)      
      N[2,] <-  start.rec * exp(-Mfinal*1:ac)
      N[3,] <-  start.rec * exp(-Mfinal*1:ac)      
      
      SSB[1:3] <- sum(N[1,] * maturity * weights)/1000       # SSB ton
      Rvariation[1:3] <- c(0,0, 0)}

      
      # Time loop, all other years
      for(t in recage:tmax) { 
        
        #Mfinish <- Mfinal*rnorm(1, 1, Mvariation)  # varying M with 10% within time
        
        Fmort <- (ifelse(t < 51, averageF, Finput1)) * Fsel
        
        # Ricker recruitment nssh & redfish
        if (r==5 & sp %in% c("nssh","red")){
          Recruitment[t+1] <- exp(alphar) * SSB[t-1] * exp(betar * SSB[t-1])
        }  # exp(alphar)*SSB1*exp(betar *SSB1) 
        
        else if (r==5.1 & sp %in% c("nssh","red")){Rvariation[t+1]<-rnorm(1,0,ricker.sd) #### Rec at age 2 == 2 years lag == SSB at t-1 determines Rec at t+1!
        Recruitment[t+1] <- exp(alphar) * SSB[t-1] * exp(betar * SSB[t-1]) *exp(Rvariation[t+1])
        while(Recruitment[t+1] > reccap) {         #can skip reccap as vals. never hit it
          Rvariation[t+1]<-rnorm(1,0,ricker.sd)
          Recruitment[t+1] <- exp(alphar) * SSB[t-1] * exp(betar * SSB[t-1]) *exp(Rvariation[t+1]) }
        #if(Recruitment[t+1]<recmin) Recruitment[t+1] <- recmin
        }
        
        else if (r==5.2 & sp %in% c("nssh","red")){Rvariation[t+1] <- rnorm(1,0,ricker.sd)+AR1par*Rvariation[t]
        Recruitment[t+1] <-exp(alphar) * SSB[t-1] * exp(betar * SSB[t-1]) *exp(Rvariation[t+1])
        while(Recruitment[t+1]> reccap) {      #can skip reccap as vals. never hit it
          Rvariation[t+1]<-rnorm(1,0,ricker.sd)
          Recruitment[t+1] <- exp(alphar) * SSB[t-1] * exp(betar * SSB[t-1]) *exp(Rvariation[t+1]) }
        #if(Recruitment[t+1]<recmin) Recruitment[t+1] <- recmin
        }
      
        ###BH recruitment algebraic solution NOT TMB
        else if (r==6 & sp %in% c("nssh","red")){Recruitment[t+1] <- 1 / (betabh + alphabh * 1/SSB[t-1])
        }
        #Recruitment[t+1] <- SSB[t-1] * (alphabh + betabh * SSB[t-1])
        #1/Rec ~ beta + alpha * (1/SSB)
        #1/Recruitment[t+1] <- betabh + alphabh * (1/SSB[t-1]) - algebraic transformation
        
        #BH rec + stochastisy
        else if (r==6.1 & sp %in% c("nssh","red")){
          Rvariation[t+1] <- rnorm(1, 0, bh.sd)
          Recruitment[t+1]<- 1/ (betabh + alphabh * 1/SSB[t-1]) * (1/Rvariation[t+1])
        }
        
        # Ricker recruitment cod
        else if (r==5 & sp=="cod"){Recruitment[t+1] <- exp(alphar) * SSB[t-2] * exp(betar * SSB[t-2])
        }
        
        else if (r==5.1 & sp=="cod"){Rvariation[t+1]<-rnorm(1,0,ricker.sd) #### Rec at age 2 == 2 years lag == SSB at t-1 determines Rec at t+1!
        Recruitment[t+1] <- exp(alphar) * SSB[t-2] * exp(betar * SSB[t-2]) *exp(Rvariation[t+1])
        
        #while(Recruitment[t+1]> reccap) {         #can skip reccap as vals. never hit it
         # Rvariation[t+1]<-rnorm(1,0,ricker.sd)
          #Recruitment[t+1] <- exp(alphar) * SSB[t-2] * exp(betar * SSB[t-2]) *exp(Rvariation[t+1]) }
        #if(Recruitment[t+1]<recmin) Recruitment[t+1] <- recmin
        }
        
        else if (r==5.2 & sp=="cod"){Rvariation[t+1] <- rnorm(1,0,ricker.sd) + AR1par *Rvariation[t]
        Recruitment[t+1] <-exp(alphar) * SSB[t-2] * exp(betar * SSB[t-2]) *exp(Rvariation[t+1])
        while(Recruitment[t+1]> reccap) {         #can skip reccap as vals. never hit it
          Rvariation[t+1]<-rnorm(1,0,ricker.sd)
          Recruitment[t+1] <- exp(alphar) * SSB[t-1] * exp(betar * SSB[t-1]) *exp(Rvariation[t+1]) }
        #if(Recruitment[t+1]<recmin) Recruitment[t+1] <- recmin
        }
        
        else if (r==6 & sp=="cod"){Recruitment[t+1] <- 1 / (betabh + alphabh * 1/SSB[t-1])
        }
        
        Recruitment[t+1]<- Recruitment[t+1]   #billions 
        N[t+1,1] <- Recruitment[t+1] 
        
        
        # All other ages 
        survival <- exp(-Mfinal-Fmort)
        
        N[t+1,2:ac] <- N[t,1:(ac-1)] * survival[1:(ac-1)]
        
        
        N[t+1,ac]<- N[t+1,ac] + N[t,ac] * survival[ac]
        
        SSB1[t+1]<- sum(N[t+1,] * maturity * weights)/1000
        
        SSB[t+1] <- ifelse(SSB1[t+1] < SSB[1]/1000, 0, SSB1[t+1])   
        #It should be SSB < SSB[1]/1000 we can also compare it with the "balanced" stock when we start the actual simulation and take, let's say, 1/100 of that. 

        
        # Catch matrix
        Ca[t,]<- N[t,] * (1-exp(-(Fmort+Mfinal)))*Fmort/(Fmort+Mfinal)  
        Ct[t]<- sum(Ca[t,] * weights)/1000                        
        if(t==tmax) {
          Ca[t+1,]<- N[t+1,] * (1-exp(-(Fmort+Mfinal)))*Fmort/(Fmort+Mfinal)  
          Ct[t+1]<- sum(Ca[t+1,] * weights)/1000         
        }
        
      }
      
      # Means of last 100 years per f
      Cmean[f] <- mean(Ct[(tmax-100):(tmax)]) 
      
      SSBmean[f] <- mean(SSB[(tmax-100):(tmax)]) 
      
      Recmean[f] <- mean(N[(tmax-100:tmax),1])
      
      
      Fbar[f] <- mean(Fmort[fbar[1]:fbar[2]])  # Fbar mean of last ages (3:amax is age group 5-13) 
      
      Mbar[f] <- mean(Mfinal[mbar[1]:mbar[2]])  ## FZ: fbar contains the upper and lower age class, defined in the parameters
      
      Blim.any[f] <- ifelse(any(SSB[(tmax-100):(tmax)] < Blim) ,1, 0) ## risk of SSB falling in any year below Blim
      Blim.mean[f] <- ifelse(SSBmean[f] < Blim, 1, 0) ## risk of mean SSB being below Blim
      
      Blim.prop[f] <- sum(SSB[(tmax-100):(tmax)] < Blim)
      
      Caamean[f] <- mean(Ca[,])
      
      
  outputTimeseries.df <- rbind(outputTimeseries.df,data.frame(Fbar=Fbar[f], Mbar=Mbar[f], Run=b, Year=1:(tmax+1), SSB=SSB, Rec= Recruitment, Rec0=c(Recruitment[(recage+1):(tmax+1)],rep(NA,recage)), Ct= Ct))
  #getting all data over for all individual years
    } 
 
    #outputlogMsd.df<-rbind(outputlogMsd.df,data.frame(Run=b, Finput=Finput, Mbar=Mbar, Fbar = Fbar, Rec=Recmean, Catch=Cmean, SSB=SSBmean, Blim.any = Blim.any, Blim.mean = Blim.mean, Blim.prop= Blim.prop))
    
     
  } ## end bootstrap
  #gc()
  #outputlogMsd.df
  
  # save output, then remove last summary table - clears up and speeds up
  #saveRDS(outputlogMsd.df,file=paste0("temp/",multiplier,"_",sp,"_temp.rds"),compress="xz") #saves one file for each multiplier 
  #remove(outputlogMsd.df)
  print(multiplier) 
  
  timeseries <- outputTimeseries.df  #storing every value (e.g. SSB & Rec) for every timeseries
  
}  

##### Run simulation function  ######

source(here('Scripts_R/NSSH_Parameters.R'), echo=TRUE) #run parameter script nssh

# can use r=5 or 5.1 or 5.2 and BH describes deterministic, stochastic or AR stochastic Rickers and BH

# map function and save temp files
system.time({seq(0.2, 2, 0.01) %>% 
  map(~ am(multiplier =.x, r=6, sp="nssh", start.rec, Fsel, tmax, maturity, weights, boot=10, fseq=200, fmax=3))}) #Mvariation,

# for running the time series output
NSSHrec6 <- seq(0.2, 2, 0.2) %>% 
    map(~ am(multiplier =.x, r=6, sp="nssh", start.rec, Fsel, tmax, maturity, weights, boot=10, fseq=200, fmax=3))

# running catch at age series 
nssh_catchatage <- seq(0.2, 2, 0.2) %>% 
  map(~ am(multiplier =.x, r=5.1, sp="nssh", start.rec, Fsel, tmax, maturity, weights, boot=5, fseq=200, fmax=3))


#### collect temp files from simulation ####
nsshMrec6 <- list.files(pattern = "nssh_temp.rds", recursive = TRUE) %>% 
  map(readRDS) %>%   # 5, 5.1, 5.2 & 6 describes deterministic, stochastic or AR stochastic Rickers, & BH
  bind_rows() 


# getting catch, ssb and N matrices for a timeseries and a F to check runs      Year=1:(tmax+1)
out_df <- nsshMrec6.1 %>% group_by(Mbar,Fbar) %>% 
  summarise(Catch_m=cummean(Catch), SSB_m=cummean(SSB), Rec_m=cummean(Rec), N=cummax(Run))

out_df %>% filter(round(Fbar,4)==0.2223 & Mbar %in% sample(out_df$Mbar,7)) %>% 
  ggplot(aes(N, SSB_m, color=as.factor(Mbar))) + geom_line()


source(here('Scripts_R/REDFISH_parameters.R'), echo=TRUE) # parameter script redfish

# map function and save temp files
system.time({seq(0.2, 2, 0.01) %>% 
    map(~ am(multiplier =.x, r=5, sp="red", start.rec, Fsel, tmax, maturity, weights, boot=10, fseq=200,fmax=0.5))})

# for running the time series output
Mrec5_ts <- seq(0.2, 2, 0.2) %>% 
  map(~ am(multiplier =.x, r=5, sp="red", start.rec, Fsel, tmax, maturity, weights, boot=10, fseq=200, fmax=0.5))

# collect temp files
Mrec5 <- list.files(pattern = "red_temp.rds",recursive = TRUE) %>%
  map(readRDS) %>% 
  bind_rows()

#checking runs 
out_df <- Mrec5.2 %>% group_by(Mbar,Fbar) %>% summarise(Catch_m=cummean(Catch), SSB_m=cummean(SSB), Rec_m=cummean(Rec), N=cummax(Run))
out_df %>% filter(round(Fbar,4)==0.0100 & Mbar %in% sample(out_df$Mbar,7)) %>% 
  ggplot(aes(N,Rec_m,color=as.factor(Mbar))) + geom_line()


source(here('Scripts_R/COD_parameters.R'), echo=TRUE) #parameter script cod 

# map function and save temp files   r= 5 or 5.1 or 5.2
system.time({seq(0.2, 2, 0.01) %>% 
    map(~ am(multiplier =.x, r=5, sp="cod", start.rec, Fsel, tmax, maturity, weights, boot=10, fseq=200, fmax=3))})

# for running the time series output
codMrec6_ts <- seq(0.2, 2, 0.2) %>% 
  map(~ am(multiplier =.x, r=6, sp="cod", start.rec, Fsel, tmax, maturity, weights, boot=5, fseq=200, fmax=3))

#collecting temp file
codMrec5 <- list.files(pattern = "cod_temp.rds", recursive = TRUE) %>%
  map(readRDS) %>% 
  bind_rows()

#check runs
out_df <- codMrec5.1 %>% group_by(Mbar,Fbar) %>% summarise(Catch_m=cummean(Catch), SSB_m=cummean(SSB), Rec_m=cummean(Rec), N=cummax(Run))
out_df %>% filter(round(Fbar,4)==0.3041 & Mbar %in% sample(out_df$Mbar,7)) %>% ggplot(aes(N,Rec_m,color=as.factor(Mbar))) + geom_line()


# saving all the data & loading # 
write.csv(codMrec5.2, here('Scripts/codMrec5.2.csv'))
read.csv('C:/Users/jte084/OneDrive - University of Bergen/Fmsy/Sim_data_2/NSSHMrec5.1.csv')


####### Summarise for Fbar & Mbar  #######

susFbars <- function(dataset) # Summarise Mean/median SSB & Catch per Mbar and Fbar
{
  varsumBRF <- dataset %>% 
    group_by(Mbar,Fbar) %>%
    summarise(Catch=mean(Catch)/1000000, #divide to get million tons
              SSB=mean(SSB)/1000000,
              Risk.Blimany=sum(Blim.any)/n(),
              Risk.Blimmean=sum(Blim.mean)/n(),
              Risk.Blimprop=sum(Blim.prop)/n()) 
  
  varsumBRF <- varsumBRF
}


varnssh_BRFrec6 <- susFbars(nsshMrec6)
varnssh_BRFrec5$multiplier <- seq(0.2, 2, 0.01) %>% rep(each=200) #multiplier for relativity bt stocks

varcod_BRFrec5 <- susFbars(codMrec5)
varcod_BRFrec5$multiplier <- seq(0.2, 2, 0.01) %>% rep(each=200)

varredfish_BRFrec5 <- susFbars(Mrec5)
varredfish_BRFrec5$multiplier <- seq(0.2, 2, 0.01) %>% rep(each=200)

# instead of each line 
speciesRFs.list <- lapply(list(varnssh_BRFrec5.1, varcod_BRFrec5.1, varredfish_BRFrec5.1),FUN=RFs)


####### REFERENCE POINTS  #######

RFs <- function(dataset)   # raw numbers. FMSY and MSY, while having a risk of B<Blim below 5%)
{
  sumBRF <- dataset %>% 
    mutate(Catch.any= ifelse(Risk.Blimany<0.05,Catch,0), Catch.prop=ifelse(Risk.Blimprop<0.05,Catch,0), Catch.mean=ifelse(Risk.Blimany<0.05,Catch,0)) %>%  
    
    group_by(Mbar) %>%       
    summarise(
      Fmsy= Fbar[Catch==max(Catch)], 
      Fpa.any= mean(Fbar[Catch.any==max(Catch.any)]),   # most precautionary term
      Fpa.prop= mean(Fbar[Catch.prop==max(Catch.prop)]), # less precautionary term
      Fpa.mean= mean(Fbar[Catch.mean==max(Catch.mean)]), # the least precautionary term                                      
      
      MSY.any=max(Catch.any),0,
      MSY.prop=max(Catch.prop),0,
      MSY.mean=max(Catch.mean),0,
      
      Bmsy=SSB[Catch==max(Catch)],
      MSY=max(Catch),0) %>% 
    
    mutate(Fmsy=ifelse(MSY==0,0,Fmsy),Fpa.any=ifelse(MSY.any==0,0,Fpa.any),Fpa.prop=ifelse(MSY.prop==0,0,Fpa.prop),Fpa.mean=ifelse(MSY.mean==0,0,Fpa.mean))
  
  sumBRF <- sumBRF
}


nssh_RFs5 <- RFs(varnssh_BRFrec5)
nssh_RFs5$multiplier <- seq(0.2, 2, 0.01)  

cod_RFs5 <- RFs(varcod_BRFrec5)
cod_RFs5$multiplier <- seq(0.2, 2, 0.01)

redfish_RFs5 <- RFs(varredfish_BRFrec5)
redfish_RFs5$multiplier <- seq(0.2, 2, 0.01)

    
    

##### Relative numbers to the FMSY produced at ICES M #####

relativeCatch <- function(dataset1, dataset2, m){
  
  # catch
  defMSY <- dataset2$MSY[dataset2$multiplier== m]   
  defFMSY <- dataset2$Fmsy[dataset2$multiplier== m]
  dataset1 <- dataset1 %>% mutate(relCatch = Catch/defMSY)
  dataset1 <- dataset1 %>% mutate(relFbar = Fbar/defFMSY)
  
  # SSB
  defBmsy <- dataset2$Bmsy[dataset2$multiplier== m]   
  dataset1 <- dataset1 %>% mutate(relSSB = SSB/defBmsy)
  
  
  #pre. any
  defMSY.any <- dataset2$MSY.any[dataset2$multiplier== m]   
  defFpa.any <- dataset2$Fpa.any[dataset2$multiplier== m]
  dataset1 <- dataset1 %>% mutate(relCatch.any = Catch/defMSY.any)
  dataset1 <- dataset1 %>% mutate(relFbar.any = Fbar/defFpa.any)
  
  #pre. prop
  defMSY.prop <- dataset2$MSY.prop[dataset2$multiplier== m]   
  defFpa.prop <- dataset2$Fpa.prop[dataset2$multiplier== m]
  dataset1 <- dataset1 %>% mutate(relCatch.prop = Catch/defMSY.prop)
  dataset1 <- dataset1 %>% mutate(relFbar.prop = Fbar/defFpa.prop)
  

  #pre. mean
  defMSY.mean <- dataset2$MSY.mean[dataset2$multiplier== m]   
  defFpa.mean <- dataset2$Fpa.mean[dataset2$multiplier== m]
  dataset1 <- dataset1 %>% mutate(relCatch.mean = Catch/defMSY.mean)
  dataset1 <- dataset1 %>% mutate(relFbar.mean = Fbar/defFpa.mean)
  
  
  output <- dataset1
} # relative Catch and F


relnssh_BRFrec5 <- relativeCatch(varnssh_BRFrec5, nssh_RFs5, m=1)

relredfish_BRFrec5.2 <- relativeCatch(varredfish_BRFrec5.2, redfish_RFs5.2, m=1)

relcod_BRFrec6 <- relativeCatch(varcod_BRFrec6, cod_RFs6, m=1)



relRFs <- function(dataset)   # relative FMSY and MSY
{
  sumBRF <- dataset %>% 
    mutate(relCatch.any= ifelse(Risk.Blimany<0.05,relCatch,0), relCatch.prop=ifelse(Risk.Blimprop<0.05,relCatch,0), relCatch.mean=ifelse(Risk.Blimany<0.05,relCatch,0)) %>%  
    
    group_by(Mbar) %>%       
    summarise(
      relFmsy= relFbar[relCatch==max(relCatch)],
      relFpa.any= mean(relFbar[relCatch.any==max(relCatch.any)]),   
      relFpa.prop= mean(relFbar[relCatch.prop==max(relCatch.prop)]), 
      relFpa.mean= mean(relFbar[relCatch.mean==max(relCatch.mean)]),                                          
      
      
      MSY.any=max(relCatch.any),0,
      MSY.prop=max(relCatch.prop),0,
      MSY.mean=max(relCatch.mean),0,
      
      Bmsy=SSB[relCatch==max(relCatch)],
      relMSY=max(relCatch),0) %>% 
    
    mutate(relFmsy=ifelse(relMSY==0,0,relFmsy), relFpa.any=ifelse(MSY.any==0,0,relFpa.any), relFpa.prop=ifelse(MSY.prop==0,0,relFpa.prop),relFpa.mean=ifelse(MSY.mean==0,0,relFpa.mean))
  
  sumBRF <- sumBRF
}

relnssh_RFs5 <- relRFs(relnssh_BRFrec5)
relnssh_RFs5 <- relnssh_RFs5 %>% filter(round(Mbar, 5) %in% round(seq(0.0300,0.3000, 0.0300), 5))


relcod_RFs6 <- relRFs(relcod_BRFrec6)
relcod_RFs6 <- relcod_RFs6  %>% filter(round(Mbar, 5) %in% round(seq(0.040, 0.400, 0.040), 5))


relredfish_RFs5.2 <- relRFs(relredfish_BRFrec5.2)
relredfish_RFs5.2 <- relredfish_RFs5.2%>% filter(round(Mbar,5) %in% round(seq(0.01, 0.1, 0.01), 5))

