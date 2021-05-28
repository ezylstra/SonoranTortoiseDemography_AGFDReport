###################################################################################################################################
#Multi-state model for capture-recapture data
###################################################################################################################################

library(plyr)
library(reshape2)
library(jagsUI)
library(runjags)

rm(list=ls())

#-----------------------------------------------------------------------------------------------# 
# Load data
#-----------------------------------------------------------------------------------------------# 
dat <- read.csv('CapRecapData.csv',header=TRUE,stringsAsFactors=FALSE,na.strings=c('',NA))
survs <- read.csv('PlotSurveySchedule.csv',header=TRUE,stringsAsFactors=FALSE)
pdsi <- read.csv('PDSI.csv',header=TRUE,stringsAsFactors=FALSE)
precip.norms <- read.csv('Precip_norms.csv',header=TRUE,stringsAsFactors=FALSE)
precip <- read.csv('Precip_Monthly.csv',header=TRUE,stringsAsFactors=FALSE)
disttocity <- read.csv('PlotDistToCity.csv',header=TRUE,stringsAsFactors=FALSE)

#-----------------------------------------------------------------------------------------------# 
# Formatting survey schedule
#-----------------------------------------------------------------------------------------------# 
#This is an excel file sent to us in late 2020, with a summary of plot survey effort from 1987-2020
#In each column labeled yYYYY, values are the number of person-days spent surveying

  surv.l <- melt(survs,id.vars=c('plot','code','area.sqmi'),value.name='persondays')
  surv.l$yr <- as.numeric(substr(surv.l$variable,2,5))
  surv.l <- surv.l[!is.na(surv.l$persondays),c('plot','code','area.sqmi','yr','persondays')]

#Adding entry for Harcuvar Mtns survey in 2020. In DataQuestions.docx file, it's noted there were 62 person days
  surv.l <- rbind(surv.l,data.frame(plot='Harcuvar Mtns',code='HM',area.sqmi=surv.l$area.sqmi[surv.l$code=='HM'][1],yr=2020,persondays=62))

#Change name of Four Peaks plot to match that in live tortoise database:
  surv.l$plot[surv.l$plot=='FS Four Peaks'] <- 'Four Peaks'

#Number of years between consecutive surveys
  surv.l <- surv.l[with(surv.l,order(code,yr)),]
  surv.l$interval <- NA
  for(i in 2:nrow(surv.l)){
    surv.l$interval[i] <- ifelse(surv.l$code[i]!=surv.l$code[i-1],NA,surv.l$yr[i]-surv.l$yr[i-1])
  }
  # summary(surv.l$interval)

#Create a wide version of the survey data table, with 1/0 indicating when surveys were done
  surv.w <- dcast(surv.l, plot + code + area.sqmi ~ yr, value.var='persondays')
  surv.w.mat <- surv.w[,4:ncol(surv.w)]
  surv.w.mat[!is.na(surv.w.mat)] <- 1
  surv.w.mat[is.na(surv.w.mat)] <- 0
  names(surv.w.mat) <- paste0('y',names(surv.w.mat))
  #Add columns for years when no plots were surveyed (1989, 2009, 2011-2014)
  surv.w.mat <- cbind(surv.w.mat,data.frame(y1989=rep(0,nrow(surv.w.mat)),y2009=0,y2011=0,y2012=0,y2013=0,y2014=0))
  surv.w.mat <- surv.w.mat[,order(names(surv.w.mat))]
  surv.w <- cbind(plot=surv.w[,2],surv.w.mat)

#Summarizing survey effort by plot
  surv.p <- ddply(surv.l,.(plot,code),summarize,area.sqmi=mean(area.sqmi),n.surveys=sum(!is.na(persondays)),
                  first.surv=min(yr[!is.na(persondays)]),last.surv=max(yr[!is.na(persondays)]))

#Summarizing survey effort by year
  surv.y <- ddply(surv.l,.(yr),summarize,nplots=length(plot),eff.per.survey=sum(persondays)/length(persondays))
  surv.y <- rbind(surv.y,data.frame(yr=c(1989,2009,2011:2014),nplots=0,eff.per.survey=NA))
  surv.y <- surv.y[order(surv.y$yr),]
  surv.y$period <- ifelse(surv.y$yr %in% 1987:2008,'1','2')
  # ddply(surv.y,.(period),summarize,nplots.mn=mean(nplots), nplots.md=median(nplots),persondays.mn=mean(eff.per.survey,na.rm=T))

#-----------------------------------------------------------------------------------------------# 
# Create capture histories for each individual
#-----------------------------------------------------------------------------------------------# 
#Format date
  dat$obsdate <- as.Date(dat$obsdate,format='%Y-%m-%d')

#Fill in MCL, where needed, to determine adult/juvenile status (cutoff = 180mm)
  dat <- dat[with(dat,order(plot,tort,obsdate)),]
  count(dat$captype[is.na(dat$MCL)]) #Only a few instances where MCL is missing on first capture ever or on year
  dat[dat$tort %in% c(dat$tort[which(is.na(dat$MCL) & dat$captype!='SecondOfYr')]),]
  #HF-122, FP-149: first listed MCL is >180, so we know they're adults. 
  #FP-209, HM-220: never had MCL listed (so unclear if adults).  Remove.
  dat <- dat[!dat$tort %in% c('HM-220','FP-209'),]

  #Using previous measurement if MCL missing.
  for(i in 2:nrow(dat)){
    dat$MCL[i] <- ifelse(dat$tort[i]==dat$tort[i-1] & is.na(dat$MCL[i]),dat$MCL[i-1],dat$MCL[i])
  }

#Clean up sex assignments (1=female, 2=male, 3=juv/unk)
  sexes <- ddply(dat,.(plot,tort),summarize,n.sex=length(unique(sex)),sex.first=sex[1],sex.last=sex[length(sex)])
  #With no other information, I'll assume sex at last capture is correct
  dat$sex <- sexes$sex.last[match(dat$tort,sexes$tort)]
  #3630 fem, 3734 male, 739 unk
  dat[dat$MCL>=180 & dat$sex==3,]  #only 7 individuals that were captured as adults but not assigned a sex

#Retain only a single capture each year (using MCL at first capture)
  datyr <- ddply(dat,.(plot,tort,sex,yr),summarize,MCL=MCL[1])
  nrow(datyr) #4146 captures
  length(unique(datyr$tort)) #2035 individuals
  #Add stage (adult/juvenile) assignment
  datyr$stage <- ifelse(datyr$MCL<180,1,2)

  #Dataframe with number of tortoises caught at least once in each year a plot was surveyed
  plotyr <- ddply(datyr,.(plot,yr),summarize,ntorts=length(tort))
  
  #Dataframe with number of tortoises caught ever at each plot
  plottort <- ddply(datyr,.(plot),summarize,ntorts=length(unique(tort)))

#Create matrix with capture histories 
#1=captured as juvenile; 2=captured as adult; 3=plot surveyed but tortoise not captured; NA=plot not surveyed)  
  cr <- dcast(datyr, plot + tort + sex ~ yr, value.var='stage')
  cr.mat <- cr[,4:ncol(cr)]
  names(cr.mat) <- paste0('y',names(cr.mat))
  #Add columns for years when no plots were surveyed (1989, 2009, 2011-2014)
  cr.mat <- cbind(cr.mat,data.frame(y1989=rep(NA,nrow(cr.mat)),y2009=NA,y2011=NA,y2012=NA,y2013=NA,y2014=NA))
  cr.mat <- cr.mat[,order(names(cr.mat))]
  #Change NAs to 3 when the plot was surveyed, but tortoise wasn't captured
  #NAs indicate plot wasn't surveyed that year
  for(i in 1:nrow(cr)){
    years <- which(surv.w.mat[surv.w$plot==cr$plot[i],]==1)
    for(j in years){
      cr.mat[i,j] <- ifelse(is.na(cr.mat[i,j]),3,cr.mat[i,j])
    }
  }
  cr <- cbind(cr[,1:3],cr.mat)

  # #Check: number of tortoises captured each year at each plot (same totals from datyr and cr dataframes?)
  # crcheck <- melt(cr,id.vars=c('plot','tort','sex'))
  # names(crcheck)[4] <- 'yr'
  # crcheck$yr <- as.numeric(substr(crcheck$yr,2,5))
  # crcheck <- ddply(crcheck,.(plot,yr),summarize,ntorts=length(which(value!=3)))
  # crcheck <- crcheck[crcheck$ntorts!=0,]
  # all.equal(crcheck$ntorts,plotyr$ntorts)

  #First and last known state of each tortoise
  firstlast <- data.frame(cr$tort)
  firstlast$state1 <- apply(cr[,4:ncol(cr)],1,function(x) head(x[!is.na(x) & x!=3],1))
  firstlast$state2 <- apply(cr[,4:ncol(cr)],1,function(x) tail(x[!is.na(x) & x!=3],1))
  #Proportion of individuals first captured as juvenile/adult
  sum(firstlast$state1==1)/nrow(firstlast); sum(firstlast$state1==2)/nrow(firstlast)
  #Proportion of juveniles that were subsequently captured as adults
  sum(firstlast$state1==1 & firstlast$state2==2)/sum(firstlast$state1==1)

#-----------------------------------------------------------------------------------------------# 
# Function to create initial values for JAGS
#-----------------------------------------------------------------------------------------------# 

#Create vector indicating the first year each tortoise was caught as juvenile:
  first1 <- rep(NA,nrow(cr.mat))
  for(i in 1:nrow(cr.mat)){
    first1[i] <- (1:ncol(cr.mat))[!is.na(cr.mat[i,]) & cr.mat[i,]==1][1]
  }

#Create vector indicating the first year each tortoise was caught as adult:  
  first2 <- rep(NA,nrow(cr.mat))
  for(i in 1:nrow(cr.mat)){
    first2[i] <- (1:ncol(cr.mat))[!is.na(cr.mat[i,]) & cr.mat[i,]==2][1]
  }

#Create a matrix of initial values for latent states (z)
#NAs up to and including the first occasion, then 1's or 2's through the remainder
  ch.init <- function(y,f1,f2){
    for (i in 1:length(f1)){
      if(!is.na(f1[i]) & f1[i]==ncol(y)) {y[i,] <- NA} else
        if(is.na(f1[i]) & !is.na(f2[i]) & f2[i]==ncol(y)) {y[i,] <- NA} else
          if(is.na(f1[i]) & !is.na(f2[i]) & f2[i]!=ncol(y)) {y[i,1:f2[i]] <- NA; y[i,(f2[i]+1):ncol(y)] <- 2} else
            if(!is.na(f1[i]) & f1[i]!=ncol(y) & is.na(f2[i])) {y[i,1:f1[i]] <- NA; y[i,(f1[i]+1):ncol(y)] <- 1} else
              if(!is.na(f1[i]) & !is.na(f2[i]) & (f2[i]-f1[i]==1)) {y[i,1:f1[i]] <- NA; y[i,f2[i]:ncol(y)] <- 2} else
              {y[i,1:f1[i]] <- NA; y[i,(f1[i]+1):(f2[i]-1)] <- 1; y[i,f2[i]:ncol(y)] <- 2}}
    return(y)
  }

#-----------------------------------------------------------------------------------------------# 
# Create vector indicating the first year each tortoise was captured
#-----------------------------------------------------------------------------------------------# 
#Create vector indicating the first year each tortoise was caught:
  first <- rep(NA,nrow(cr.mat))
  for(i in 1:nrow(cr.mat)){
    first[i] <- (1:ncol(cr.mat))[!is.na(cr.mat[i,]) & cr.mat[i,]!=3][1]
  }

#-----------------------------------------------------------------------------------------------# 
# Formatting covariates
#-----------------------------------------------------------------------------------------------#
#Formatting individual covariates
  sex <- cr$sex
  sex[sex==3] <- NA
  male.ind <- sex-1

#Formatting site covariate: distance to nearest "major" city (population > 10,000)
  names(disttocity)[names(disttocity)=='code'] <- 'plot'
  disttocity <- disttocity[order(disttocity$plot),]
  # #Check that plots are in the same order they're found in capture histories
  # all.equal(unique(cr$plot),unique(disttocity$plot))
  dist.mn <- mean(disttocity$dist.km)
  dist.sd <- sd(disttocity$dist.km)
  distance <- (disttocity$dist.km - dist.mn)/dist.sd

#Formatting site covariate: precipitation normals
  precip.norms <- precip.norms[order(precip.norms$plot),]
  # #Check that plots are in the same order they're found in capture histories
  # all.equal(unique(cr$plot),unique(precip.norms$plot))
  precipnorm.mn <- mean(precip.norms$ppt.mm)
  precipnorm.sd <- sd(precip.norms$ppt.mm)
  precip.norm <- (precip.norms$ppt.mm - precipnorm.mn)/precipnorm.sd

#Formatting site*year covariate: drought
  #Going to calculate mean PDSI averaged over the previous 12 and 24 months (used 24-mon index in 2013 paper)
  pdsi <- pdsi[with(pdsi,order(div,yr,mon)),]
  pdsi$pdsi.12 <- NA
  pdsi$pdsi.24 <- NA
  jul.index <- which(pdsi$yr %in% 1988:2020 & pdsi$mon==7)
  for(i in jul.index){
    pdsi$pdsi.12[i] <- mean(pdsi$pdsi[(i-11):i])
    pdsi$pdsi.24[i] <- mean(pdsi$pdsi[(i-23):i])
  }
  pdsi.mn <- pdsi[!is.na(pdsi$pdsi.12),c('div','yr','pdsi.12','pdsi.24')]
  #Link plot to the appropriate climate division
  plotdiv <- data.frame(plot=unique(cr$plot),div=c(1,1,3,1,6,6,6,1,3,6,3,6,5,6,6,3,7),stringsAsFactors=FALSE)
  pdsi.plot <- expand.grid(yr=1988:2020,plot=plotdiv$plot)
  pdsi.plot$div <- plotdiv$div[match(pdsi.plot$plot,plotdiv$plot)]
  pdsi.plot <- join(pdsi.plot,pdsi.mn,by=c('div','yr'),type='left')
  pdsi.12 <- dcast(pdsi.plot,plot~yr,mean,value.var='pdsi.12')
  pdsi.24 <- dcast(pdsi.plot,plot~yr,mean,value.var='pdsi.24')
  # #Check that plots are in the same order they're found in capture histories
  # all.equal(unique(cr$plot),as.character(pdsi.12$plot))	
  # all.equal(unique(cr$plot),as.character(pdsi.24$plot))	
  pdsi.12.mat <- as.matrix(pdsi.12[,2:ncol(pdsi.12)])
  pdsi.24.mat <- as.matrix(pdsi.24[,2:ncol(pdsi.24)])
  pdsi12.mn <- mean(pdsi.12.mat)
  pdsi12.sd <- sd(pdsi.12.mat)
  pdsi24.mn <- mean(pdsi.24.mat)
  pdsi24.sd <- sd(pdsi.24.mat)	
  pdsi12.z <- (pdsi.12.mat - pdsi12.mn)/pdsi12.sd
  pdsi24.z <- (pdsi.24.mat - pdsi24.mn)/pdsi24.sd

#Formatting site*year covariate: precipitation	
  precip <- precip[with(precip,order(plot,season)),]
  # #Check that plots are in the same order they're found in capture histories
  # all.equal(unique(cr$plot),unique(precip$plot))
  #Cumulative precipitation (mm) at each plot from Aug-Jul
  precip.aj <- ddply(precip[precip$season %in% 1988:2020,],.(plot,season),summarize,ppt=sum(ppt))
  #Standardizing values by plot mean/SDs (can then compare effects of precipitation that's 1-SD above average)
  precipbyplot <- ddply(precip.aj,.(plot),summarize,ppt.mn=mean(ppt),ppt.sd=sd(ppt))
  precip.aj <- join(precip.aj,precipbyplot,by='plot',type='left')
  precip.aj$ppt.z <- (precip.aj$ppt - precip.aj$ppt.mn)/precip.aj$ppt.sd
  # #check:
  # ddply(precip.aj,.(plot),summarize,mn=mean(ppt.z),sd=sd(ppt.z))
  ppt.df <- dcast(precip.aj,plot~season,value.var='ppt.z')
  ppt.mat <- as.matrix(ppt.df[,2:ncol(ppt.df)])

#Formatting site*year covariate: survey effort (person days or person days scaled by plot area)
  surv.l$area.sqkm <- surv.l$area.sqmi*2.59
  surv.l$effort.sc <- surv.l$persondays/surv.l$area.sqkm
  #Standarize values
  surv.l$effort.z <- (surv.l$persondays - mean(surv.l$persondays))/sd(surv.l$persondays)
  surv.l$effort.sc.z <- (surv.l$effort.sc - mean(surv.l$effort.sc))/sd(surv.l$effort.sc)
  surv.l$plot <- NULL
  names(surv.l)[names(surv.l)=='code'] <- 'plot'
  surv.l <- surv.l[order(surv.l$plot),]	
  #Put into wide form (plot ~ yr)
  eff <- dcast(surv.l, plot ~ yr, value.var='effort.z')
  names(eff)[2:ncol(eff)] <- paste0('y',names(eff[,2:ncol(eff)]))
  effsc <- dcast(surv.l, plot ~ yr, value.var='effort.sc.z')
  names(effsc)[2:ncol(effsc)] <- paste0('y',names(effsc[,2:ncol(effsc)]))
  # #Check that plots are in the same order they're found in capture histories
  # all.equal(unique(cr$plot),unique(eff$plot))	
  #Add columns for years when no plots were surveyed (1989, 2009, 2011-2014)
  eff <- cbind(eff,data.frame(y1989=rep(NA,nrow(eff)),y2009=NA,y2011=NA,y2012=NA,y2013=NA,y2014=NA))
  eff <- eff[,order(names(eff))]
  effsc <- cbind(effsc,data.frame(y1989=rep(NA,nrow(eff)),y2009=NA,y2011=NA,y2012=NA,y2013=NA,y2014=NA))
  effsc <- effsc[,order(names(effsc))]
  #Convert from dataframe to matrix and replace all NAs with 0s
  #Don't need 1987 data, so removing
  eff <- as.matrix(eff[,3:ncol(eff)])
  eff[is.na(eff)] <- 0
  effsc <- as.matrix(effsc[,3:ncol(effsc)])
  effsc[is.na(effsc)] <- 0	

#Plot index for each tortoise
  plots <- data.frame(plot=unique(cr$plot),plot.index=1:length(unique(cr$plot)))
  plot.index <- plots$plot.index[match(cr$plot,plots$plot)]	

#-----------------------------------------------------------------------------------------------# 
# Run multistate model in JAGS
#-----------------------------------------------------------------------------------------------#	
#Prep data objects for JAGS
  ntorts <- nrow(cr.mat)             #number of tortoises
  nyears <- ncol(cr.mat)             #number of occasions
  nplots <- length(unique(cr$plot))  #number of plots
  yr.trend <- seq(0,nyears-2)        #sequence to evaluate a (logit) linear trend in adult survival
  yr.trendz <- (yr.trend - mean(yr.trend))/sd(yr.trend)
  yr.trend2z <- yr.trendz*yr.trendz

  tortdata <- list(y=as.matrix(cr.mat),
                   ntorts=ntorts,
                   nyears=nyears,
                   nplots=nplots,
                   first=first,
                   male=male.ind,
                   plot=plot.index,
                   distance=distance,
                   mean.precip=precip.norm,
                   drought=pdsi24.z,
                   yr.trend=yr.trendz,
                   yr.trend2=yr.trend2z,
                   precip=ppt.mat,
                   effort=eff)

#JAGS model
  # sink("MS_siteREtrend.txt")
  # cat("
  #   model{  
  # 
  #     #-- Priors and constraints
  # 
  #     alpha.p1 ~ dlogis(0,1)
  #     alpha.p2 ~ dlogis(0,1)
  #     beta.phi1 ~ dlogis(0,1)
  #     beta.phi2 ~ dlogis(0,1)
  #     gamma.psi ~ dlogis(0,1)
  # 
  #     a1.precip ~ dnorm(0,0.1)
  #     a1.effort ~ dnorm(0,0.1)
  #     a2.male ~ dnorm(0,0.1)
  #     a2.precip ~ dnorm(0,0.1)
  #     a2.effort ~ dnorm(0,0.1)
  # 
  #     b1.distance ~ dnorm(0,0.1)
  #     b1.mnprecip ~ dnorm(0,0.1)
  #     b1.drought ~ dnorm(0,0.1)
  #     b1.int ~ dnorm(0,0.1)
  #     b2.male ~ dnorm(0,0.1)
  #     b2.distance ~ dnorm(0,0.1)
  #     b2.mnprecip ~ dnorm(0,0.1)
  #     b2.drought ~ dnorm(0,0.1)
  #     b2.int ~ dnorm(0,0.1)
  #     b2.trend ~ dnorm(0,0.1)
  #     b2.trend2 ~ dnorm(0,0.1)
  # 
  #     c.mnprecip ~ dnorm(0,0.1)
  # 
  #     omega ~ dunif(0,1)
  # 
  #     sigma.site.2 ~ dt(0,pow(2.5,-2),1)T(0,)
  #     tau.site.2 <- 1/(sigma.site.2*sigma.site.2)
  # 
  #     for(p in 1:nplots){
  #       e.site.2[p] ~ dnorm(0,tau.site.2)
  #     }
  # 
  #     for (i in 1:ntorts){
  #       for(t in first[i]:(nyears-1)){
  # 
  #         logit(p1[i,t]) <- alpha.p1 + a1.precip*precip[plot[i],t] + a1.effort*effort[plot[i],t]
  #         logit(phi1[i,t]) <- beta.phi1 + b1.distance*distance[plot[i]] +
  #                             b1.mnprecip*mean.precip[plot[i]] + b1.drought*drought[plot[i],t] +
  #                             b1.int*mean.precip[plot[i]]*drought[plot[i],t]
  # 
  #         logit(psi12[i,t]) <- gamma.psi + c.mnprecip*mean.precip[plot[i]]
  # 
  #         logit(p2[i,t]) <- alpha.p2 + a2.male*male[i] + a2.precip*precip[plot[i],t] + a2.effort*effort[plot[i],t]
  #         logit(phi2[i,t]) <- beta.phi2 + b2.male*male[i] + b2.distance*distance[plot[i]] + b2.mnprecip*mean.precip[plot[i]] + 
  #                             b2.drought*drought[plot[i],t] + b2.int*mean.precip[plot[i]]*drought[plot[i],t] + 
  #                             b2.trend*yr.trend[t] + b2.trend2*yr.trend2[t] + e.site.2[plot[i]]
  # 
  #         #Define state transition probabilities
  #         #First index = states at time t-1, last index = states at time t
  #         ps[1,i,t,1] <- phi1[i,t] * (1-psi12[i,t])
  #         ps[1,i,t,2] <- phi1[i,t] * psi12[i,t]
  #         ps[1,i,t,3] <- 1-phi1[i,t]
  #         ps[2,i,t,1] <- 0
  #         ps[2,i,t,2] <- phi2[i,t]
  #         ps[2,i,t,3] <- 1-phi2[i,t]
  #         ps[3,i,t,1] <- 0
  #         ps[3,i,t,2] <- 0
  #         ps[3,i,t,3] <- 1
  # 
  #         #Define stage-dependent detection probabilities
  #         #First index = states at time t, last index = detection type at time t
  #         po[1,i,t,1] <- p1[i,t]
  #         po[1,i,t,2] <- 0
  #         po[1,i,t,3] <- 1-p1[i,t]
  #         po[2,i,t,1] <- 0
  #         po[2,i,t,2] <- p2[i,t]
  #         po[2,i,t,3] <- 1-p2[i,t]
  #         po[3,i,t,1] <- 0
  #         po[3,i,t,2] <- 0
  #         po[3,i,t,3] <- 1
  # 
  #       } #t
  #     } #i
  # 
  #     #-- Likelihood
  # 
  #     for(i in 1:ntorts){
  #       z[i,first[i]] <- y[i,first[i]]
  #       male[i] ~ dbern(omega)
  # 
  #       for (t in (first[i]+1):nyears){
  # 
  #         #State process: draw State[t] given State[t-1]
  #         z[i,t] ~ dcat(ps[z[i,t-1],i,t-1,])
  # 
  #         #Observation process: draw Obs[t] given State[t]
  #         y[i,t] ~ dcat(po[z[i,t],i,t-1,])
  # 
  #       } #t
  #     } #i
  # 
  #     #-- Derived parameters
  # 
  #     logit(psi12.mn) <- gamma.psi
  #     logit(phi1.mn) <- beta.phi1
  #     logit(p1.mn) <- alpha.p1
  # 
  #     phi2.f <- exp(beta.phi2)/(1 + exp(beta.phi2))
  #     phi2.m <- exp(beta.phi2 + b2.male)/(1 + exp(beta.phi2 + b2.male))
  #     p2.f <- exp(alpha.p2)/(1 + exp(alpha.p2))
  #     p2.m <- exp(alpha.p2 + a2.male)/(1 + exp(alpha.p2 + a2.male))
  # 
  #   } #model
  # ",fill=TRUE)
  # sink()

#MCMC settings, parameters, initial values  
  ni <- 15000; na <- 2000; nb <- 10000; nc <- 3; nt <- 15; ni.tot <- ni + nb

  params <- c('alpha.p1','a1.precip','a1.effort',
              'beta.phi1','b1.distance','b1.mnprecip','b1.drought','b1.int',
              'gamma.psi','c.mnprecip',
              'alpha.p2','a2.male','a2.precip','a2.effort',
              'beta.phi2','b2.male','b2.distance','b2.mnprecip','b2.drought','b2.int','b2.trend','b2.trend2',
              'omega','sigma.site.2','e.site.2',
              'p1.mn','phi1.mn','psi12.mn',
              'phi2.f','phi2.m','p2.f','p2.m')
  
  inits <- function() {list(alpha.p1=runif(1,-1,1),
                            alpha.p2=runif(1,-1,2),
                            beta.phi1=runif(1,-1,1),
                            beta.phi2=runif(1,1,3),
                            gamma.psi=runif(1,-2,0),
                            a1.precip=runif(1,-0.5,0.5),
                            a1.effort=runif(1,-0.5,0.5),
                            a2.male=runif(1,-0.5,0.5),
                            a2.precip=runif(1,-0.5,0.5),
                            a2.effort=runif(1,-0.5,0.5),
                            b1.distance=runif(1,-0.5,0.5),
                            b1.mnprecip=runif(1,-0.5,0.5),
                            b1.drought=runif(1,-0.5,0.5),
                            b1.int=runif(1,-0.5,0.5),
                            b2.male=runif(1,-0.5,0.5),
                            b2.distance=runif(1,-0.5,0.5),
                            b2.mnprecip=runif(1,-0.5,0.5),
                            b2.drought=runif(1,-0.5,0.5),
                            b2.int=runif(1,-0.5,0.5),
                            b2.trend=runif(1,-0.5,0.5),
                            b2.trend2=runif(1,-0.5,0.5),
                            c.mnprecip=runif(1,-0.5,0.5),
                            omega=runif(1,0,1),
                            sigma.site.2=runif(1,0,3),
                            male=ifelse(is.na(male.ind),1,NA),
                            z=ch.init(as.matrix(cr.mat),first1,first2))}

#Run model
  fit.ms <- jags(data=tortdata, inits=inits, parameters.to.save=params,
                 model.file='MS_siteREtrend.txt',
                 n.chains=nc, n.adapt=na, n.iter=ni.tot, n.burnin=nb, n.thin=nt,
                 parallel=T, n.cores=3, DIC=FALSE)
  
  #load(file.choose())
  print(fit.ms,digits=2)	

  #Create a matrix of posterior samples
  out <- fit.ms$samples
  comb <- combine.mcmc(out)
  phi1.s <- comb[,c('beta.phi1',colnames(comb)[grep('b1.',colnames(comb))])]
  phi2.s <- comb[,c('beta.phi2',colnames(comb)[grep('b2.',colnames(comb))])]
  phi2RE.s <- comb[,grep('e.site.2',colnames(comb))]
  psi12.s <- comb[,c('gamma.psi','c.mnprecip')]

#-----------------------------------------------------------------------------------------------# 
# Post-processing: covariate effects
#-----------------------------------------------------------------------------------------------#	

#For all figures, decide whether to use mean or median for measure of central tendency
  ctend <- mean
  #ctend <- median
#And 90% or 95% credible intervals
  #qprobs <- c(0.05,0.95)   #90% CIs
  qprobs <- c(0.025,0.975) #95% CIs

#Marginal effects of covariates on adult survival: mean precipitation and drought
  #Use survival estimates for last year (2019-2020)
  #Assume average distance from city (standardized distance = 0)
  phi2 <- phi2.s[,c('beta.phi2','b2.male','b2.mnprecip','b2.drought','b2.int','b2.trend','b2.trend2')]
  xmin <- min(pdsi24.z); xmax <- max(pdsi24.z)

  #Generate covariate values for figure
  #Use 3 values of mnprecip (values range from 177-409; mean = 280)
  mnprecip3 <- c(180,280,380)
  mnprecip3.z <- (mnprecip3-precipnorm.mn)/precipnorm.sd
  predx <- cbind(int=1,male=rep(0:1,each=300),mnprecip=rep(rep(mnprecip3.z,2),each=100),
                 drought=rep(seq(xmin,xmax,length=100),6))
  predx <- cbind(predx,interact=predx[,3]*predx[,4],trend=tail(yr.trendz,1),trend2=tail(yr.trend2z,1))
  predl <- predx %*% t(phi2)  #[600*7] %*% [7*3000] = [600*3000]
  pred <- exp(predl)/(1+exp(predl))  
  #Calculate mean/median and CRIs:  
  mean.f.dry <- apply(pred[1:100,],1,ctend)
  mean.f.avg <- apply(pred[101:200,],1,ctend)
  mean.f.wet <- apply(pred[201:300,],1,ctend)
  mean.m.dry <- apply(pred[301:400,],1,ctend)
  mean.m.avg <- apply(pred[401:500,],1,ctend)
  mean.m.wet <- apply(pred[501:600,],1,ctend)
  ci.f.dry <- apply(pred[1:100,],1,quantile,probs=qprobs)
  ci.f.avg <- apply(pred[101:200,],1,quantile,probs=qprobs)
  ci.f.wet <- apply(pred[201:300,],1,quantile,probs=qprobs)
  ci.m.dry <- apply(pred[301:400,],1,quantile,probs=qprobs)
  ci.m.avg <- apply(pred[401:500,],1,quantile,probs=qprobs)
  ci.m.wet <- apply(pred[501:600,],1,quantile,probs=qprobs)
  
  mycol <- col2rgb(c('salmon3','steelblue4'))
  col1 <- rgb(mycol[1,1],mycol[2,1],mycol[3,1],alpha=255,max=255)
  col1p <- rgb(mycol[1,1],mycol[2,1],mycol[3,1],alpha=0.2*255,max=255)
  col2 <- rgb(mycol[1,2],mycol[2,2],mycol[3,2],alpha=255,max=255)
  col2p <- rgb(mycol[1,2],mycol[2,2],mycol[3,2],alpha=0.2*255,max=255)

  plotx <- predx[1:100,'drought']*pdsi24.sd + pdsi24.mn

  #Figure with M/F adult survival at extreme mnprecip values (wet/dry plots)
  #jpeg('AdultSurvival_DroughtPrecipExtremes.jpg',width=80,height=70,units='mm',res=600)
    par(mar=c(2.5,3.5,0.5,0.6),cex=0.8)
    plot(mean.f.dry~plotx,type='l',lty=1,xaxt='n',yaxt='n',xlab='',ylab='', ylim=c(0.78,1),
         bty='n',yaxs='i',col=col1)
    axis(1,at=c(par('usr')[1],par('usr')[2]),tck=F,labels=F)
    axis(1,at=seq(-4,4,by=2),labels=seq(-4,4,by=2),tcl=-0.25,mgp=c(1.5,0.4,0))
    axis(2,at=c(par('usr')[3],par('usr')[4]),tck=F,labels=F)
    axis(2,at=seq(0.8,1,by=0.05),labels=c('0.80','0.85','0.90','0.95','1.00'),
         tcl=-0.25,las=1,mgp=c(1.5,0.5,0))
    #polygon(c(plotx,rev(plotx)),c(ci.f.dry[1,],rev(ci.f.dry[2,])),col=col1p,border=NA)
    lines(mean.f.wet~plotx,type='l',lty=2,col=col1)
    #polygon(c(plotx,rev(plotx)),c(ci.f.wet[1,],rev(ci.f.wet[2,])),col=col1p,border=NA)  
    lines(mean.m.dry~plotx,type='l',lty=1,col=col2)
    #polygon(c(plotx,rev(plotx)),c(ci.m.dry[1,],rev(ci.m.dry[2,])),col=col2p,border=NA)  
    lines(mean.m.wet~plotx,type='l',lty=2,col=col2)
    #polygon(c(plotx,rev(plotx)),c(ci.m.wet[1,],rev(ci.m.wet[2,])),col=col2p,border=NA)
    arrows(x0=0,x1=0,y0=0.68,y1=1,length=0,col='gray50',lty=3)
    mtext('Adult survival',side=2,las=0,line=2.5,cex=0.8)
    mtext('PDSI (24-month)',side=1,line=1.5,cex=0.8)
    legend('bottomright',c('Arid:F','Arid:M','Semiarid:F','Semiarid:M'),
           lty=c(1,1,2,2),col=c(col1,col2,col1,col2),bty='n')
  #dev.off()
  
  #Figure with M/F adult survival at plot with mnprecip = mean
  #jpeg('AdultSurvival_DroughtPrecipMean.jpg',width=80,height=70,units='mm',res=600)
    par(mar=c(2.5,3.5,0.5,0.6),cex=0.8)
    plot(mean.f.avg~plotx,type='l',lty=1,xaxt='n',yaxt='n',xlab='',ylab='', ylim=c(0.75,1),
         bty='n',yaxs='i',col=col1)
    axis(1,at=c(par('usr')[1],par('usr')[2]),tck=F,labels=F)
    axis(1,at=seq(-4,4,by=2),labels=seq(-4,4,by=2),tcl=-0.25,mgp=c(1.5,0.4,0))
    axis(2,at=c(par('usr')[3],par('usr')[4]),tck=F,labels=F)
    axis(2,at=seq(0.8,1,by=0.05),labels=c('0.80','0.85','0.90','0.95','1.00'),
         tcl=-0.25,las=1,mgp=c(1.5,0.5,0))
    polygon(c(plotx,rev(plotx)),c(ci.f.avg[1,],rev(ci.f.avg[2,])),col=col1p,border=NA)  
    lines(mean.m.avg~plotx,type='l',lty=1,col=col2)
    polygon(c(plotx,rev(plotx)),c(ci.m.avg[1,],rev(ci.m.avg[2,])),col=col2p,border=NA)  
    arrows(x0=0,x1=0,y0=0.68,y1=1,length=0,col='gray50',lty=3)
    mtext('Adult survival (95% CI)',side=2,las=0,line=2.5,cex=0.8)
    mtext('PDSI (24-month)',side=1,line=1.5,cex=0.8)
    legend('bottomright',c('Female','Male'),lty=1,col=c(col1,col2),bty='n')
  #dev.off()  

#Marginal effects of covariates on juvenile survival: mean precipitation and drought
  #Assume average distance from city (standardized distance = 0)
  phi1 <- phi1.s[,c('beta.phi1','b1.mnprecip','b1.drought','b1.int')]

  #Generate covariate values for figure
  pred1x <- cbind(int=1,mnprecip=rep(mnprecip3.z,each=100),drought=rep(seq(xmin,xmax,length=100),3))
  pred1x <- cbind(pred1x,interact=pred1x[,2]*pred1x[,3])
  pred1l <- pred1x %*% t(phi1)
  pred1 <- exp(pred1l)/(1+exp(pred1l))  
  #Calculate mean/median and CRIs:  
  mean.j.dry <- apply(pred1[1:100,],1,ctend)
  mean.j.avg <- apply(pred1[101:200,],1,ctend)
  mean.j.wet <- apply(pred1[201:300,],1,ctend)
  ci.j.dry <- apply(pred1[1:100,],1,quantile,probs=qprobs)
  ci.j.avg <- apply(pred1[101:200,],1,quantile,probs=qprobs)
  ci.j.wet <- apply(pred1[201:300,],1,quantile,probs=qprobs)
  
  #Figure with juvenile survival at extreme mnprecip values (wet/dry plots)
  #jpeg('JuvSurvival_DroughtPrecipExtremes.jpg',width=80,height=70,units='mm',res=600)
    par(mar=c(2.5,3.5,0.5,0.6),cex=0.8)
    plot(mean.j.dry~plotx,type='l',lty=1,xaxt='n',yaxt='n',xlab='',ylab='', ylim=c(0.45,1),
         bty='n',yaxs='i',col='black')
    axis(1,at=c(par('usr')[1],par('usr')[2]),tck=F,labels=F)
    axis(1,at=seq(-4,4,by=2),labels=seq(-4,4,by=2),tcl=-0.25,mgp=c(1.5,0.4,0))
    axis(2,at=c(par('usr')[3],par('usr')[4]),tck=F,labels=F)
    axis(2,at=seq(0.5,1,by=0.1),labels=c('0.50','0.60','0.70','0.80','0.90','1.00'),
         tcl=-0.25,las=1,mgp=c(1.5,0.5,0))
    #polygon(c(plotx,rev(plotx)),c(ci.j.dry[1,],rev(ci.j.dry[2,])),col=rgb(0,0,0,0.2),border=NA)
    lines(mean.j.wet~plotx,type='l',lty=2,col='black')
    #polygon(c(plotx,rev(plotx)),c(ci.j.wet[1,],rev(ci.j.wet[2,])),col=rgb(0,0,0,0.2),border=NA)  
    arrows(x0=0,x1=0,y0=0.45,y1=1,length=0,col='gray50',lty=3)
    mtext('Juvenile survival',side=2,las=0,line=2.5,cex=0.8)
    mtext('PDSI (24-month)',side=1,line=1.5,cex=0.8)
    legend('bottomright',c('Arid','Semiarid'),lty=c(1,2),col='black',bty='n')
  #dev.off() 
  
  #Figure with juvenile survival at plot with mnprecip = mean
  #jpeg('JuvSurvival_DroughtPrecipMean.jpg',width=80,height=70,units='mm',res=600)
    par(mar=c(2.5,3.5,0.5,0.6),cex=0.8)
    plot(mean.j.avg~plotx,type='l',lty=1,xaxt='n',yaxt='n',xlab='',ylab='', ylim=c(0.45,1),
         bty='n',yaxs='i',col='black')
    axis(1,at=c(par('usr')[1],par('usr')[2]),tck=F,labels=F)
    axis(1,at=seq(-4,4,by=2),labels=seq(-4,4,by=2),tcl=-0.25,mgp=c(1.5,0.4,0))
    axis(2,at=c(par('usr')[3],par('usr')[4]),tck=F,labels=F)
    axis(2,at=seq(0.5,1,by=0.1),labels=c('0.50','0.60','0.70','0.80','0.90','1.00'),
         tcl=-0.25,las=1,mgp=c(1.5,0.5,0))
    polygon(c(plotx,rev(plotx)),c(ci.j.avg[1,],rev(ci.j.avg[2,])),col=rgb(0,0,0,0.2),border=NA)  
    arrows(x0=0,x1=0,y0=0.45,y1=1,length=0,col='gray50',lty=3)
    mtext('Juvenile survival (95% CI)',side=2,las=0,line=2.5,cex=0.8)
    mtext('PDSI (24-month)',side=1,line=1.5,cex=0.8)
  #dev.off()  	

  #Stacked figure with juvenile and adult survival, extreme mnprecip values
  #jpeg('Survival_DroughtPrecipExtremes.jpg',width=80,height=135,units='mm',res=600)
    par(mfrow=c(2,1),mar=c(1.2,3.5,0.5,0.6),oma=c(1.3,0,0,0),cex=0.8)
    plot(mean.j.dry~plotx,type='l',lty=1,xaxt='n',yaxt='n',xlab='',ylab='', ylim=c(0.45,1),
         bty='n',yaxs='i',col='black')
    axis(1,at=c(par('usr')[1],par('usr')[2]),tck=F,labels=F)
    axis(1,at=seq(-4,4,by=2),labels=NA,tcl=-0.25,mgp=c(1.5,0.4,0))
    axis(2,at=c(par('usr')[3],par('usr')[4]),tck=F,labels=F)
    axis(2,at=seq(0.5,1,by=0.1),labels=c('0.50','0.60','0.70','0.80','0.90','1.00'),
         tcl=-0.25,las=1,mgp=c(1.5,0.5,0))
    #polygon(c(plotx,rev(plotx)),c(ci.j.dry[1,],rev(ci.j.dry[2,])),col=rgb(0,0,0,0.2),border=NA)
    lines(mean.j.wet~plotx,type='l',lty=2,col='black')
    #polygon(c(plotx,rev(plotx)),c(ci.j.wet[1,],rev(ci.j.wet[2,])),col=rgb(0,0,0,0.2),border=NA)  
    arrows(x0=0,x1=0,y0=0.45,y1=1,length=0,col='gray50',lty=3)
    mtext('Juvenile survival',side=2,las=0,line=2.5,cex=0.8)
    legend('bottomright',c('Arid','Semiarid'),lty=c(1,2),col='black',bty='n') 
    plot(mean.f.dry~plotx,type='l',lty=1,xaxt='n',yaxt='n',xlab='',ylab='', ylim=c(0.78,1),
         bty='n',yaxs='i',col=col1)
    axis(1,at=c(par('usr')[1],par('usr')[2]),tck=F,labels=F)
    axis(1,at=seq(-4,4,by=2),labels=seq(-4,4,by=2),tcl=-0.25,mgp=c(1.5,0.4,0))
    axis(2,at=c(par('usr')[3],par('usr')[4]),tck=F,labels=F)
    axis(2,at=seq(0.8,1,by=0.05),labels=c('0.80','0.85','0.90','0.95','1.00'),
         tcl=-0.25,las=1,mgp=c(1.5,0.5,0))
    #polygon(c(plotx,rev(plotx)),c(ci.f.dry[1,],rev(ci.f.dry[2,])),col=col1p,border=NA)
    lines(mean.f.wet~plotx,type='l',lty=2,col=col1)
    #polygon(c(plotx,rev(plotx)),c(ci.f.wet[1,],rev(ci.f.wet[2,])),col=col1p,border=NA)  
    lines(mean.m.dry~plotx,type='l',lty=1,col=col2)
    #polygon(c(plotx,rev(plotx)),c(ci.m.dry[1,],rev(ci.m.dry[2,])),col=col2p,border=NA)  
    lines(mean.m.wet~plotx,type='l',lty=2,col=col2)
    #polygon(c(plotx,rev(plotx)),c(ci.m.wet[1,],rev(ci.m.wet[2,])),col=col2p,border=NA)
    arrows(x0=0,x1=0,y0=0.68,y1=1,length=0,col='gray50',lty=3)
    mtext('Adult survival',side=2,las=0,line=2.5,cex=0.8)
    mtext('PDSI (24-month)',side=1,line=1.5,cex=0.8)
    legend('bottomright',c('Arid:F','Arid:M','Semiarid:F','Semiarid:M'),
           lty=c(1,1,2,2),col=c(col1,col2,col1,col2),bty='n')    
  #dev.off()  

#Temporal trends in adult survival
  #For an average site with PDSI = 0
  phi2t <- phi2.s[,c('beta.phi2','b2.male','b2.trend','b2.trend2')]
  
  predtx <- cbind(int=1,male=rep(0:1,each=33),trend=rep(yr.trendz,2),trend2=rep(yr.trend2z,2))
  predtl <- predtx %*% t(phi2t)
  predt <- exp(predtl)/(1+exp(predtl))  
  
  predt.both <- data.frame(endyr=1988:2020,interval=rep(paste(1987:2019,1988:2020,sep='-')))
  predt.both$female <- round(apply(predt[1:33,],1,ctend),3)
  predt.both$female.lcl <- round(apply(predt[1:33,],1,quantile,probs=qprobs[1]),3)
  predt.both$female.ucl <- round(apply(predt[1:33,],1,quantile,probs=qprobs[2]),3)
  predt.both$male <- round(apply(predt[34:66,],1,ctend),3)
  predt.both$male.lcl <- round(apply(predt[34:66,],1,quantile,probs=qprobs[1]),3)
  predt.both$male.ucl <- round(apply(predt[34:66,],1,quantile,probs=qprobs[2]),3)  
  
  #Figure with survival rates by year and sex
  #jpeg('AdultSurvival_Trend.jpg',width=80,height=70,units='mm',res=600)
    par(mar=c(2.5,3.5,0.5,0.6),cex=0.8)  
    plot(female~endyr,predt.both,type='l',lty=1,col=col1,ylim=c(0.73,1),xaxt='n',yaxt='n',
         bty='n',xlab='',ylab='',yaxs='i')
    axis(1,at=c(par('usr')[1],par('usr')[2]),tck=F,labels=F)
    axis(1,at=seq(1990,2020,by=5),labels=seq(1990,2020,by=5),tcl=-0.25,mgp=c(1.5,0.4,0))
    axis(2,at=c(par('usr')[3],par('usr')[4]),tck=F,labels=F)
    axis(2,at=seq(0.75,1,by=0.05),labels=c('0.75','0.80','0.85','0.90','0.95','1.00'),
         tcl=-0.25,las=1,mgp=c(1.5,0.5,0)) 
    polygon(c(predt.both$endyr,rev(predt.both$endyr)),c(predt.both$female.lcl,rev(predt.both$female.ucl)),
            border=NA,col=col1p)
    lines(male~endyr,predt.both,type='l',lty=1,col=col2)
    polygon(c(predt.both$endyr,rev(predt.both$endyr)),c(predt.both$male.lcl,rev(predt.both$male.ucl)),
            border=NA,col=col2p) 
    mtext('Adult survival (95% CI)',side=2,las=0,line=2.5,cex=0.8)
    mtext('Year',side=1,line=1.5,cex=0.8)
    legend('bottomright',c('Female','Male'),lty=1,col=c(col1,col2),bty='n') 
  #dev.off()  

#Temporal trend in PDSI values?
  pdsi24t <- unique(pdsi.plot[,c('yr','div','pdsi.24')])
  pdsi24t$yr0 <- pdsi24t$yr - min(pdsi24t$yr) 

  #Quick run of ML linear regression models
  summary(lm.trend1 <- lm(pdsi.24~yr0,data=pdsi24t))  
  summary(lm.trend5 <- lm(pdsi.24~yr0*factor(div),data=pdsi24t))
  AIC(lm.trend1); AIC(lm.trend5) #AIC 10 points lower for the simpler model
  
  #Bayesian linear regression
  # sink("LinearRegression.txt")
  # cat("
  #   model{
  # 
  #     for(i in 1:nobs){
  #       y[i] ~ dnorm(mu[i], tau)
  #       mu[i] <- b0 + b1*x[i]
  #     }
  # 
  #     b0 ~ dnorm(0,0.0000001)
  #     b1 ~ dnorm(0,0.0000001)
  #     sigma ~ dunif(0,1000)
  #     tau <- 1/(sigma*sigma)
  # 
  #   } #model
  # ",fill=TRUE)
  # sink()
  
  ni <- 2000; na <- 1000; nb <- 8000; nc <- 3; ni.tot <- ni + nb
  params <- c('b0','b1','sigma')
  inits <- function() {list(b0=runif(1,-30,0),
                            b1=runif(1,-2,2),
                            sigma=runif(1,1,100))}  
  pdsi.dat <- list(x=pdsi24t$yr0,
                   y=pdsi24t$pdsi.24,
                   nobs=nrow(pdsi24t))
  set.seed(123)
  fit.pdsi <- jags(data=pdsi.dat, inits=inits, parameters.to.save=params,
                   model.file='LinearRegression.txt',
                   n.chains=nc, n.adapt=na, n.iter=ni.tot, n.burnin=nb,
                   parallel=T, n.cores=3, DIC=F)
  print(fit.pdsi)
  
  X <- data.frame(int=1,yr0=seq(0,32,length=100))
  fit.s <- fit.pdsi$samples
  fit.mat <- combine.mcmc(fit.s)
  betas <- fit.mat[,1:2]
  pdsipreds <- as.matrix(X) %*% t(betas)
  cent.pdsi <- apply(pdsipreds,1,ctend)
  cri.pdsi <- apply(pdsipreds,1,quantile,probs=qprobs)
  
  col5 <- col2rgb(c('mediumpurple4','steelblue4','darkseagreen4','goldenrod4','salmon4'))
  trans <- 0.5
  col5.1p <- rgb(col5[1,1],col5[2,1],col5[3,1],alpha=trans*255,max=255)
  col5.2p <- rgb(col5[1,2],col5[2,2],col5[3,2],alpha=trans*255,max=255)
  col5.3p <- rgb(col5[1,3],col5[2,3],col5[3,3],alpha=trans*255,max=255)
  col5.4p <- rgb(col5[1,4],col5[2,4],col5[3,4],alpha=trans*255,max=255)
  col5.5p <- rgb(col5[1,5],col5[2,5],col5[3,5],alpha=trans*255,max=255)
  
  #Figure with juvenile survival at plot with mnprecip = mean
  #jpeg('PDSI_trend.jpg',width=6.5,height=3,units='in',res=600)
    par(mar=c(2.5,3.0,0.5,0.6),cex=0.7)
    plot(pdsi.24~yr0,data=pdsi24t[pdsi24t$div==1,],type='b',pch=19,xaxt='n',yaxt='n',xlab='',
         ylab='', ylim=c(-4.5,4.8),bty='n',yaxs='i',col=col5.1p)
    axis(1,at=c(par('usr')[1],par('usr')[2]),tck=F,labels=F)
    axis(1,at=seq(0,32,by=8),labels=seq(1988,2020,by=8),tcl=-0.25,mgp=c(1.5,0.4,0))
    axis(2,at=c(par('usr')[3],par('usr')[4]),tck=F,labels=F)
    axis(2,at=seq(-4,4,by=2),labels=seq(-4,4,by=2),tcl=-0.25,las=1,mgp=c(1.5,0.5,0))
    #arrows(x0=par('usr')[1],x1=par('usr')[2],y0=0,y1=0,length=0,col='gray50',lty=2)
    points(pdsi.24~yr0,data=pdsi24t[pdsi24t$div==3,],type='b',pch=19,col=col5.2p)
    points(pdsi.24~yr0,data=pdsi24t[pdsi24t$div==5,],type='b',pch=19,col=col5.3p)
    points(pdsi.24~yr0,data=pdsi24t[pdsi24t$div==6,],type='b',pch=19,col=col5.4p)
    points(pdsi.24~yr0,data=pdsi24t[pdsi24t$div==7,],type='b',pch=19,col=col5.5p)    
    lines(cent.pdsi~X$yr0)
    polygon(c(X$yr0,rev(X$yr0)),c(cri.pdsi[1,],rev(cri.pdsi[2,])),col=rgb(0,0,0,0.2),border=NA)
    mtext('PDSI (24-month)',side=2,las=0,line=2.1,cex=0.7)
    mtext('Year',side=1,line=1.5,cex=0.7)
    legend(x=27,y=4.8,c('1','3','5','6','7'),title='Division',pch=19,y.intersp=0.9,adj=c(0,0.5),
           col=c(col5.1p,col5.2p,col5.3p,col5.4p,col5.5p),bty='n')
  #dev.off() 

#-----------------------------------------------------------------------------------------------# 
# Post-processing: plot-specific estimates
#-----------------------------------------------------------------------------------------------#	

#Adult survival
  #Generating values for 2019-2020, with drought = 0
  phi2p <- phi2.s[,c('beta.phi2','b2.male','b2.distance','b2.mnprecip','b2.trend','b2.trend2')]
  
  predpx <- cbind(int=1,male=rep(0:1,each=17),distance=rep(distance,2),mnprecip=rep(precip.norm,2),
                  trend=tail(yr.trendz,1),trend2=tail(yr.trend2z,1))
  predpl <- predpx %*% t(phi2p)
  RE.mat <- rbind(t(phi2RE.s),t(phi2RE.s))
  predpl <- predpl + RE.mat
  predp <- exp(predpl)/(1+exp(predpl)) 	
  
  plotests <- data.frame(plot=disttocity$plot)
  plotests$ad.fem <- round(apply(predp[1:17,],1,ctend),2)
  plotests$ad.fem.lcl <- round(apply(predp[1:17,],1,quantile,qprobs[1]),2)
  plotests$ad.fem.ucl <- round(apply(predp[1:17,],1,quantile,qprobs[2]),2)
  plotests$ad.male <- round(apply(predp[18:34,],1,ctend),2)
  plotests$ad.male.lcl <- round(apply(predp[18:34,],1,quantile,qprobs[1]),2)
  plotests$ad.male.ucl <- round(apply(predp[18:34,],1,quantile,qprobs[2]),2)
  
  #Adult survival averaged over plots for 2019-2020
  phi2pA <- phi2.s[,c('beta.phi2','b2.male','b2.trend','b2.trend2')]
  
  predpxA <- cbind(int=1,male=0:1,trend=tail(yr.trendz,1),trend2=tail(yr.trend2z,1))
  predplA <- predpxA %*% t(phi2pA)
  predpA <- exp(predplA)/(1+exp(predplA)) 
  adsurvivalests <- data.frame(sex=c('F','M'))
  adsurvivalests$mn <- round(apply(predpA,1,ctend),2)
  adsurvivalests$lcl <- round(apply(predpA,1,quantile,qprobs[1]),2)
  adsurvivalests$ucl <- round(apply(predpA,1,quantile,qprobs[2]),2)
  adsurvivalests

#Juvenile survival
  phi1p <- phi1.s[,c('beta.phi1','b1.distance','b1.mnprecip')]
  predp1x <- cbind(int=1,distance=distance,mnprecip=precip.norm)
  predp1l <- predp1x %*% t(phi1p)
  predp1 <- exp(predp1l)/(1+exp(predp1l)) 
  plotests$juv <- round(apply(predp1,1,ctend),2)
  plotests$juv.lcl <- round(apply(predp1,1,quantile,qprobs[1]),2)
  plotests$juv.ucl <- round(apply(predp1,1,quantile,qprobs[2]),2)

#Transition rates
  predpsix <- cbind(int=1,mnprecip=precip.norm)
  predpsil <- predpsix %*% t(psi12.s)
  predpsi <- exp(predpsil)/(1+exp(predpsil))
  plotests$transition <- round(apply(predpsi,1,ctend),2)
  plotests$transition.lcl <- round(apply(predpsi,1,quantile,qprobs[1]),2)
  plotests$transition.ucl <- round(apply(predpsi,1,quantile,qprobs[2]),2)

plotests

#-----------------------------------------------------------------------------------------------# 
# Post-processing: population growth rates
#-----------------------------------------------------------------------------------------------#	 

#Gather plot-specific covariate values (standardized)
  plotcovs <- data.frame(plot=disttocity$plot,dist.z=distance,mnprecip.z=precip.norm)
  #For each plot, using 3 drought values (-3, 0, +3)
  phi.Xdf <- plotcovs[rep(seq_len(nrow(plotcovs)),each=3),] 
  drought3 <- rep(c(-3,0,3),nrow(plotcovs))
  phi.Xdf$drought.z <- (drought3 - pdsi24.mn)/pdsi24.sd
  phi.Xdf$int.z <- phi.Xdf$mnprecip.z*phi.Xdf$drought.z
  
  #Create a matrix of covariate values for juvenile survival  
  phi1.X <- as.matrix(phi.Xdf[,c('dist.z','mnprecip.z','drought.z','int.z')])
  phi1.X <- cbind(rep(1,nrow(phi1.X)),phi1.X)  
  
  #Create a matrix of covariate values for adult survival 
  #Last two columns are for the trend effect (using 2019-2020 estimates)
  phi2.X <- as.matrix(phi.Xdf[,c('dist.z','mnprecip.z','drought.z','int.z')])
  phi2.X <- cbind(rep(1,nrow(phi2.X)),phi2.X,rep(tail(yr.trendz,1),nrow(phi2.X)),rep(tail(yr.trend2z,1),nrow(phi2.X)))
  
  #Create a matrix of covariate values for transition rate  
  psi12.X <- as.matrix(phi.Xdf[,'mnprecip.z'])
  psi12.X <- cbind(rep(1,nrow(psi12.X)),psi12.X)    

#Calculate survival, transition estimates for each combination of covariates, iteration

  #Juvenile survival
  lphi1 <- phi1.X %*% t(phi1.s)
  phi1 <- exp(lphi1)/(1+exp(lphi1))
  
  #Adult survival (need to add in random site effects)
  phi2.s.female <- phi2.s[,colnames(phi2.s)!='b2.male']
  lphi2 <- phi2.X %*% t(phi2.s.female)
  REs <- t(phi2RE.s)
  REs <- REs[rep(seq_len(nrow(REs)),each=3),]
  lphi2RE <- lphi2 + REs
  phi2 <- exp(lphi2RE)/(1+exp(lphi2RE))	
  
  #Transition rates
  lpsi12 <- psi12.X %*% t(psi12.s)
  psi12 <- exp(lpsi12)/(1+exp(lpsi12))

#Create population projection matrices, and estimate lambda values
#Assuming recruitment = 0.32 f/f/yr 

  lambda <- matrix(NA,nrow=nrow(phi2),ncol=ncol(phi2))
  
  for(i in 1:nrow(phi2)){ 
    for (j in 1:ncol(phi2)){
      proj.mat <- matrix(c(phi1[i,j]*(1-psi12[j]), 0.32,
                           phi1[i,j]*psi12[j], phi2[i,j]),
                         nrow=2,ncol=2,byrow=TRUE)
      lambda[i,j] <- eigen(proj.mat)$values[1]
    }
  } 

#Summarize distributions of lambda values for each plot-drought combination
  lambda.df <- data.frame(plot=phi.Xdf$plot,drought=drought3)
  lambda.df$mn <- apply(lambda,1,mean)
  lambda.df$q0.025 <- apply(lambda,1,quantile,0.025)
  lambda.df$q0.05 <- apply(lambda,1,quantile,0.05)
  lambda.df$q0.5 <- apply(lambda,1,quantile,0.5)
  lambda.df$q0.95 <- apply(lambda,1,quantile,0.95)
  lambda.df$q0.975 <- apply(lambda,1,quantile,0.975)
  lambda.df$probdecline <- apply(lambda,1,function(x) sum(x<1)/length(x))
  (lambda.df <- lambda.df[with(lambda.df,order(drought,plot)),])
  ddply(lambda.df,.(drought),summarize,mn.prob=mean(probdecline),min.prob=min(probdecline),max.prob=max(probdecline))
  #write.table(lambda.df,'clipboard',sep='\t',row.names=FALSE,col.names=TRUE)

#Calculate overall lambda values (city, mnprecip = 0):
  #Juvenile survival
  phi1O.X <- phi1.X[1:3,c(1,4)]
  phi1O.s <- phi1.s[,c('beta.phi1','b1.drought')]
  lphi1O <- as.matrix(phi1O.X) %*% t(phi1O.s)
  phi1O <- exp(lphi1O)/(1+exp(lphi1O))
  #Adult female survival
  phi2O.X <- phi2.X[1:3,c(1,4,6,7)]
  phi2O.s <- phi2.s.female[,c('beta.phi2','b2.drought','b2.trend','b2.trend2')]
  lphi2O <- as.matrix(phi2O.X) %*% t(phi2O.s)
  phi2O <- exp(lphi2O)/(1+exp(lphi2O))
  #Transition
  psiO <- comb[,'psi12.mn']
  
  lambdaO <- matrix(NA,nrow=nrow(phi2O),ncol=ncol(phi2O))
  for(i in 1:nrow(phi2O)){ 
    for (j in 1:ncol(phi2O)){
      proj.mat <- matrix(c(phi1O[i,j]*(1-psiO[j]), 0.32,
                           phi1O[i,j]*psiO[j], phi2O[i,j]),
                         nrow=2,ncol=2,byrow=TRUE)
      lambdaO[i,j] <- eigen(proj.mat)$values[1]
    }
  }   
  lambdaO.df <- data.frame(drought=drought3[1:3])
  lambdaO.df$mn <- apply(lambdaO,1,mean)
  lambdaO.df$q0.025 <- apply(lambdaO,1,quantile,0.025)
  lambdaO.df$q0.05 <- apply(lambdaO,1,quantile,0.05)
  lambdaO.df$q0.5 <- apply(lambdaO,1,quantile,0.5)
  lambdaO.df$q0.95 <- apply(lambdaO,1,quantile,0.95)
  lambdaO.df$q0.975 <- apply(lambdaO,1,quantile,0.975)
  lambdaO.df$probdecline <- apply(lambdaO,1,function(x) sum(x<1)/length(x))
  lambdaO.df

#Correlations between lambda values and latitude/longitude
  lambda.df <- join(lambda.df,disttocity[,c('plot','lat','long')],by='plot',type='left')
  cor.test(lambda.df$mn[lambda.df$drought==-3],lambda.df$lat[lambda.df$drought==-3])
  cor.test(lambda.df$mn[lambda.df$drought==-3],lambda.df$long[lambda.df$drought==-3])
