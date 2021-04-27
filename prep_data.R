library(plyr)
kb       <- 8.62e-5
Tref     <- 15
T0       <- 273.15
T.K      <- function(x) 1/((Tref+T0)*kb) - 1/((x + T0)*kb)

#Get species name
get_spname <- function(string) {
  tx <- gregexpr(" ", string)
  outstring <- character(length(string))
  for (i in 1:length(string)) {
    tmp <- tx[[i]]
    NCH <- nchar(string[i])
    if (length(tmp) > 1) {
      nn <- tmp[2] - 1
      outstring[i] <- substr(string[i], 1, nn)
    } else if (tmp[1] == NCH) {
      outstring[i] <- as.character(NA)
    } else{
      outstring[i] <- string[i]
    }
  }
  return(outstring)
}

culldata <- function(newdat){
    newdat$X     <- T.K(newdat$Temperature)
    newdat$xbar  <- NA  #mean X (below optimal temperature)
    newdat$incl  <- NA  #Whether to include this in nonlinear fitting
    newdat$right <- NA  #Whether the data points fall to the right of Topt
    newdat$XoptL <- NA  #Optimal temperature without nonlinear fitting
    newdat$umL   <- NA  #Maximal growth rate without nonlinear fitting
    
    PHYu <- ddply(newdat, .(ID), summarize,
                         ID       = ID[1],
                     #    Habitat  = Habitat[1],
                         Group    = Group[1],
                     #    Genus    = Genus[1],
                         Species  = Species[1]
                     #    lnVolume = mean(log(Volume), na.rm=T)
                  )
    
    N  <- nrow(PHYu)

    for (i in 1:N){
      code <- PHYu$ID[i]
      tmp  <- newdat[newdat$ID == code,]
      tmp  <- tmp[order(tmp$Temperature,tmp$Growth),]
      
      #Whether include the data in later analysis
      incl <- TRUE
      
      # At least 3 data points 
      if (length(tmp[tmp$Growth > 0, ]$Temperature) < 4) incl <- FALSE
      
      # At least 2 unique temperature
      if (length(unique(tmp[tmp$Growth > 0, ]$Temperature)) < 2) incl <- FALSE
      
     #  Average Optimal temperature:
      tmp$r.temp <- round(tmp$Temperature, 1)
      tmp.mean   <- tapply(tmp$Growth, tmp$r.temp, mean)
      
      #Distinct temperatures
      dtemp <- attributes(tmp.mean)$dimnames[[1]]
      dtemp <- as.numeric(dtemp)
      optT  <- dtemp[which.max(tmp.mean)] #Optimal temperature
      newdat[newdat$ID == code, ]$XoptL <- T.K(optT)
      
      #Maximal growth rate at optimal temp
      newdat[newdat$ID == code, ]$umL   <- as.numeric(tmp.mean[which.max(tmp.mean)])
      
     # temperatures lower than optT
      Tl   <- tmp$Temperature[tmp$Temperature <= optT & tmp$Growth > 0]
    
      #At least 2 data points below or equal to optT    
      if (length(unique(Tl)) < 2)  incl <- FALSE
      
      # temperatures higher than optT
      Th   <- tmp$Temperature[tmp$Temperature > optT]
      #if (length(Th) < 1)  incl <- FALSE
      
      #Xbar for each taxon
      xbar <- mean(tmp[tmp$Temperature <= optT, 'X'])
      newdat[newdat$ID == code,]$xbar  <- rep(xbar, nrow(tmp))
      
      newdat[newdat$ID == code, ]$incl <- incl
      RIGHT                            <- rep(F, nrow(tmp))
      RIGHT[tmp$Temperature > optT]    <- T
      newdat[newdat$ID == code,]$right <- RIGHT
    }
    newdat$Sp    <- as.character(get_spname(newdat$Species))
    newdat       <- subset(newdat, incl ==TRUE)
    newdat2      <- subset(newdat, right==FALSE)
    return(list(oridat = newdat, newdat = newdat2))
}

#Original phytoplankton data
#source('Phy_merge.R')
load('NEWphy24Mar2020.Rdata')
phydat       <- culldata(newdat)

#Final oridat of phyto used for analysis
phydat2      <- phydat$newdat  #Contains zero growth rate
phydat       <- phydat$oridat

#Microzooplankton data from David
dat <- read.csv('Microzoo_23Mar2020.csv')
dat[dat$Habitat == 'freshwater','Habitat'] <- 'Freshwater'
zoodat         <- culldata(dat)
zoodat2        <- zoodat$newdat  #Final oridat of zoo used for analysis

#The microzoo data in Chen & Laws 2017
dat2017 <- read.csv('../Global_PP/microz.csv')
names(dat2017)[names(dat2017) == 'u']      <- 'Growth' 
names(dat2017)[names(dat2017) == 'Temp']   <- 'Temperature'
names(dat2017)[names(dat2017) == 'Classification'] <- 'Group'
names(dat2017)[names(dat2017) == 'Strain'] <- 'ID' 
dat2017 <- subset(dat2017, Habitat == 'Marine')
dat2017 <- culldata(dat2017)
dat2017 <- dat2017$newdat

#The Algae data in Wang et al. 2019
wang_Algae2019 <- read.csv('Wang2019Algae.csv')
wang_Algae2019$ID <- factor(wang_Algae2019$ID)
wang_Algae2019$Temperature <- wang_Algae2019$Temperature - T0 
wangAlgae <- culldata(wang_Algae2019)
wangAlgae <- wangAlgae$newdat

#Using nls following Wang et al. (2019)
wang_Algae2019$X <- T.K(wang_Algae2019$Temperature)

#The microzoo data in Wang et al. 2019
wang2019 <- read.csv('Wang2019Protist.csv')
wang2019$Temperature <- wang2019$Temperature - T0 
wang <- culldata(wang2019)
wang <- wang$newdat
