setwd('~/OneDrive/Kremer/')
library(plyr)
kb       <- 8.62e-5
Tref     <- 15
T0       <- 273.15
T.K      <- function(x) 1/((Tref+T0)*kb) - 1/((x + T0)*kb)

culldata <- function(newdat){
    newdat$X     <- T.K(newdat$Temperature)
    newdat$incl  <- NA  #Whether to include this in nonlinear fitting
    newdat$right <- NA  #Whether the data points fall to the right of Topt
    newdat$Tbar  <- NA  #Avg. temperature below optT
    newdat$bj    <- NA  #growth rate normalized to the reference temperature T0.
    newdat$EaL   <- NA  #Intraspecific E based on linear regression
    newdat$var1  <- NA  #sum of (xij - x__)^2
    newdat$XoptL <- NA  #Optimal temperature without nonlinear fitting
    newdat$umL   <- NA  #Maximal growth rate without nonlinear fitting
    
    PHYu <- ddply(newdat, .(ID), summarize,
                         ID       = ID[1],
                         Habitat  = Habitat[1],
                         Group    = Group[1],
                         Genus    = Genus[1],
                         Species  = Species[1],
                         lnVolume = mean(log(Volume), na.rm=T))
    
    N  <- nrow(PHYu)

    for (i in 1:N){
      code <- PHYu$ID[i]
      tmp  <- newdat[newdat$ID == code,]
      tmp  <- tmp[order(tmp$Temperature,tmp$Growth),]
      
      #Whether include the data in later analysis
      incl <- TRUE
      
      # At least 5 data points for nonlinear fitting
      if (length(tmp$Temperature) < 5) incl <- FALSE
      
      # At least 3 unique data points for nonlinear fitting
      if (length(unique(tmp$Temperature)) < 3) incl <- FALSE
      
     #  Average Optimal temperature:
      tmp$r.temp <- round(tmp$Temperature, 1)
      tmp.mean   <- tapply(tmp$Growth, tmp$r.temp, mean)
      
      #Distinct temperatures
      dtemp <- attributes(tmp.mean)$dimnames[[1]]
      dtemp <- as.numeric(dtemp)
      optT  <- dtemp[which.max(tmp.mean)] #Optimal temperature
      newdat[newdat$ID == code, ]$XoptL <- T.K(optT)
      newdat[newdat$ID == code, ]$umL   <- as.numeric(tmp.mean[which.max(tmp.mean)])
      
     # temperatures lower than optT
      Tl   <- tmp$Temperature[tmp$Temperature <= optT]
      if (length(unique(Tl)) < 3)  incl <- FALSE
      
      #Calculate Tbar (average temperature below optT)
      newdat[newdat$ID == code, ]$Tbar <- mean(T.K(Tl), na.rm=T)
      
      # temperatures higher than optT
      Th   <- tmp$Temperature[tmp$Temperature > optT]
      if (length(Th) < 1)  incl <- FALSE
      
      newdat[newdat$ID == code, ]$incl <- incl

      #Calculate bj
      if (incl){
        tmp1 <- tmp[tmp$Temperature <= optT & tmp$Growth > 0, ]
        LM1  <- lm(log(Growth) ~ T.K(Temperature), tmp1)
        newdat[newdat$ID == code, ]$bj  <- as.numeric(coef(LM1)[1])
        newdat[newdat$ID == code, ]$EaL <- as.numeric(coef(LM1)[2])
      }
      
      RIGHT <- rep(F, nrow(tmp))
      RIGHT[tmp$Temperature > optT]    <- T
      newdat[newdat$ID == code,]$right <- RIGHT
    }
    
    newdat2 <- newdat[!(newdat$right) & newdat$incl, ]
    GrandX  <- mean(newdat2$X, na.rm=T) #Calculate x__
    
    for (i in 1:N){
      code <- PHYu$ID[i]
      tmp  <- newdat[newdat$ID == code,]
      tmp  <- tmp[order(tmp$Temperature,tmp$Growth),]
      tmp  <- tmp[!tmp$right, ]
      if (tmp$incl[1]){
        newdat[newdat$ID == code, 'var1'] <- sum((tmp$X - GrandX)**2)
      }
    }
    return(newdat)
}

#Original zooplankton data
dat <- read.csv('Microzoo_23Mar2020.csv')
dat[dat$Habitat == 'freshwater','Habitat'] <- 'Freshwater'

zoodat <- culldata(dat)
zoodat <- zoodat[zoodat$incl & !is.na(zoodat$EaL), ]

zoodat2 <- zoodat[!(zoodat$right), ]

M.z <- nrow(zoodat2)  #Total number of valid observations



#Original phytoplankton data
#source('Phy_merge.R')
load('NEWphy24Mar2020.Rdata')
newdat$Sp    <- newdat$Species

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

newdat$Sp    <- as.character(get_spname(newdat$Species))
phydat       <- culldata(newdat)
phydat       <- phydat[phydat$incl & !is.na(phydat$EaL), ]

phydat2 <- phydat[!(phydat$right) & phydat$incl, ]
phydat2$mu <- phydat2$Growth
phydat2    <- phydat2[phydat2$mu > 0, ]

M.p <- nrow(phydat2)  #Total number of valid observations

#Var(X)
(VARx.p <- var(phydat2$X))


zoodat$Species <- as.character(zoodat$Species)
zoodat$Sp    <- as.character(get_spname(zoodat$Species))
zoodat$mu    <- zoodat$Growth
zoodat       <- zoodat[zoodat$mu > 0, ]
zoodat2      <- zoodat[!zoodat$right,]
#Var(X)
(VARx.z <- var(zoodat2$X))


