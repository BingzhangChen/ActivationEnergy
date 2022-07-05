library(plyr) #Version 1.8.7
library(dplyr) #Version 1.0.9
kb       <- 8.62e-5
Tref     <- 15
T0       <- 273.15
T.K      <- function(x) 1/((Tref+T0)*kb) - 1/((x + T0)*kb)
revTK    <- function(x) 1/(1/(Tref + T0) - kb*x) - T0

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
    newdat       <- newdat[!is.na(newdat$Temperature),]
    newdat$X     <- T.K(newdat$Temperature)
    newdat$xbar  <- NA  #mean X (below optimal temperature)
    newdat$incl  <- NA  #Whether to include this in nonlinear fitting
    newdat$right <- NA  #Whether the data points fall to the right of Topt
    newdat$XoptL <- NA  #Optimal temperature without nonlinear fitting
    newdat$umL   <- NA  #Maximal growth rate without nonlinear fitting
    
    PHYu <- ddply(newdat, .(ID), summarize,
                         ID       = ID[1],
                         Group    = Group[1],
                         Species  = Species[1]
                     #    lnVolume = mean(log(Volume), na.rm=T)
                  )
    
    N  <- nrow(PHYu)

    for (i in 1:N){
      code <- PHYu$ID[i]
      tmp  <- newdat[newdat$ID == code,]
      tmp  <- tmp[!is.na(tmp$X), ]
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
      wx    <- which(newdat$ID == code & !is.na(newdat$X))
      newdat[wx, ]$XoptL <- T.K(optT)
      
      #Maximal growth rate at optimal temp
      newdat[wx, ]$umL   <- as.numeric(tmp.mean[which.max(tmp.mean)])
      
     # temperatures lower than optT
      Tl   <- tmp$Temperature[tmp$Temperature <= optT & tmp$Growth > 0]
    
      #At least 2 data points below or equal to optT    
      if (length(unique(Tl)) < 2)  incl <- FALSE
      
      # temperatures higher than optT
      Th   <- tmp$Temperature[tmp$Temperature > optT]
      #if (length(Th) < 1)  incl <- FALSE
      
      #Xbar for each taxon
      xbar <- mean(tmp[tmp$Temperature <= optT, 'X'])
      newdat[wx,]$xbar  <- rep(xbar, nrow(tmp))
      
      newdat[wx, ]$incl <- incl
      RIGHT                            <- rep(F, nrow(tmp))
      RIGHT[tmp$Temperature > optT]    <- T
      newdat[wx,]$right <- RIGHT
    }
    newdat$Sp    <- as.character(get_spname(newdat$Species))
    newdat       <- subset(newdat, incl ==TRUE)
    newdat2      <- subset(newdat, right==FALSE) #Remove supraoptimal temperature
    return(list(oridat = newdat, newdat = newdat2))
}

#Process phytoplankton data (autotrophic protists)
load('Merged_PHY.Rdata')
phydat       <- culldata(newdat)

#Final phyto data used for analysis
phydat2      <- phydat$newdat  #Contains zero growth rate
phydat       <- phydat$oridat  #Contains the data of supraoptimal temperatures

#Remove cyanobacteria (only eukaryotes)
PEuk         <-  phydat[ phydat$Group != 'Cyan',] #All Euk phyto. data
PEuk2        <- subset(phydat2, Group != 'Cyan')  #Eukaryotic phytoplankton after superoptimal data have been removed

#Cyanobacteria only
Cyn          <- phydat[ phydat$Group == 'Cyan',]  #All prokaryotic phyto. data
Cyn2         <- subset(phydat2, Group== 'Cyan')

#Processing heterotrophic protist data provided by David J. S. Montagnes
dat <- read.csv('HProtist.csv')
dat[dat$Habitat == 'freshwater','Habitat'] <- 'Freshwater'
zoodat       <- culldata(dat)
zoodat2      <- zoodat$newdat  #Final oridat of zoo used for analysis
zoodat       <- zoodat$oridat

#Insect data from Rezende and Bozinovic (2019)
Insect <- 'Insects.csv'
Insect <- read.csv(Insect)
Insect$Group <- NA
Insect$Growth<- Insect$Performance
Insect$ID    <- as.factor(Insect$ID)
Insect1      <- culldata(Insect)
Insect2      <- Insect1$newdat

#Bacteria data from Smith et al. NC 2019
Bac <- read.csv('Smith2019Bac.csv')

#Restrict data to only Specific Growth rate
Bac <- Bac %>% subset(StandardisedTraitName == 'Specific Growth Rate'
          )%>% subset(is.na(ConTrophic) | ConTrophic != 'producer' #Remove photosynthetical autrophic bacteria
          )%>% select(Longitude,
                      Latitude,
                      StandardisedTraitName, 
                      StandardisedTraitValue,
                      StandardisedTraitUnit,
                      AmbientTemp,
                      ConKingdom,
                      ConTemp,
                      ConGenus,
                      ConSpecies,
                      ConTrophic
          )%>% subset(!(ConGenus %in% c('Trichodesmium',  #Remove autotrophs
                                        'Synechococcus',
                                        'Rhodomicrobium',
                                        'Prochlorococcus',
                                        'Rhodospirillum',
                                        'Rhodopila',
                                        'Rhodobacter'
                                        )))

#Assign ID
Bac.ag <- Bac %>% select(Longitude,
                         Latitude,
                         StandardisedTraitName, 
                         StandardisedTraitUnit,
                         ConKingdom,
                         ConTemp,
                         ConGenus,
                         ConSpecies,
                         ConTrophic)

Bac.ag <- ddply(Bac.ag, .(Longitude, 
                          Latitude, 
                          StandardisedTraitName, 
                          StandardisedTraitUnit,
                          ConKingdom,
                          ConGenus,
                          ConSpecies,
                          ConTrophic), 
                summarize,
                Longitude             = Longitude[1],
                Latitude              = Latitude[1],
                ConGenus              = ConGenus[1],
                ConSpecies            = ConSpecies[1],
                StandardisedTraitName = StandardisedTraitName[1])

Bac.ag$ID <- row.names(Bac.ag)

#Assign ID back to Bac
Bac$ID <- NA
for (i in 1:nrow(Bac.ag)){
  w <- which(Bac$ConGenus  == Bac.ag[i, 'ConGenus'  ]&
             Bac$ConSpecies== Bac.ag[i, 'ConSpecies']&
           #  Bac$Longitude == Bac.ag[i, 'Longitude' ]&
           #  Bac$Latitude  == Bac.ag[i, 'Latitude'  ]&
             Bac$StandardisedTraitName == Bac.ag[i, 'StandardisedTraitName'])
  Bac[w, 'ID'] <- Bac.ag[i, 'ID']
}

Bac         <- Bac[!is.na(Bac$ID),]
Bac$ID      <- as.factor(Bac$ID)
Bac$Group   <- 'Bacteria'
Bac$Species <- NA
for (i in 1:nrow(Bac)){
  Bac$Species[i] <- paste(Bac$ConGenus[i], Bac$ConSpecies[i])  
}
Bac$Temperature <- Bac$ConTemp
Bac$Growth      <- Bac$StandardisedTraitValue * 86400

#Remove extremephiles
Bac             <- Bac[Bac$Temperature <= 50, ]

#Clean the dataset and calculate the several attributes needed
HBac            <- culldata(Bac)
HBac2           <- HBac$newdat
HBac            <- HBac$oridat    

#Again, remove the extremephiles
HBac2           <- HBac2[HBac2$XoptL <= T.K(35), ]
HBac            <-  HBac[HBac$XoptL  <= T.K(35), ]
  