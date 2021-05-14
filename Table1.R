OS <-  Sys.info()[['sysname']]
if (OS == 'Windows'){
  prefix <- "C:/Users/Dell/OneDrive - University of Strathclyde/"
}else{
  prefix <- "~/OneDrive/"
}

setwd(paste0(prefix,"Kremer"))
library(tidyverse)
library(foreach)
library(nlme)
source('prep_data.R')

#Remove cyanobacteria
phydat3 <- subset(phydat2, Group != 'Cyan')

#Compute Ea using nonlinear mixed-effect model
Ea_nlme <- function(dat){
  dat$ID <- factor(dat$ID)
  
  # #Linear mixed-effect model
   sdat <- subset(dat, Growth > 0)
  # 
   zL2 <- lme(log(Growth) ~ X,
             random = ~ 1 + X|ID,
             data   = sdat, method = 'REML')
  
   z2 <- nlme(Growth ~ A*exp(E*X), 
            data  = dat,
            fixed = A + E ~ 1,
            random= A + E ~ 1,
            groups= ~ID,
            start = c(A = 1, E = 0.65), method = 'REML')


  return(list(nlme = coefficients(summary(z2)), 
               lme = coefficients(summary(zL2)) ))
  #separate within-species and across-species
  # z3 <- nlme(Growth ~ A*exp(Ea*(X-xbar))*exp(Ee*xbar), 
  #            data  = dat,
  #            fixed = A + Ea + Ee ~ 1,
  #            random= A + Ea      ~ 1,
  #            groups= ~ID,
  #            start = c(A = 1, Ea = 0.65, Ee = -0.1), method = 'REML')
  
}

res_QQplot <- function(mod){
  
  #Examine residuals
  resid <- residuals(mod)
  
  #Plot Q-Q norm of standardized residuals of nlme model:
  stan_resid <- (resid-mean(resid))/sd(resid)
  qqnorm(stan_resid, pch=16, col=2, cex=.5)
  abline(0,1)
  st = shapiro.test(resid)
  return(list(w=st$statistic, p=st$p.value, SSE=sum(resid**2)))
}

#Weighed Covariance function
wcov <- function(x, y, w){
  stopifnot(length(x) == length(y) && length(x) == length(w))
  stopifnot(round(sum(w),0) == 1L)   
  xbar <- sum(x*w)
  ybar <- sum(y*w)
  
  return(sum(w*(x-xbar)*(y-ybar)))
}

decomp <- function(oridat, remove.nonpositive = T, epsilon=0.01){
  if (remove.nonpositive) {
    oridat <- oridat %>% subset(Growth > 0)
  }else{
    oridat[oridat$Growth <= 0, 'Growth'] <- epsilon
  }
  
  Eapp       <- numeric(5L)
  oridat$eps <- NA
  unidat     <- ddply(oridat, .(ID), summarize,
                    #  Habitat  = Habitat[1],
                      ID       = ID[1],
                      Species  = Sp[1],
                      XoptL    = XoptL[1]     )

  #grandmean of x
  GrandX <- mean(oridat$X)

  for (i in 1:nrow(unidat)){
    code              <- unidat$ID[i]
    tmp               <- oridat[oridat$ID == code,]
    unidat[i, 'mj'  ] <- nrow(tmp)
    unidat[i, 'Tbar'] <- mean(tmp$X)
    unidat[i, 'var1'] <- sum((tmp$X - GrandX)**2)
    
    #Calculate bj
    LM1  <- lm(log(Growth) ~ X, tmp)
    unidat[i, 'bj'] <- as.numeric(coef(LM1)[1])
    unidat[i, 'Ea'] <- as.numeric(coef(LM1)[2])
    
    #Add epsilon to oridat
    oridat[oridat$ID == code, 'eps'] <- LM1$residuals
  }
 
  #Median Ea
  Ea_median <- median(unidat$Ea)
  
  #Total number of observations of the original data
  M <- sum(unidat$mj)
  
  #Number of taxa
  n <- nrow(unidat)
  
  #Calculate probabilities for each taxon
  p <- unidat$mj/M

  #variance of x in the original data
  VARx  <- var(oridat$X)
  
  #variance of x_bar with unequal probabilities
  VARxbar <- sum(p*(unidat$Tbar - GrandX)^2)
  
  #First term
  Eapp[1] <- sum(unidat$Ea * unidat$var1)/M/VARx

  #Calculate Ee for phytoplankton
  Ee <- lm(bj ~ Tbar, unidat, weights = p)  #weighed linear regression
  
  #Obtain residuals (beta)
  beta <- unidat$bj - predict(Ee, newdata=unidat)

  #check whether mean of weighed residuals is 0
  #sum(p*beta)

  Ee  <- as.numeric(coef(Ee)[2])
  
  # weighed variance of x
  # VARxbar1 <- wcov(unidat$Tbar, unidat$Tbar, p) confirms the calculation above
  
  #Second term
  Eapp[2]  <- Ee * VARxbar/VARx
  
  #3rd term
  COVEXbar <-  wcov(unidat$Ea, unidat$Tbar, p) #weighed covariance between Ea and xbar 
  Eapp[3]  <-  GrandX  * COVEXbar  / VARx

  #4th term
  COVbetaxbar <- wcov(beta, unidat$Tbar, p)
  Eapp[4]     <- COVbetaxbar/VARx
  
  #5th term 
  COVepsX  <- cov(oridat$X, oridat$eps)
  Eapp[5]  <- COVepsX / VARx
  
  #run OLS regression for full data
  fullLM <- lm(log(Growth) ~ X, oridat) 
  
  return(list(n       = n,
              M       = M,
              VARx    = VARx,
              VARxbar = VARxbar,
              VARtheta= wcov(unidat$XoptL, unidat$XoptL, p),
              Ee      = Ee, 
              XMEAN   = GrandX,
              COVEXbar= COVEXbar,
              Ea_med  = Ea_median,
              EappLM  = as.numeric(coef(fullLM)[2]), 
              EappCal = Eapp))
} 

#phyto
boot <- function(dat, Nrep = 100){
  result <- foreach(i=1:Nrep, .combine='rbind') %do% {
    
    #Randomly sample 60% of the taxa
    IDs    <- dat %>% select(ID) %>% unique()
    x      <- sample(IDs[,1], 0.6*nrow(IDs))
    x      <- subset(dat, ID %in% x)   #subsample
    x      <- decomp(x)
    
    c(x$EappLM, x$Ea_med, x$Ee, sum(x$EappCal), x$EappCal[1], x$EappCal[2])
  }
  MEAN <- apply(result,2,mean)
  SE   <- apply(result,2,sd)
  
  return(list(mean=MEAN, se=SE))
}

#OLS regression by considering cell size (Table S1)
OLSsize <- function(dat){
  sdat <- dat %>% subset(Growth > 0)
  Z <- lm(log(Growth) ~ X + log(Volume), sdat)
  
  Z <- summary(Z)
  return(list(N=sum(Z$df[1:2]), R2=Z$r.squared, 
         Eapp   = coefficients(Z)[2,1],
         SE_Eapp= coefficients(Z)[2,2],
         alpha  = coefficients(Z)[3,1],
         SE_alpha = coefficients(Z)[3,2]))
}

aut0   <- decomp(phydat2, remove.nonpositive = F, epsilon = 0.001)
aut0SE <-   boot(phydat2) 

aut    <- decomp(phydat3)
autSE  <-   boot(phydat3)
autnlme<- Ea_nlme(phydat3)
OLSsize(phydat3)

autw   <- decomp(phydat4)
autwSE <-   boot(phydat4)

#MicroZooplankton
het    <- decomp(zoodat2)
hetSE  <-   boot(zoodat2)
hetnlme<- Ea_nlme(zoodat2)
OLSsize(zoodat2)

#micrzoo data in Chen & Laws 2017
het2017 <- decomp(dat2017)

#algae data in Wang et al. 2019
WangAlgae2019 <- decomp(wangAlgae)


#micrzoo data in Wang et al. 2019
Wang2019 <- decomp(wang)

#Plot residual plots
QQplot <- function(dat, method='nlme', label='a'){
  #NLME residual plots
  dat$ID <- factor(dat$ID)
  
  if (method == 'nlme'){
    mod <- nlme(
      Growth ~ A * exp(E * X),
      data  = dat,
      fixed = A + E ~ 1,
      random = A + E ~ 1,
      groups = ~ ID,
      start = c(A = 1, E = 0.65),
      method = 'REML'
    )
  }else if (method == 'lme'){
    #LME residual plots
    sdat   <- dat %>% subset(Growth > 0)
    mod    <- lme(log(Growth) ~ X, 
                   random = ~ 1 + X|ID,
                   data   = sdat, method = 'REML')
  }else if (method == 'lfe'){
    #OLS residual plots
    sdat <- dat %>% subset(Growth > 0)
    mod  <- lm(log(Growth) ~ ID * X, data = sdat)
  }else{
    stop('Method wrong!')
  }
  
  #Examine residuals
  resid <- residuals(mod)
  
  #Plot Q-Q norm of standardized residuals of nlme model:
  stan_resid <- (resid-mean(resid))/sd(resid)
  qqnorm(stan_resid, xlim=c(-4,4), ylim=c(-4,4),pch=16, col=2, cex=.5, main='')
  abline(0,1)
  mtext(label, adj=0, cex=.8)
  
}

#Plot both auto- and heterotrophs
pdf('Residual_QQplot.pdf', width = 6, height = 4)
par(font.lab  = 1,
    family    = "serif",
    mgp       = c(2.2,1,0),
    mfrow     = c(2,3),
    oma       = c(4,2,2,2),
    mar       = c(4,4,1,.2))

QQplot(phydat3, 'nlme', 'a) NLME, autotroph')
QQplot(phydat3, 'lme',  'b) LME,  autotroph')
QQplot(phydat3, 'lfe',  'c) LFE,  autotroph')

QQplot(zoodat2, 'nlme', 'd) NLME, heterotroph')
QQplot(zoodat2, 'lme',  'e) LME,  heterotroph')
QQplot(zoodat2, 'lfe',  'f) LFE,  heterotroph')

dev.off()

#Plot each taxon
TODAY      <- Sys.Date()
phypdffile <- paste0('phyto_bytaxon_linear',    TODAY, '.pdf')
zoopdffile <- paste0('microzoo_bytaxon_linear', TODAY, '.pdf')

plot_taxon <- function(dat, filename, caption = ''){

  dat$ID <- factor(dat$ID)
  
  #NLME
  NLME <- nlme(
      Growth ~ A * exp(E * X),
      data   = dat,
      fixed  = A + E ~ 1,
      random = A + E ~ 1,
      groups = ~ ID,
      start  = c(A = 1, E = 0.65),
      method = 'REML')
   
  result.nlme <- coefficients(NLME)
  sdat   <- dat %>% subset(Growth > 0)
  
  #Linear mixed-effect model
  LME    <- lme(log(Growth) ~ X, 
                  random = ~ 1 + X|ID,
                  data   = sdat, method = 'REML')
  result.lme <- coefficients(LME)
  Ntaxon <- length(levels(dat$ID))
    
  YLAB <- expression(paste("Growth rate (" * d ^ -1 * ")"))
  
  pdf(filename,
      width = 7,
      height = 9,
      paper = 'a4')
  par(font.lab  = 1,
      family    = "serif",
      mgp       = c(2.2, 1, 0),
      mfrow     = c(4, 3),
      mar       = c(4, 4, 3.5, .2),
      oma       = c(4,4,0,0))
  

  for (i in 1:Ntaxon) {
 
      id   <- levels(dat$ID)[i]
      tmp  <- dat %>% subset(ID == id)
      
      #Linear fixed-effect model 
      LFE  <- lm(log(Growth) ~ X, data = tmp[tmp$Growth > 0,])
      
      #Essential to reset the newx
      minx <- min(tmp$Temperature)
      maxx <- max(tmp$Temperature)
      newx <- data.frame(Temperature = seq(minx, maxx, 0.01))
      newy <- data.frame('NLME'=numeric(nrow(newx)),
                         'LME' =numeric(nrow(newx)),
                         'LFE' =numeric(nrow(newx)))
      #Find index for this taxon
      w <- which(rownames(result.nlme)==id)
      newy$NLME <- result.nlme[w, 'A'] * exp(result.nlme[w, 'E'] * T.K(newx$Temperature))
      
      w <- which(rownames(result.lme)==id)
      newy$LME  <- exp(result.lme[w, 1] + result.lme[w, 2] * T.K(newx$Temperature))
      
      newy$LFE  <- exp(coef(LFE)[1] + coef(LFE)[2] * T.K(newx$Temperature))
      
      #Plotting original data points
      umax <- max(tmp$Growth, na.rm = T)

      plot(tmp$Temperature, tmp$Growth,
           las = 1,
           xlim = c(-2,   35),
           ylim = c(-0.05, umax + .2),
           xlab = 'Temperature (ºC)',
           ylab = YLAB,
           cex.lab = 1.1)
      lines(newx$Temperature, newy$NLME, col=2, lwd=1.3)
      lines(newx$Temperature, newy$LME,  col=3, lwd=1.3)
      lines(newx$Temperature, newy$LFE,  col='blue', lwd=1.3)
      txt  <- paste0(tmp$Sp[1])
      mtext(txt, adj = 0, cex=.7)  #Add caption
      text(34.4, 0, paste('Ea.nlme = ', round(result.nlme[w, 'E'],2)), pos=2, col=2)
      text(34.4, umax*0.15,paste('Ea.lme = ', round( result.lme[w, 2],2)), pos=2, col=3)
      text(34.4, umax*0.3, paste('Ea.lfe = ', round( coef(LFE)[2], 2)), pos=2, col='blue')
  }
  #Add figure caption
  mtext(caption, side = 1, line = 0, outer = T, adj = 0)
  dev.off()
}

plot_taxon(phydat3, phypdffile, 'Fig. S1. Growth rate versus temperature plots for each autotrophic taxon')
plot_taxon(zoodat2, zoopdffile, 'Fig. S2. Growth rate versus temperature plots for each heterotrophic taxon')
