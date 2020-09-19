############################################
setwd('~/OneDrive/Kremer')
#Apply NLS for each strain:
Johnson <-  function(temp, Topt, Lnmumax, LnEa, LnEh){
      T0 <- 273.15
      kb <- 8.62E-5
      Ea <- exp(LnEa)
      Eh <- exp(LnEh)
      Ed <- Eh - Ea
      b  <- T.K(temp) - T.K(Topt)   #Standardized temperature (ivT)
      h  <- (Ea/Ed + 1) * exp(Ea * b)/(1+ Ea/Ed * exp(Eh * b) )
      h  <- exp(Lnmumax)* h 
      return(h)
}

source('prep_new.R')
library(plyr)
nls.control(maxiter = 2000, tol = 1e-03, 
    minFactor = 1/204800,
    printEval = FALSE, warnOnly = TRUE)
TODAY      <- Sys.Date()
phypdffile <- paste0('phyto_bystrain', TODAY, '.pdf')
zoopdffile <- paste0('microzoo_bystrain', TODAY, '.pdf')

get_NLS <- function(newdat, filename){
  newdat[newdat$Growth < 0, 'Growth']=0

  PHYu <- ddply(newdat, .(ID), summarize,
                     ID       = ID[1],
                     Habitat  = Habitat[1],
                     Group    = Group[1],
		                 Genus    = Genus[1],
                     Species  = Sp[1],
                     lnVolume = mean(log(Volume), na.rm=T),
                     Tbar     = Tbar[1],
                     bj       = bj[1],
                     EaL      = EaL[1],
                     var1     = var1[1],
                     XoptL    = XoptL[1],
                     umL      = umL[1]
                )

  #Initialize Estimate:
  PHYu$Lnum   <- NA
  PHYu$LnEa   <- NA
  PHYu$LnEh   <- NA
  PHYu$Topt   <- NA
  
  N    <- nrow(PHYu)
  umin <- 0.01
  YLAB <- expression(paste("Growth rate (" * d ^ -1 * ")"))
  
  pdf(filename,
      width = 7,
      height = 9,
      paper = 'a4')
  par(
    font.lab  = 1,
    family    = "serif",
    mgp       = c(2.2, 1, 0),
    mfrow     = c(4, 3),
    mar       = c(4, 4, 3.5, .2)
  )
  
  for (i in 1:N) {
    code <- PHYu$ID[i]
    tmp  <- newdat[newdat$ID == code, ]
    tmp  <- tmp[order(tmp$Temperature, tmp$Growth), ]
    
    #Johnson nls regression:
    if (tmp$incl[1]) {
      # Run nls:
      mumin  <- 1e-2
      muMAX  <- 10
      Eamin  <- 1e-2
      Eamax  <- 4
      Ehmin  <- .1
      Ehmax  <- 20
      Tomin  <- -6
      Tomax  <- 50
      algs   <- c('port', 'default', 'plinear')
      z   <- NULL
      for (k in 1:length(algs)) {
        try(z <- nls(
          Growth ~ Johnson(Temperature, Topt, Lnum, LnEa, LnEh),
          data = tmp,
          start = list(
            Lnum = max(log(tmp$Growth), na.rm = T),
            LnEa = log(0.6),
            LnEh = log(3),
            Topt = optT
          ),
          lower = list(
            Lnum = log(mumin),
            LnEa = log(Eamin),
            LnEh = log(Ehmin),
            Topt = Tomin
          ),
          upper = list(
            Lnum = log(muMAX),
            LnEa = log(Eamax),
            LnEh = log(Ehmax),
            Topt = Tomax
          ),
          algorithm = algs[k]
        ))
        if (!is.null(z))
          break
      }
      
      #Essential to reset the newx
      minx <- min(tmp$Temperature)
      maxx <- max(tmp$Temperature)
      newx <- data.frame(Temperature = seq(minx, maxx, 0.01))
      
      optT  <- tmp$Temperature[which.max(tmp$Growth)] #Optimal temperature
      newx1 <- data.frame(Temperature = seq(minx, optT, 0.01))
      newx2 <- T.K(newx1$Temperature)
      newyL <- exp(tmp$bj[1] + tmp$EaL[1] * newx2)
        
      if (!is.null(z)) {
        newy <-
          predict(z, newdata = data.frame(Temperature = newx$Temperature))
        
        #Plotting original data points
        umax <- max(c(tmp$Growth, newy), na.rm = T)
        
        for (p in c('Topt', 'Lnum', 'LnEa', 'LnEh')) {
          PHYu[PHYu$ID == code, p] <- as.numeric(coef(z)[p])
        }
      }else{
        umax <- max(tmp$Growth, na.rm = T)
      }
      plot(tmp$Temperature, tmp$Growth,
          las = 1,
          xlim = c(-2,   35),
          ylim = c(umin, umax),
          xlab = 'Temperature (ºC)',
          ylab = YLAB,
          cex.lab = 1.4)
      lines(newx1$Temperature, newyL, col=3)
      
      if (!is.null(z)) {
        points(newx$Temperature, newy, type = 'l', col = 2)
        EaNL <- exp(as.numeric(coef(z)['LnEa']))
        text(-2, umax*0.6, 
             paste('EaNL = ', round(EaNL,2)), pos=4)
      }
      txt  <- paste0(tmp$Sp[1])
      mtext(txt, adj = 0)  #Add caption
      text(-2, umax*0.9, 
           paste('EaL = ', round(tmp$EaL[1],2)), pos=4)
      
    }
  }
  
  dev.off()
  return(PHYu)
}

#nonlinear regression results:
PHYnls = get_NLS(phydat, phypdffile)
PHYnls$Habitat = as.factor(as.character(PHYnls$Habitat))
PHYnls$Group   = as.factor(as.character(PHYnls$Group))
save(PHYnls, file = paste0('PHYnls', TODAY, '.Rdata'))
PHYnls$EaNL <- exp(PHYnls$LnEa)

#median (EaL)
(EaLmed.p <- median(PHYnls$EaL))

#Calculate SE of median using bootstrapping
SE_Ea <- function(DAT, vname='EaL', Nrep = 500, Nsample=50){
  
  Ea <- numeric(Nrep)
  for (i in 1:Nrep){
    db <- sample(1:nrow(DAT), Nsample, replace=T)
    dd <- DAT[db, vname]
    Ea[i] <- median(dd, na.rm=T)
  }
  return(list(EaMed = mean(Ea), EaSE = sd(Ea)))
}

(EaL.p <- SE_Ea(PHYnls))

#median Ea (from nonlinear regression)
(ENLmed.p <- median(exp(PHYnls$LnEa), na.rm=T))
(EaNL.p    <- SE_Ea(PHYnls, 'EaNL'))

#Var(X_bar)
(VARxbar.p <- var(PHYnls$Tbar))

#Grand mean of X
(GrandX.p <- mean(phydat2$X))

#Cov(E, x_bar)
(COVEXbar.p <- cov(PHYnls$EaL, PHYnls$Tbar))

#First term:
F1.p <- sum(PHYnls$EaL * PHYnls$var1)
(F1.p <- F1.p/M.p/VARx.p)

#Second term:
(F2.p <-  GrandX.p * COVEXbar.p / VARx.p)

#Calculate K for phytoplankton
Kp.p  <- lm(bj ~ Tbar, PHYnls)
(Kp.p  <- coef(Kp.p)[2])
Kp1.p <- lm(bj ~ T.K(Topt), PHYnls)

#Calculate Ki using linear calculations
Ki1.p <- lm(log(umL) ~ XoptL, PHYnls)

#Calculate Ki using nonlinear calculations
(Ki2.p <- lm(Lnum ~ T.K(Topt), PHYnls))

#Third term:
(F3.p <- Kp.p * VARxbar.p/VARx.p)

plot(T.K(PHYnls$Topt), PHYnls$Tbar)

# #Get mean and covariance matrices
# summary(PHYnls)
# vnames <- c('Topt', 'Lnum', 'LnEa', 'LnEh')
# P <- PHYnls[, vnames]
# round(var(P), 2)
# CI95 = function(x) {
#   return(c(mean(x) - 1.96 * sd(x), mean(x) + 1.96 * sd(x)))
# }
# exp(quantile(PHYnls$LnEa, probs = c(0.025, 0.5, 0.975)))

#MicroZooplankton
ZOOnls <- get_NLS(zoodat, zoopdffile)
ZOOnls$EaNL <- exp(ZOOnls$LnEa)
save(ZOOnls, file = paste0('ZOOnls', TODAY, '.Rdata'))

#median (EaL)
(ELmed.z <- median(ZOOnls$EaL))

(EaL.z <- SE_Ea(ZOOnls))

#Compare EaL between phyto and zoo using Mann-Whitney test
wilcox.test(ZOOnls$EaL, PHYnls$EaL) #Zoo EaL significantly higher

#median Ea (from nonlinear regression)
(Eamed.z <- median(exp(ZOOnls$LnEa), na.rm=T))
(EaNL.z  <- SE_Ea(ZOOnls, 'EaNL'))
wilcox.test(ZOOnls$EaNL, PHYnls$EaNL)  #Zooplankton EaNL marginally higher

#Grand mean of X
(GrandX.z <- mean(zoodat2$X))

#Var(X_bar)
(VARxbar.z <- var(ZOOnls$Tbar))

#Cov(E, x_bar)
(COVEXbar.z <- cov(ZOOnls$EaL, ZOOnls$Tbar))

#First term:
F1.z <- sum(ZOOnls$EaL * ZOOnls$var1)
(F1.z <- F1.z/M.z/VARx.z)

#Second term:
(F2.z <-  GrandX.z * COVEXbar.z / VARx.z)

#Calculate K for zooplankton
Kp.z  <- lm(bj ~ Tbar,      ZOOnls)
(Kp.z <- coef(Kp.z)[2])

(Kp1.z <- lm(bj ~ T.K(Topt), ZOOnls))

#Calculate Ki using nonlinear calculations
(Ki2.z <- lm(Lnum ~ T.K(Topt), ZOOnls))

#Third term:
(F3.z <- Kp.z * VARxbar.z/VARx.z)

Z <- ZOOnls[,vnames]
round(var(Z),2)
exp(CI95(ZOOnls$LnEa))


#Perform Wilcoxon 2 sample tests to compare medians
wilcox.test(exp(ZOOnls$LnEa), exp(PHYnls$LnEa))

#Plot a boxplot of Ea of both Phyto and microzoo
fname   <- paste0('Fig4Boxplot_Ea.pdf')
pdf(fname, width = 9, height = 5)
par(font.lab  = 1,
    family    = "serif",
    mgp       = c(2.2,1,0),
    mfrow     = c(1,2),
    oma       = c(4,2,2,2),
    mar       = c(4,4,1,.2))
Y <- boxplot(PHYnls$EaL, ZOOnls$EaL, 
        outline=FALSE,
        names = c('Phyto', 'Microzoo'),
        ylab = expression(paste(italic(E)[aL]*' (eV)')))
mtext('a)', adj = 0)
X <- boxplot(PHYnls$EaNL, ZOOnls$EaNL, 
          outline=FALSE,
          names = c('Phyto', 'Microzoo'),
          ylab = expression(paste(italic(E)[aNL]*' (eV)')))
mtext('b)', adj = 0)
dev.off()
