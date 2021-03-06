setwd('~/Onedrive/Kremer')
#Load prep.R
source('prep_new.R')
#Try phytoplankton data as a whole
phy.lm0      <- lm(log(mu) ~ X, phydat2)
EpOLS        <- coef(phy.lm0)[2] #Simplest Eapp
sdEpOLS      <- coef(summary(phy.lm0))[2,2]

phy.lm1      <- lm(log(mu) ~ X + log(Volume), phydat2)
#Check the effect of size
library(lme4)
phy.lm2 <- lmer(log(mu) ~ X + log(Volume) + (1|Group), phydat2)
phy.lm3 <- lmer(log(mu) ~ X + log(Volume) + (1|Group/Species), phydat2) #Best
phy.lm4 <- lmer(log(mu) ~ X + log(Volume) + (1|Species), phydat2)
AIC(phy.lm0, phy.lm1, phy.lm2, phy.lm3, phy.lm4)



#Try zooplankton data as a whole
zoo.avgVol   <- mean(log(zoodat2$Volume), na.rm=T)
zoo.lm0      <- lm(log(mu) ~ X, zoodat2)
EzOLS        <- coef(zoo.lm0)[2]
sdEzOLS      <- coef(summary(zoo.lm0))[2,2]

#Check the effect of size
zoo.lm1 <-   lm(log(mu) ~ X + log(Volume),                     zoodat2)
zoo.lm2 <- lmer(log(mu) ~ X + log(Volume) + (1|Group),         zoodat2)
zoo.lm3 <- lmer(log(mu) ~ X + log(Volume) + (1|Group/Species), zoodat2)
zoo.lm4 <- lmer(log(mu) ~ X + log(Volume) + (1|Species),       zoodat2) #best


AIC(zoo.lm0, zoo.lm1, zoo.lm2, zoo.lm3, zoo.lm4)

fname   <- paste0('Fig3OLS_both.pdf')
pdf(fname, width = 5, height = 7)
par(font.lab  = 1,
    family    = "serif",
    mgp       = c(2.2,1,0),
    mfrow     = c(2,1),
    oma       = c(4,2,2,2),
    mar       = c(2,4,1,.2))
x1   <- seq(0,30,10)
y1   <- T.K(x1)
best.results <- summary(phy.lm3)$coefficients
LMEslope <- best.results[2,1]
LMEint   <- best.results[1,1]
phy.avgVol <- mean(log(phydat2$Volume), na.rm=T)
LMEint2  <- LMEint + best.results[3,1] * phy.avgVol
plot(phydat$X, log(phydat$Growth), pch = 16, cex=.5,
    xlim = c(-2.5,4),
    ylim = c(-6, 2),
    xlab = '',
    ylab = expression(paste('Ln growth rate ('*d^-1*')')))
points(phydat[phydat$right,]$X, log(phydat[phydat$right,]$Growth), pch = 16, col = 'gray', cex=.5)
abline(phy.lm0)

#Add regression line of LME model
abline(a = LMEint2, b = LMEslope, col = 2)

text(-2.6, 1.85, pos=4, 'a) Phytoplankton')
x1[length(x1)]=paste0(x1[length(x1)],' ºC')
axis(3, at= y1, label = x1)
legend('bottomright', c('OLS', 'LME'), lty = 1, col = 1:2)

par(mar = c(4,4,1,.2))
plot(zoodat$X, log(zoodat$mu), pch = 16,  cex=.5,
    xlim = c(-2.5,4),
    ylim = c(-6, 2),
    xlab = expression(paste('1/'*k[b]*'(1/'*T[r]*' - 1/T) or '*italic(x))),
    ylab = expression(paste('Ln growth rate ('*d^-1*')')))
points(zoodat[zoodat$right,]$X, log(zoodat[zoodat$right,]$mu), pch = 16, col = 'gray', cex=.5)
abline(zoo.lm0)

best.results <- summary(zoo.lm4)$coefficients
LMEslope <- best.results[2,1]
LMEint   <- best.results[1,1]
LMEint2  <- LMEint + best.results[3,1] * zoo.avgVol
abline(a = LMEint2, b = LMEslope, col = 2)
#zoodat3 = zoodat2[zoodat2$Temperature <= 25,]
#qr2 = rq(log(mu) ~ X, zoodat3, tau = 0.95)
text(-2.6, 1.8, pos=4, 'b) Microzooplankton')
dev.off()

