setwd('~/OneDrive/Kremer')
source('../Global_PP/Rscripts/Johnson.R')
T.K <- function(x, Tref = 15,kb=8.62E-5) 1/((Tref+273)*kb) - 1/((x + 273)*kb)

temp    <- seq(-2, 35, 0.1)
temp    <- T.K(temp)
numtemp <- length(temp)
x0    <- seq(10,30,10)
x1    <- x0
x1[length(x1)] <- paste0(x1[length(x1)],' ºC')

#Simulate an artificial dataset
set.seed(12345)
randata <- function(N1=1000){
  X1 <- matrix(NA, nr=N1, nc=2)
  colnames(X1) <- c('X', 'mu')
  X1 <- as.data.frame(X1)
  
  X1$X <- rnorm(N1, 0, .8)
  err1 <- rnorm(N1, 0, .8)
  X1$mu<- X1$X * 0.65 + err1
  return(X1)
}

Z1   <- randata()
ols1 <- lm(mu ~ X, Z1)

#quantile regression
library(quantreg)
qr1 = rq(mu ~ X, Z1, tau = 0.99)

pdf('Fig1_E_definition.pdf', height=4, width=9)
op <- par(font.lab = 1,
           family  = "serif",
           cex.lab = 1.2,
           oma     = c(4,4,4,0),
           cex.axis= 1.2,
           mgp     = c(2.2,1,0),
           mar     = c(2,3,3,1),
           mfrow   = c(1,3))

plot(Z1$X, Z1$mu, 
    xlim = c(-2.5, 3),
    ylim = c(-3, 3),
    xlab = '', ylab = '', pch = 16, cex=.8)
abline(ols1, lwd=2, col=2)
#abline(qr1,  lwd=2, lty=3)
txt <- expression(paste('a) '* italic(E[app])))
mtext(txt, adj=0, cex=1.2)
#legend('topleft', legend = c('OLS', 'QR'), lty=c(1,3), lwd=2, col=c(2,1))
axis(3, at= T.K(x0), label = x1)

#Ea: intraspecific activation energy
Z2 = data.frame(X=seq(-2, 2, length.out=10), mu=0)
Z2$mu <- sapply(1:nrow(Z2),function(i)Johnson2(Z2$X[i], Topt=T.K(25),
                mumax=1, Ea =0.8,Eh=5))
err3  <- rnorm(nrow(Z2), 0, .01)
Z2$mu <- Z2$mu + err3

plot(Z2$X, log(Z2$mu), xlab='', ylab='',
    xlim = c(-2.5, 3),
    ylim = c(-3, 3))
Z2b = Z2[Z2$X <= T.K(20),]
abline(lm(log(mu) ~ X, Z2b), lwd=2, col=2)
axis(3, at= T.K(x0), label = x1)

txt <- expression(paste('b) '* italic(E[a])))
mtext(txt, adj=0, cex=1.2)
#
#Ei: interspecific activation energy
mumax_ =function(mu0=1,Topt, Ei = 0, kb = 8.62e-5, Tr = 15, T0=273.15) mu0 * exp(Ei/kb *(1/(Tr+T0) - 1/(T0+Topt)))
#For phytoplankton, lower Ei
Topto <- c(5, 10, 20, 25)
mumaxP = mumax_(Topt=Topto, Ei = 0.2)

plot(temp,temp**2,
     xlim = c(-2.5, 3), 
     ylim = c(-3, 2), 
     xlab = '',
     ylab = '',
     type = 'n')

axis(3, at= T.K(x0), label = x1)
txt <- expression(paste('c) '* italic(E[i])))
mtext(txt, adj=0, cex=1.2)

for (j in 1:length(Topto)){
  mu_ <- sapply(1:numtemp,function(i)Johnson2(temp[i], Topt=T.K(Topto)[j],
                mumax=mumaxP[j], Ea =0.8,Eh=5))
  points(temp, log(mu_), type = 'l', col=j)
}
points(T.K(Topto), log(mumaxP))
Eapp1 <- lm(log(mumaxP) ~ T.K(Topto))
abline(Eapp1, lwd=2, col=2)

txt1 <- bquote('Normalized Boltzmann temperature ('~eV^-1~')')
txt2 <- bquote('Ln growth rate ('~d^-1~')')
txt3 <- 'Temperature (ºC)'
mtext(side=1, txt1, outer=T,line=1)
mtext(side=2, txt2, outer=T,line=.1)
mtext(side=3, txt3, outer=T,line=.1)
#Add histogram of Ea
Eas <- rnorm(100, 0.8, .3)
par(fig = c(.38,.55,.45,.9), new=T)
hist(Eas, xlab='Activation energy (eV)', ylab='', xlim=c(0,1.5), 
    main=expression(paste('Distribution of '* italic(E[a]))))

dev.off()

