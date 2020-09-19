setwd('~/OneDrive/Kremer') 
T.K <- function(x, Tref = 15,kb=8.62E-5) 1/((Tref+273)*kb) - 1/((x + 273)*kb)

#Plot four typical linear growth~temperature relationships and the fitness landscape under a given environmental temperature
#different xbar & Topt, 
dx   <- 7.5
xbar <- seq(dx, 30-dx, length.out = 4)
theta<- xbar + dx/2
xmin <- xbar - dx/2
xbar <- T.K(xbar)
theta<- T.K(theta)
xmin <- T.K(xmin)
Ea   <- 0.65
b0   <- log(2)
b    <- b0 - Ea * theta

pdf('Fig2Eapp.pdf', height=6, width=9)
op <- par(font.lab = 1,
           family  = "serif",
           cex.lab = 1.2,
           oma     = c(2,2,0.1,0.1),
           cex.axis= 1.2,
           mgp     = c(1,.1,0),
           mar     = c(3,3,2,0.1),
           mfrow   = c(1,2))

temp    <- seq(-2, 35, 0.1)
temp    <- T.K(temp)
yrange  <- c(0.5, 3)

plot(temp,temp**2, xlim=c(T.K(0), T.K(30)), ylim=log(yrange), 
     xlab = expression(paste('x ('*eV^-1*')')),
     ylab = expression(paste('y (ln '*d^-1*')')),
     xaxs = "i",
     yaxs = "i",
     xaxt = 'n',
     yaxt = 'n',
     type = 'n')

pool <- data.frame(X = c(), Y = c())
#Plot each TPC
for (i in 1:length(xbar)){
  curve(Ea*x + b[i], from = xmin[i], to = theta[i], add = T)
  #add data points
  X <- seq(xmin[i], theta[i], length.out = 4)
  Y <- Ea*X + b[i]
  points(X,Y)
  pool <- rbind(pool, data.frame(X,Y))
}

LMp <- lm(Y ~ X, pool)
abline(LMp, col=2, lwd=1.5)
abline(v=0, lty = 3, col = 2)
#add xbar [j]
segments(x0 = xbar[1], y0 = log(0.5), 
         x1 = xbar[1], y1 = Ea*xbar[1] + b[1],
         lty = 2, col = 'blue')

#add b[j]
segments(x0 = theta[1], y0 = Ea*theta[1] + b[1], 
         x1 = 0,   y1 = b[1],
         lty = 2, col = 2)

segments(x0 = 0,      y0 = b[1], 
         x1 = T.K(0), y1 = b[1],
         lty = 2, col = 2)

#add theta[j] and ym[j]
segments(x0 = theta[1], y0 = Ea*theta[1] + b[1], 
         x1 = theta[1], y1 = log(yrange[1]),
         lty = 2, col = 3)
segments(x0 = theta[1], y0 = Ea*theta[1] + b[1], 
         x1 = T.K(0),   y1 = Ea*theta[1] + b[1],
         lty = 2, col = 3)

mtext(expression(paste('a) '*italic(E[app])*', '*italic(E[aL])*', & '*italic(E[i]))), adj=0, cex=1.4)
mtext(expression(paste(bar(x[j]))), side = 1, adj = .25, col = 'blue', line = .1)
mtext(expression(paste(theta[j])),  side = 1, adj = .4,  col = 3,      line = .1)
mtext('0', side = 1, adj = .52, col = 2)
mtext(expression(paste(b[j])),      side = 2, adj = .97, col = 2)
mtext(expression(paste(y[mj])),     side = 2, adj = .77,  col = 3)

#b) plot K (b ~ xbar)
plot(xbar, b, 
     xlab = expression(paste(bar(x)*' ('*eV^-1*')')),
     ylab = expression(paste('b (ln '*d^-1*')')),
     xaxt = 'n',
     yaxt = 'n')
abline(lm(b ~ xbar), lwd = 1.5, col = 2 )
mtext(expression(paste('b) '* italic(K))), adj=0, cex=1.4)

dev.off()
