#Load prep.R
source('prep_data.R')
set.seed(10)

#Plot Ea examples and distributions of Ea
#Randomly select n (n = 8) taxa to plot Ea
YLIM <- c(-6, 2.2)
plot_Ea_example <- function(dat, XLAB = ''){
  dat <- dat %>% subset(Growth > 0)
  plot(dat$X, log(dat$Growth), pch = 16, cex=.5,
       xlim = c(-2.5,4),
       ylim = YLIM,
       type = 'n',
       xlab = XLAB,
       ylab = expression(paste(italic(ln)*' growth rate ') ) )
  
  N   <- 8
  IDs <- unique(dat$ID)
  IDs <- sample(IDs, N)
  for (i in 1:N){
    tmp <- dat %>% subset(ID == IDs[i])
    points(tmp$X, log(tmp$Growth), cex=.5, col=i, pch=i)
    z <- lm(log(Growth) ~ X, tmp[tmp$Growth > 0, ])
    newx <- data.frame(X=seq(min(tmp$X), max(tmp$X), .01))
    newy <- predict(z, newdata=newx)
    lines(newx$X, newy, col=i, lwd=.5, lty=i)
  }
}

#Plot histogram of Ea
#First, a function to generate summarized data
plot_Ea_hist <- function(oridat){
  oridat <- oridat %>% subset(Growth > 0)
  unidat <- ddply(oridat, .(ID), summarize,
                  ID       = ID[1],
                  Species  = Sp[1],
                  XoptL    = XoptL[1]     )
  
  for (i in 1:nrow(unidat)){
    code            <- unidat$ID[i]
    tmp             <- oridat[oridat$ID == code, ]
    LM1             <- lm(log(Growth) ~ X, tmp)
    unidat[i, 'bj'] <- as.numeric(coef(LM1)[1])
    unidat[i, 'Ea'] <- as.numeric(coef(LM1)[2])
  }
  
  hist(unidat$Ea, 
       xlim  =c(0, 1.5),
      # xlab  =expression(paste(E[a] * ' (eV)')),
       xlab  ='',
       ylab  ='',
       freq  =T,
       las   =1,
       yaxs  ='i',
       yaxt  ='n',
       breaks=seq(min(unidat$Ea)-.1, max(unidat$Ea)+.1, .1),
       cex.lab = 0.5,
       cex.axis= 0.8,
       main  = '')
  box()
}

#Function plotting out Einter (Ee)
Ncol  <- 10  #10 color codes
Topts <- seq(min(PEuk2$XoptL)-.01, max(PEuk2$XoptL)+.01, length.out=Ncol)
ToptC <- round(revTK(Topts),0) 

#Point colors
jet.colors   <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                    "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))  #A good color


#Symbol shapes
Shapes <- 16:(15+Ncol)

plot_Einter <- function(oridat, XLAB = '', plot.legend = FALSE){
  oridat <- oridat %>% subset(Growth > 0)
  unidat <- ddply(oridat, .(ID), summarize,
                  ID       = ID[1],
                  Species  = Sp[1],
                  XoptL    = XoptL[1]     )
  
 #Assign color code for optimal temperature (XoptL)
 unidat$Col_code <- NA
 
 for (i in 1:nrow(unidat)){
       code              <- unidat$ID[i]
       tmp               <- oridat[oridat$ID == code,]
       unidat[i, 'mj'  ] <- nrow(tmp)
       unidat[i, 'Tbar'] <- mean(tmp$X)
       LM1               <- lm(log(Growth) ~ X, tmp)
       unidat[i, 'bj']   <- as.numeric(coef(LM1)[1])
       
       wcol <- which(Topts > unidat$XoptL[i])
       unidat[i, 'Col_code'] <- wcol[1]-1
  }
    
  #Total number of observations of the original data
  M <- sum(unidat$mj)
  
  #Calculate probabilities for each taxon
  p <- unidat$mj/M

  #Calculate Ee 
  Ee <- lm(bj ~ Tbar,  unidat, weights = p)  #weighed linear regression
 
  #Plot
  plot(unidat$Tbar, unidat$bj, 
       pch  = Shapes[unidat$Col_code],
       col  = jet.colors(Ncol)[unidat$Col_code],
       cex  = 0.6,
       xaxs ='i',
       yaxs ='i',
       xlim = c(-2.5,4),
       ylim = YLIM,
       xlab = XLAB,
       ylab = expression(paste(italic(b) * ', dimensionless')))
  abline(Ee)
  
  #Add color legend
  if (plot.legend){
      txt <- paste(ToptC, 'ºC')
      legend('bottomleft', 
             legend = txt, 
             col = jet.colors(Ncol), 
             cex = .6,
             pch = Shapes)
      text(-2.4, -0.31, bquote(italic(T[opt])), cex=.6, pos=4)   
  }
}

plot_EL <- function(oridat, XLAB = '', plot.legend = FALSE){
  oridat <- oridat %>% subset(Growth > 0)
  unidat <- ddply(oridat, .(ID), summarize,
                  ID       = ID[1],
                  Species  = Sp[1],
                  XoptL    = XoptL[1]     )
  
 #Assign color code for optimal temperature (XoptL)
 unidat$Col_code <- NA
 
  for (i in 1:nrow(unidat)){
       code              <- unidat$ID[i]
       tmp               <- oridat[oridat$ID == code,]
       unidat[i, 'mj']   <- nrow(tmp)
       tmp$X1            <- tmp$X - tmp$XoptL[1]
       LM2               <- lm(log(Growth) ~ X1, tmp)
       unidat[i, 'ym']   <- as.numeric(coef(LM2)[1])
       
       wcol <- which(Topts > unidat$XoptL[i])
       unidat[i, 'Col_code'] <- wcol[1]-1
  }
    
  #Total number of observations of the original data
  M <- sum(unidat$mj)
  
  #Calculate probabilities for each taxon
  p <- unidat$mj/M

  #Calculate Ee 
  EL <- lm(ym ~ XoptL, unidat, weights = p)  #weighed linear regression
 
  #Plot
  plot(unidat$XoptL, unidat$ym, 
       pch  = Shapes[unidat$Col_code],
       col  = jet.colors(Ncol)[unidat$Col_code],
       cex  = 0.6,
       xaxs ='i',
       yaxs ='i',
       xlim = c(-2.5,4),
       ylim = YLIM,
       xlab = XLAB,
       ylab = expression(paste(italic(y[m]) * ', dimensionless')))
  abline(EL)
  
  #Add color legend
  if (plot.legend){
      txt <- paste(ToptC, 'ºC')
      legend('bottomleft', 
             legend = txt, 
             col = jet.colors(Ncol), 
             cex = .6,
             pch = Shapes)
      text(-2.2, -0.68, bquote(T[opt]), cex=.6, pos=4)   
  }
}

fname <- paste0('Fig2OLS_both.pdf')
pdf(fname, width = 7, height = 12, useDingbats=FALSE)
par(font.lab  = 1,
    family    = "serif",
    mgp       = c(2.2,1,0),
    oma       = c(3,2,2,2),
    mar       = c(3,4,1,.2),
    fig       = c(0,1,0,1)
   )
x1 <- c(0, 15)
y1 <- T.K(x1)

#All eukaryotic autotrophic protists
par(fig = c(0, .5, 0.75, 1))  
plot(PEuk2$X, log(PEuk2$Growth), pch = 16, cex=.5,
    xlim = c(-2.5,4),
    ylim = YLIM,
    xlab = expression(paste('1/'*italic(k[b])*'(1/'*italic(T[r])*' - 1/'*italic(T)*') or '*italic(x) *' ('*eV^-1*')')),
    ylab = expression(paste(italic(ln)*' growth rate (dimensionless)') ) )
points(PEuk[PEuk$right,]$X, 
       log(PEuk[PEuk$right,]$Growth), pch = 16, col = 'gray', cex=.5)
z <- lm(log(Growth) ~ X, PEuk2[PEuk2$Growth > 0, ])
abline(z)

text(-2.6, 1.85, pos=4, expression('A) Autotrophs '*italic(E[app])))
x1[length(x1)] <- paste0(x1[length(x1)],' ºC')
axis(3, at= y1, label = x1)  #Add comparisons of Celcius

par(fig = c(0.5, 1, 0.75, 1), new = T)  
plot(zoodat2$X, log(zoodat2$Growth), pch = 16,  cex=.5,
    xlim = c(-2.5,4),
    ylim = YLIM,
    xlab = expression(paste(italic(x) *' ('*eV^-1*')')),
    ylab = expression(paste(italic(ln)*' growth rate ') ) )
points(zoodat[zoodat$right,]$X, log(zoodat[zoodat$right,]$Growth), pch = 16, col = 'gray', cex=.5)
zoo.lm0 <- lm(log(Growth) ~ X, zoodat2[zoodat2$Growth > 0, ])
abline(zoo.lm0)

mtext(expression('B) Heterotrophs '*italic(E[app])), side=3, adj=0)

par(fig = c(0, .5, 0.5, .75), new = T)  
plot_Ea_example(PEuk2, XLAB = expression(paste(italic(x) *' ('*eV^-1*')' ) ) )

mtext(expression('C) Autotrophs '*italic(E[intra])), adj=0)
par(fig = c(0.22,.48, 0.53, .7), new = T)  #Plot inset
par(mgp = c(1,.1,0))
plot_Ea_hist(PEuk2)

par(mgp=c(2.2,1,0))
par(fig = c(0.5, 1, 0.5,.75), new = T)  
plot_Ea_example(zoodat2, XLAB = expression(paste(italic(x) *' ('*eV^-1*')') ) ) 
mtext(expression('D) Heterotrophs '*italic(E[intra])), adj=0)
par(fig = c(0.72,.98, 0.53, .7), new = T)  #Plot inset
par(mgp=c(1,.1,0))
plot_Ea_hist(zoodat2)

#Plot out Einter
par(mgp=c(2.2,1,0))
par(fig = c(0, 0.5, 0.25, 0.5), new = T)  
plot_Einter(PEuk2, XLAB = expression(paste(italic(bar(x)) * ' ('*eV^-1*')' )) )
mtext(expression('E) Autotrophs '*italic(E[inter])), adj=0)

par(fig = c(0.5, 1, 0.25, 0.5), new = T)  
plot_Einter(zoodat2, XLAB = expression(paste(italic(bar(x)) * ' ('*eV^-1*')')),
            plot.legend = T )
mtext(expression('F) Heterotrophs '*italic(E[inter])), adj=0)

#Plot out EL
par(fig = c(0, 0.5, 0, 0.25), new = T)  
plot_EL(PEuk2, XLAB = expression(paste(italic(x[m]) *' ('*eV^-1*')' ))  )

mtext(expression(bold(G)*') Autotrophs '*italic(E[L])), adj=0)

par(fig = c(0.5, 1, 0, 0.25), new = T)  
plot_EL(zoodat2, XLAB = expression(paste(italic(x[m]) *' ('*eV^-1*')' )) ) 

mtext(expression(bold(H)*') Heterotrophs '*italic(E[L])), adj=0)

dev.off()