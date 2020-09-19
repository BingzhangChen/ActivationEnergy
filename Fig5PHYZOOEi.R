setwd('~/OneDrive/Kremer/')
library(lme4)
load('PHYnls2020-09-17.Rdata')
load('ZOOnls2020-09-17.Rdata')
kb       <- 8.62e-5
Tref     <- 15
T0       <- 273.15
T.K      <- function(x, Tref = 15,kb=8.62E-5) 1/((Tref+T0)*kb) - 1/((x + T0)*kb)

pdf('Fig5PHYZOOEi.pdf', height = 7, width = 7)
op <- par(font.lab = 1,
           family  = "serif",
           cex.lab = 1.2,
           cex.axis= 1.2,
           mgp     = c(2.2,1,0),
           oma     = c(2,2,0,0),
           mar     = c(3.1,4,2.3,.5),
           mfrow   = c(2,2))

XLAB = expression(paste('1/'*k[b]*'(1/'*T[0]*' - 1/'*T[opt]*') or '*theta*' (1/eV)'))
YLAB = expression(paste('Size corrected Ln '*italic(mu)[m]*' ('*d^-1*')') )
YLAB1= expression(paste(T[opt]*' corrected Ln '*italic(mu)[m]*' ('*d^-1*')') )
#Take into account size & species using lme

PHYnls$Species <- as.factor(PHYnls$Species)
PHYnls<- PHYnls[!is.na(PHYnls$Lnum),]
Ym0   <- lm(Lnum  ~ T.K(Topt), PHYnls) 
Ym1   <- lm(Lnum  ~ T.K(Topt)+ lnVolume, PHYnls)
Ymer1 <- lmer(Lnum  ~ T.K(Topt) + (1|Species),  data=PHYnls)
Ymer2 <- lmer(Lnum  ~ T.K(Topt) + (1|Group/Species),  data=PHYnls)
Ymer3 <- lmer(Lnum  ~ T.K(Topt) + lnVolume + (1|Group), data=PHYnls)
Ymer4 <- lmer(Lnum  ~ T.K(Topt) + lnVolume + (1|Group/Species), data=PHYnls) #Best
Ymer5 <- lmer(Lnum  ~ T.K(Topt) + lnVolume + (1|Species), data=PHYnls)
AIC(Ym0, Ym1,  Ymer3,Ymer4, Ymer5)

#Obtain coefficients
CF = coef(summary(Ymer4))

#Size coefficient
alpha = CF[3,1]
beta  = CF[2,1]
intp  = CF[1,1]
Ei    = round(beta,2)
seEi  = round(CF[2,2],2)
sealpha = round(CF[3,2],2)

#Calculate size-corrected Lnum
PHYnls$corrLnum = PHYnls$Lnum - alpha*PHYnls$lnVolume

#Calcualte Topt-corrected Lnum
PHYnls$corrLnum2 = PHYnls$Lnum - beta*T.K(PHYnls$Topt)

plot(PHYnls$lnVolume, PHYnls$corrLnum2, 
     pch = 16, cex = 0.5,
     xlim = c(-5, 17),
     ylim = c(log(0.02), 2),
     xlab = expression(paste('Ln volume ('*µm^3*')')),
     ylab = YLAB1,
     type = 'n')
txt = bquote('a) Phyto, ' 
             * italic(alpha[µm]) * ' = ' 
             * .(round(alpha, 2)) * '±' 
             * .(sealpha))
mtext(txt, adj = 0, cex = 1.2)

Groups = unique(PHYnls$Group)
for (i in 1:length(unique(PHYnls$Group))){
    tmp = PHYnls[PHYnls$Group == Groups[i], ]
    points(tmp$lnVolume, tmp$corrLnum2, 
           pch = 16, cex = 0.5, col = i)
}

abline(intp, alpha)

plot(T.K(PHYnls$Topt), PHYnls$corrLnum, 
     pch = 16, 
     cex = 0.5,
    xlim = T.K(c(-5,40)),
    ylim = c(log(0.1), 3.1),
    xlab = XLAB,
    ylab = YLAB,
    type = 'n')
txt = bquote('b) Phyto, '
             * italic(E[i])*' = '
             * .(Ei) * '±' * .(seEi)* ' eV' )

mtext(txt, adj=0, cex=1.2)
Groups = unique(PHYnls$Group)
for (i in 1:length(unique(PHYnls$Group))){
    tmp = PHYnls[PHYnls$Group == Groups[i], ]
    points(T.K(tmp$Topt), tmp$corrLnum, 
           pch = 16, cex = 0.5, col = i)
}

abline(intp, Ei)
#x1   <- seq(0,30,10)
#y1   <- T.K(x1)
#x1[length(x1)]=paste0(x1[length(x1)],' ºC')
#axis(3, at= y1, label = x1)

legend('topleft',pch=16, cex=1.2,
        col=1:5, legend=as.character(Groups))

ZOOnls$ivTopt <- T.K(ZOOnls$Topt)
ZOOnls        <- ZOOnls[!is.na(ZOOnls$Lnum),]
Ym1  <-  lm(Lnum~ ivTopt, ZOOnls)  
Ym2  <-  lm(Lnum~ ivTopt + lnVolume, ZOOnls)
#Yme1 <- lmer(Lnum~ ivTopt + (1|Species), ZOOnls)  #Better
#Yme2 <- lmer(Lnum~ ivTopt + (1|Group/Species), ZOOnls)  
Yme3 <- lmer(Lnum~ ivTopt + lnVolume + (1|Group), ZOOnls)
Yme4 <- lmer(Lnum~ ivTopt + lnVolume + (1|Group/Species), ZOOnls)  
Yme5 <- lmer(Lnum~ ivTopt + lnVolume + (1|Species), ZOOnls) #best
AIC(Ym1, Ym2, Yme3, Yme4, Yme5) 

#Obtain coefficients
CF = coef(summary(Yme5))

#Size coefficient
alpha = CF[3,1]
beta  = CF[2,1]
intp  = CF[1,1]
sealpha = round(CF[3,2],2)

#Calcualte Topt-corrected Lnum
ZOOnls$corrLnum2 = ZOOnls$Lnum - beta*T.K(ZOOnls$Topt)

#Calculate size-corrected Lnum
ZOOnls$corrLnum = ZOOnls$Lnum - alpha*ZOOnls$lnVolume

plot(ZOOnls$lnVolume, ZOOnls$corrLnum2, 
     pch = 16, cex = 0.5,
     xlim = c(-5, 17),
     ylim = c(-3, 2),
     xlab = expression(paste('Ln volume ('*µm^3*')')),
     ylab = YLAB1,
     type = 'n')
txt = bquote('c) Microzoo, ' 
             * italic(alpha[µm]) * ' = ' 
             * .(round(alpha, 2)) * '±' 
             * .(sealpha))
mtext(txt, adj=0, cex=1.2)

Groups = unique(ZOOnls$Group)
for (i in 1:length(unique(ZOOnls$Group))){
    tmp = ZOOnls[ZOOnls$Group == Groups[i], ]
    points(tmp$lnVolume, tmp$corrLnum2, 
           pch = 16, cex = 0.5, col = i)
}

abline(intp, alpha)

Ei  = round(CF[2,1],2)
seEi= round(CF[2,2],2)
txt = bquote('d) Microzoo, '*italic(E[i])*' = '* .(Ei) * '±' * .(seEi)* ' eV' )

plot(ZOOnls$ivTopt, ZOOnls$corrLnum, pch = 16, cex = 0.5,
     xlim = T.K(c(-5,40)),
     ylim = c(log(1), 3.1),
     xlab = XLAB,
     ylab = YLAB,
     type = 'n')

mtext(txt, adj=0, cex=1.2)
Groups = unique(ZOOnls$Group)
for (i in 1:length(unique(ZOOnls$Group))){
    tmp = ZOOnls[ZOOnls$Group == Groups[i], ]
    points(T.K(tmp$Topt), tmp$corrLnum, 
           pch = 16, cex = 0.5, col = i)
}
legend('topleft',pch=16, cex=1.2,
       col=1:5, legend=as.character(Groups))
abline(lm(corrLnum ~ ivTopt, ZOOnls))

dev.off()

#Plot size vs. Topt
# pdf('Fig6PHYZOOTopt_vol.pdf', height = 4, width = 8)
# op <- par(font.lab = 1,
#           family  = "serif",
#           cex.lab = 1.2,
#           cex.axis= 1.2,
#           mgp     = c(2.2,1,0),
#           oma     = c(2,2,0,0),
#           mar     = c(3,3,2.3,.5),
#           mfrow   = c(1,2))
# Groups = unique(PHYnls$Group)
# plot(PHYnls$Topt, PHYnls$lnVolume, 
#      xlim = c(-5, 40),
#      ylim = c(-10, 15.5),
#      pch  = 16, cex = 0.5,
#      xlab = '',
#      ylab = '',
#      type = 'n')
# mtext('a) Phytoplankton',  cex=1.2, adj=0)
# for (i in 1:length(unique(PHYnls$Group))){
#     tmp = PHYnls[PHYnls$Group == Groups[i], ]
#     points(tmp$Topt, tmp$lnVolume, 
#            pch = 16, cex = 0.5, col = i)
# }
# lm1 <- lm(lnVolume ~ Topt, PHYnls)
# abline(lm1)
# legend('bottomleft',pch=16,
#        col=1:5, legend=as.character(Groups))
# 
# Groups = unique(ZOOnls$Group)
# plot(ZOOnls$Topt, ZOOnls$lnVolume,  
#      pch  = 16, cex = 0.5,
#      xlim = c(-5, 40),
#      ylim = c(-10, 15.5),
#      xlab = '',
#      ylab = '',
#      type = 'n')
# for (i in 1:length(unique(ZOOnls$Group))){
#     tmp = ZOOnls[ZOOnls$Group == Groups[i], ]
#     points(tmp$Topt, tmp$lnVolume, 
#            pch = 16, cex = 0.5, col = i)
# }
# lm1 <- lm(lnVolume ~ Topt, ZOOnls)
# abline(lm1)
# legend('bottomright',pch=16, cex=1.2,
#        col=1:5, legend=as.character(Groups))
# mtext('b) Microzooplankton', cex=1.2, adj=0)
# XLAB = expression(paste(T[opt]*' (ºC)'))
# YLAB = expression(paste('Ln volume ('*µm^3*')') )  
# mtext(XLAB, side = 1, outer=T, adj = .5, cex=1.2)
# mtext(YLAB, side = 2, outer=T, adj = .5, cex=1.2)
# dev.off()

#Check Pro and Syn
# Prodat <- PHYnls[PHYnls$Genus == 'Prochlorococcus', ]
# ProLM1 <- lm(Lnum ~ T.K(Topt), Prodat)
# Syndat <- PHYnls[PHYnls$Genus == 'Synechococcus', ]
# SynLM1 <- lm(Lnum ~ T.K(Topt), Syndat)

# pdf('ProSyn.pdf', height = 4, width = 4)
# op <- par(font.lab = 1,
#           family  = "serif",
#           cex.lab = 1.2,
#           cex.axis= 1.2,
#           mgp     = c(2.2,1,0),
#           oma     = c(2,2,0,0),
#           mar     = c(3,3,2.3,.5),
#           mfrow   = c(1,1))
# plot(T.K(Prodat$Topt), Prodat$Lnum, 
#      xlim = T.K(c(15,30)),
#      ylim = log(c(0.2, 2)),
#      xlab = 'Normalized optimal temperature',
#      ylab = 'Ln maximal growth rate (d-1)',
#      pch  = 16,
#      col  = 'blue')
# abline(ProLM1, col = 'blue')
# points(T.K(Syndat$Topt), Syndat$Lnum, pch = 16, col=2)
# abline(SynLM1, col = 2)
# legend('bottomleft', c('Pro','Syn'), lwd=1, col=c('blue','red'))
# dev.off()
