library(here)
source('carbonate.r')

alks  <- seq(2100,2500,length.out=3)
temps <- seq(0,30,length.out=50)
DICs  <- seq(1800,2300,length.out=3)
#phs   <- seq(4,9,length.out=3)

n <- length(alks)*length(temps)*length(DICs)

PCO2 <- array(NA, dim=c(length(alks),length(temps),length(DICs)))

for(i in 1:length(alks)){
  for(j in 1:length(temps)){
    for(k in 1:length(DICs)){
         PCO2[i,j,k] <- carbonate(TEMP=temps[j],alk=alks[i],DIC=DICs[k])$pCO2        
      }
    }
  }

pdf('slopes.pdf',height=4,width=8)
par(mfrow=c(1,2),mar=c(2,2,2,2),oma=c(2,2,2,2))
plot(temps,log(PCO2[2,,2]),type='l',ylim=c(4,9),ylab='',xlab='')
  mtext(side=2,line=2.5,expression('log'~italic(p)*'CO'['2']))
  mtext(side=1,line=2.5,expression('Temperature ['*degree*'C]'))
lines(temps,log(PCO2[1,,2]),lty=2)
lines(temps,log(PCO2[3,,2]),lty=3)
legend('topleft', lty=c(2,1,3), legend=c('ALK = 2100 [umol/kg]', 'ALK = 2300','ALK = 2500'),cex=0.6,bty='n')
text(25,8,paste(round(summary(lm(log(PCO2[1,,2])~temps))$coefficients[2,1],5)))
text(25,6.75,paste(round(summary(lm(log(PCO2[2,,2])~temps))$coefficients[2,1],5)))
text(25,5,paste(round(summary(lm(log(PCO2[3,,2])~temps))$coefficients[2,1],5)))

plot(temps,log(PCO2[2,,2]),type='l',ylim=c(4,9),ylab='',xlab='')
  mtext(side=1,line=2.5,expression('Temperature ['*degree*'C]'))
lines(temps,log(PCO2[2,,1]),lty=2)
lines(temps,log(PCO2[2,,3]),lty=3)
legend('topleft', lty=c(2,1,3), legend=c('DIC = 1800 [umol/kg]', 'DIC = 2050','DIC = 2300'),cex=0.6,bty='n')
text(25,8.25,paste(round(summary(lm(log(PCO2[2,,1])~temps))$coefficients[2,1],5)))
text(25,6.75,paste(round(summary(lm(log(PCO2[2,,2])~temps))$coefficients[2,1],5)))
text(25,4.5,paste(round(summary(lm(log(PCO2[2,,3])~temps))$coefficients[2,1],5)))

mtext(side=3,outer=TRUE,'Chen et al. 2007 slope: 0.0408',line=-1)
dev.off()

