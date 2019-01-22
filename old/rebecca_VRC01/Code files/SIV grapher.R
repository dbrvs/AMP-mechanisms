library(ggplot2)
source("getdata SIV.R")

#term glossary:
#fa = fraction of infected cells affected/killed by ADCC
#fu = fraction of infection events unaffected by NAB
#IIP = instantaneous inhibitory potential (measure of NAB activity)
#IAP = instantaneous ADCC potential (measure of ADCC activity)
#Rpot = R potential

#setting up starting conditions
initial = 0 #mcg/mL
starttime = -50
timestep = 1
endtime = 300
infusionday = 0
times = seq(starttime, endtime, timestep)

#getting data from other file
nofa = fa=numeric(length(times)) #vector of 0's, used to simulate NO ADCC
nofu = fu=seq(from=1, to=1, length.out=length(times)) #vector of 1's, used to simulate NO NAB
datAb = getdatAb(times, timestep, initial, infusionday)

datHIVboth = getdatHIV(times, timestep, datAb$fu, datAb$fa) #both mechanisms scenario
datHIVNAB = getdatHIV(times, timestep, datAb$fu, nofa) #NAB only (no ADCC) scenario
datHIVADCC = getdatHIV(times, timestep, nofu, datAb$fa) #ADCC only (no NAB) scenario
Rboth = getR(fu=datAb$fu, fa=datAb$fa, S=datHIVboth$S) #both mechanisms scenario
RNAB = getR(fu=datAb$fu, nofa, S=datHIVNAB$S) #NAB only (no ADCC) scenario
RADCC = getR(fu=nofu, fa=datAb$fa, S=datHIVADCC$S) #ADCC only (no NAB) scenario
Rpot = getR(fu=nofu, fa=nofa, S=datHIVboth$S) #R potential is a measure of the effective reproductive number, 
                                              #if there were no antibody present

#condensing all data into one data frame
alldata = data.frame(datAb, S=datHIVboth$S, I=datHIVboth$I, V=datHIVboth$V,
                     SNAB=datHIVNAB$S, INAB=datHIVNAB$I, VNAB=datHIVNAB$V,
                     SADCC=datHIVADCC$S, IADCC=datHIVADCC$I, VADCC=datHIVADCC$V,
                     Rboth=Rboth, RNAB=RNAB, RADCC=RADCC, Rpot)

#Graphing susceptible
plot(alldata$times, log10(alldata$S), type = 'l', xlab="time (days)", ylab="log10(cells/uL)",
     ylim=c(-4,5))
title("Susceptible")


#Graphing infected
plot(alldata$times, log10(alldata$I), type = 'l', xlab="time (days)", ylab="log10(cells/uL)", 
     ylim=c(-4,5))
title("Infected")

#Graphing total cell density
celltot = alldata$S + alldata$I
plot(alldata$times, log10(celltot), type = 'l', xlab="time (days)", ylab="log10(cells/uL)",
     ylim=c(-4,5))
title("Total cell density")

#Graphing IIP and Fu
plot(alldata$times, alldata$IIP, type='l', col=2, xlab="time (days)", ylab="IIP")
par(new=T)
plot(alldata$times, log10(alldata$fu), axes=F, xlab="", ylab="", type='l', ylim=c(-5,3), col=4)
axis(4)
mtext("log 10 Fu (fraction unaffected)", side=4, line=2.5)
title("Ab free virus neutralizing potency after 1000 mcg/mL infusion on day 50")
legend(x=100, y=2, "IIP = red, Fu = blue")

#Graphing virus and Ab concentration
par(mar=c(3,4,4,5))
plot(alldata$times, log10(alldata$V*1000), xlab="days after VRC01 infusion", 
     ylab="log 10 viral load", type='l', col=4, ylim=c(-4,7))
par(new=T)
plot(alldata$times, log10(alldata$conc), type='l', axes=F, col=2, xlab="", ylab="")
axis(4, ylim=c(-5,5))
mtext("log 10 Ab concentration", side=4, line=2.5)
title("")
legend(x=-50, y=-1, "Ab = red, Virus = blue")

#Plotting "detectable" levels and patient data
source("finding slopes.R")
patient20 = make20()
patient21 = make21()
patient22 = make22()
patient23 = make23()
patient24 = make24()
patient25 = make25()
patient26 = make26()
patient27 = make27()
par(mfrow=c(1,1))
detectable = log10(alldata$V*1000)
for(i in 1:length(detectable)){
  if(detectable[i]<1){
    detectable[i]=1 
  }
}
plot(alldata$times, detectable, type='l', col=1, xlab="days after VRC01 infusion",
     ylab="log 10 viral load", ylim=c(1, 7), xlim=c(-20,300), lwd=2)
par(new=T)
plot(alldata$times, log10(alldata$conc), type='l', axes=F, lty=2, xlab="", ylab="", ylim=c(-3,4), xlim=c(-20,300))
axis(4)
mtext("log 10 Ab concentration", side=4, line=2.5)
title("")
legend(x=50, y=4, c("Virus","Ab"), lty=c(1,2), col=c(1,1), cex=.8)
par(new=T)
plot(patient20$studytime, log10(patient20$viral), axes=F, type='o', 
     xlab="", ylab="", xlim=c(-20,300), ylim=c(1,6), lwd=.5, pch=20, col="blue")
par(new=T)
plot(patient21$studytime, log10(patient21$viral), axes=F, type='o', 
     xlab="", ylab="", xlim=c(-20,300), ylim=c(1,6), lwd=.5, pch=20, col="red")
par(new=T)
plot(patient22$studytime, log10(patient22$viral), axes=F, type='o', 
     xlab="", ylab="", xlim=c(-20,300), ylim=c(1,6), lwd=.5, pch=20, col="green")
par(new=T)
plot(patient23$studytime, log10(patient23$viral), axes=F, type='o', 
     xlab="", ylab="", xlim=c(-20,300), ylim=c(1,6), lwd=.5, pch=20, col="purple")
par(new=T)
plot(patient24$studytime, log10(patient24$viral), axes=F, type='o', 
     xlab="", ylab="", xlim=c(-20,300), ylim=c(1,6), lwd=.5, pch=20, col="orange")
par(new=T)
plot(patient25$studytime, log10(patient25$viral), axes=F, type='o', 
     xlab="", ylab="", xlim=c(-20,300), ylim=c(1,6), lwd=.5, pch=20, col="brown")
par(new=T)
plot(patient26$studytime, log10(patient26$viral), axes=F, type='o', 
     xlab="", ylab="", xlim=c(-20,300), ylim=c(1,6), lwd=.5, pch=20, col="pink")
par(new=T)
plot(patient27$studytime, log10(patient27$viral), axes=F, type='o', 
     xlab="", ylab="", xlim=c(-20,300), ylim=c(1,6), lwd=.5, pch=20, col="yellow")

#Graphing viral loads for different scenarios
par(mar=c(5,4,4,4))
plot(alldata$time, log10(alldata$V*1000), type='l', col=4, 
     xlab="days after VRC01 infusion", ylab="log 10 viral load", lwd='1.5', ylim=c(1, 8), xlim=c(-10,20))
lines(alldata$times, log10(alldata$VNAB*1000), type='l', col=2, lwd='1.5', ylim=c(1, 8), xlim=c(-10,20))
lines(alldata$times, log10(alldata$VADCC*1000), type='l', col=3, lwd='1.5', ylim=c(1, 8), xlim=c(-10,20))
legend(x=-5, y=4, c("NAb and ADCC","NAb only (or ART only)","ADCC only"),
       col=c(4,2,3), lty=c(1,1,1), bty='n', cex=.8)

#Graphing R for different scenarios
plot(alldata$times, (alldata$Rboth), type='l', col=4, xlim=c(-25,100), ylim=c(0,8),
     lty=5, lwd='1.5', ylab="effective reproductive number", xlab="days after VRC01 infusion")
lines(alldata$times, (alldata$RNAB), lty=5, col=2, lwd='1.5', xlim=c(-25,100))
lines(alldata$times, (alldata$RADCC), lty=5, col=3, lwd='1.5', xlim=c(-25,100))
legend(x=-10, y=7, c("NAb and ADCC","NAb only","ADCC only"), col=c(4,2,3), lty=c(2,2,2), bty='n', cex=.8)

#comparing r with r potential
plot(alldata$times, log10(alldata$Rboth), type='l', col=2, ylab="log 10 viruses produced per infected cell", xlab="days")
lines(alldata$times, log10(alldata$Rpot), col=4)
abline(lty="dashed", h=0)  
title("Effective reproductive number vs potential reproductive number")  
  
  
  
  

