#this file plots dose-response curves given potency parameters

conc = seq(0,1000,.001)

#ADCC POTENCY MEASUREMENTS
EC50 = 2.3
h = .5846
ADCC = 100-100/(1+(conc/EC50)^h)

#NAB POTENCY MEASUREMENTS
IC50 = .33
m = 1.329
NAb = 100-100/(1+(conc/IC50)^m)

plot(log10(conc), ADCC, type='l', ylab= "% effect", xlab = "log 10 conc", col='green')
lines(log10(conc), NAb, col=2)
legend(x=-2.5, y=60, c("Neutralization","ADCC"), col=c(2,3), bty='n', cex=.8, lty=c(1,1))

