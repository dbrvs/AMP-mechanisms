#data from Bruel

#For NLAD8
#EC50=2.3 mcg/mL
concs1 = c(0.3189788489448251, 0.36366015559943676, 0.3970640717565601) #log10(mcl/ml)
ADCCs1 = c(log10(82.27927006257289), log10(86.87016725750038), log10(91.45249934572098)) #log10(%adcc/day)
plot(concs1, ADCCs1)

linfit1 = lm(ADCCs1 ~ concs1)
lines(concs1, predict(linfit1))
linfit1

#For NL4.3
#EC50=2.1 mcg/mL
concs2 = c(0.3263157894736839, 0.3578947368421055, 0.2842105263157899)
ADCCs2 = c(log10(152.39380804953566), log10(157.34984520123845), log10(145.31517027863782))
plot(concs2, ADCCs2)

linfit2 = lm(ADCCs2 ~ concs2)
lines(concs2, predict(linfit2))
linfit2

##GRAPHING CONC AND FRACTION AFFECTED BY ADCC OVER CONCENTRATIONS (for NLAD8)
EC50 = 2.3 #from NLAD8
m = linfit1$coefficients[2] #hill slope for NLAD8
concs = seq(0, 1000, .01)
IAP = log10(1+(concs/EC50)^m) #instantaneous ADCC potential
fa = (1-1/(10^IAP))
plot(log10(concs), fa, type='l')

