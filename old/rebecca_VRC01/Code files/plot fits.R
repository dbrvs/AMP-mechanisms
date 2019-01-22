library(ggplot2)
source("fitting Ab curve.R")
source("finding slopes.R")

#IMPORTING PATIENT VIRAL LOAD DATA
patient20 = make20()
patient21 = make21()
patient22 = make22()
patient23 = make23()
patient24 = make24()
patient25 = make25()
patient26 = make26()
patient27 = make27()

#IMPORTING PATIENT ANTIBODY CONCENTRATION DATA
ab20 = get20()
ab21 = get21()
ab22 = get22()
ab23 = get23()
ab24 = get24()
ab25 = get25()
ab26 = get26()
ab27 = get27()

#GETTING ANTIBODY FITS FOR EACH PATIENT
decay20 = getfit(ab20, 20)
decay21 = getfit(ab21, 21)
decay22 = getfit(ab22, 22)
decay23 = getfit(ab23, 23)
decay24 = getfit(ab24, 24)
decay25 = getfit(ab25, 25)
decay26 = getfit(ab26, 26)
decay27 = getfit(ab27, 27)

#GETTING ADCC ONLY FITS FOR EACH PATIENT
source("fit only ADCC.R")
par(mfrow=c(1,1))
fit20ADCC = completeADCCfit20(patient20, decay20)
fit21ADCC = completeADCCfit21(patient21, decay21)
fit22ADCC = completeADCCfit22(patient22, decay22)
fit23ADCC = completeADCCfit23(patient23, decay23)
fit24ADCC = completeADCCfit24(patient24, decay24)
#fit25ADCC = completeADCCfit25(patient25, decay25)
fit26ADCC = completeADCCfit26(patient26, decay26)
fit27ADCC = completeADCCfit27(patient27, decay27)

#GETTING NAB ONLY FITS FOR EACH PATIENT
source("fit only NAB.R")
par(mfrow=c(1,1))
fit20NAB = completeNABfit20(patient20, decay20)
fit21NAB = completeNABfit21(patient21, decay21)
fit22NAB = completeNABfit22(patient22, decay22)
fit23NAB = completeNABfit23(patient23, decay23)
fit24NAB = completeNABfit24(patient24, decay24)
#fit25NAB = completeNABfit25(patient25, decay25)
fit26NAB = completeNABfit26(patient26, decay26)
fit27NAB = completeNABfit27(patient27, decay27)

#GETTING BOTH ADDC AND NAB FITS FOR EACH PATIENT
source("fit both NAB and ADCC.R")
par(mfrow=c(1,1))
fit20BOTH = completebothfit20(patient20, decay20)
fit21BOTH = completebothfit21(patient21, decay21)
fit22BOTH = completebothfit22(patient22, decay22)
fit23BOTH = completebothfit23(patient23, decay23)
fit24BOTH = completebothfit24(patient24, decay24)
#fit25BOTH = completebothfit25(patient25, decay25)
fit26BOTH = completebothfit26(patient26, decay26)
fit27BOTH = completebothfit27(patient27, decay27)


########################################################
#This part of code is intended to be used to graph fits
#After getting all of the fits (above), you could then use the following functions
#to get data sets from the fits, and could then graph them (code to actually plot
#the data needs to written)

#THESE FUNCTIONS GET S,I,V OVER TIME FOR DIFFERENT SCENARIOS
getmodelADCC = function(fit, decay, known){
  parameters = fit$solution
  alpha = unname(parameters[1])
  deltas = unname(parameters[2])
  pi = unname(parameters[3])
  EC50 = unname(parameters[4])
  h = unname(parameters[5])
  V = unname(known[1])
  deltai = unname(known[2])
  beta = unname(known[3])
  S=alpha/(deltas+beta*unname(V))
  I=beta*unname(S)*unname(V)/deltai
  c = pi*I/V
  #browser()
  state = c(S=S, I=I, V=V)
  model = as.data.frame(ode(y=state, times=times, func=ODEADCC, 
                            parms=c(alpha=alpha,deltas=deltas,beta=beta,deltai=deltai,
                                    pi=pi,c=c,EC50=EC50,h=h), 
                            decay=decay, maxsteps=1e6))
}
getmodelNAB = function(fit, decay, known){
  parameters = fit$solution
  alpha = unname(parameters[1])
  deltas = unname(parameters[2])
  pi = unname(parameters[3])
  IC50 = unname(parameters[4])
  m = unname(parameters[5])
  V = unname(known[1])
  deltai = unname(known[2])
  beta = unname(known[3])
  S=alpha/(deltas+beta*unname(V))
  I=beta*unname(S)*unname(V)/deltai
  c = pi*I/V
  state = c(S=S, I=I, V=V)
  model = as.data.frame(ode(y=state, times=times, func=ODENAB, 
                            parms=c(alpha=alpha,deltas=deltas,beta=beta,deltai=deltai,
                                    pi=pi,c=c,IC50=IC50,m=m), 
                            decay=decay, maxsteps=1e6))
}
getmodelboth = function(fit, decay, known){
  parameters = fit$solution
  alpha = unname(parameters[1])
  deltas = unname(parameters[2])
  pi = unname(parameters[3])
  IC50 = unname(parameters[4])
  EC50 = unname(parameters[5])
  m = unname(parameters[6])
  h = unname(parameters[7])
  V = unname(known[1])
  deltai = unname(known[2])
  beta = unname(known[3])
  S=alpha/(deltas+beta*unname(V))
  I=beta*unname(S)*unname(V)/deltai
  c = pi*I/V
  state = c(S=S, I=I, V=V)
  model = as.data.frame(ode(y=state, times=times, func=ODEboth, 
                            parms=c(alpha=alpha,deltas=deltas,beta=beta,deltai=deltai,
                                    pi=pi,c=c,IC50=IC50,EC50=EC50,m=m,h=h), 
                            decay=decay, maxsteps=1e6))
}

#DIFF EQ'S FOR DIFFERENT SCENARIOS
ODEADCC = function(times, state, parms, decay) {
  with(as.list(c(state, parms)), {
    if(times<0){
      conc=0
    }
    else{
      r1 = coef(decay)[1]
      r2 = coef(decay)[2]
      a = coef(decay)[3]
      b = coef(decay)[4]
      conc = a*exp(r1*times)+a*b*exp(r2*times)
    }
    fu = 1
    IAP = log10(1+(conc/EC50)^h)
    fa = 1-1/(10^IAP)
    if(V<.05 || is.na(V)){
      V=.05
    }
    dS = alpha - deltas*S -fu*beta*S*V
    dI = fu*beta*S*V - (deltai+fa)*I 
    dV = pi*I - c*V
    list(c(dS, dI, dV))
  })
}
ODENAB = function(times, state, parms, decay) {
  with(as.list(c(state, parms)), {
    if(times<0){
      conc=0
    }
    else{
      r1 = coef(decay)[1]
      r2 = coef(decay)[2]
      a = coef(decay)[3]
      b = coef(decay)[4]
      conc = a*exp(r1*times)+a*b*exp(r2*times)
    }
    IIP = log10(1+(conc/IC50)^m)
    fu = 1/(10^IIP)
    fa = 0
    if(V<.05 || is.na(V)){
      V=.05
    }
    dS = alpha - deltas*S -fu*beta*S*V
    dI = fu*beta*S*V - (deltai+fa)*I 
    dV = pi*I - c*V
    #browser()
    list(c(dS, dI, dV))
  })
}
ODEboth = function(times, state, parms, decay) {
  with(as.list(c(state, parms)), {
    if(times<0){
      conc=0
    }
    else{
      r1 = coef(decay)[1]
      r2 = coef(decay)[2]
      a = coef(decay)[3]
      b = coef(decay)[4]
      conc = a*exp(r1*times)+a*b*exp(r2*times)
    }
    IIP = log10(1+(conc/IC50)^m)
    fu = 1/(10^IIP)
    IAP = log10(1+(conc/EC50)^h)
    fa = 1-1/(10^IAP)
    if(V<.05 || is.na(V)){
      V=.05
    }
    dS = alpha - deltas*S -fu*beta*S*V
    dI = fu*beta*S*V - (deltai+fa)*I 
    dV = pi*I - c*V
    #browser()
    list(c(dS, dI, dV))
  })
}
