#This file executes all of the necessary calculations to examine HIV and the effective of the antibody over time

parameters = c(alpha=.15, deltas=.01, beta=6.5*10^-4, deltai=.39, pi=850, c=3) #Stafford HIV parameters

#this function returns a matrix with the antibody data
getdatAb = function(times, timestep, initial, infusionday) {
  
  infusionv = seq(times[1], infusionday-timestep, timestep)
  afterv = seq(infusionday, endtime, timestep)
  
  conc = NULL
  for(i in 1:length(infusionv)/timestep) {
    conc[i] = 0
  }
  for(i in 1:length(afterv)) {
    #these parameters are averages from the parameters of the the biphasic decay fits of patient data
    conc[i+(infusionday-starttime)/timestep] = 1050*exp(-1.1*(afterv[i]-infusionday)) + 1050*.414*exp(-.1*(afterv[i]-infusionday))
  }
  IC50 = .33 #mcg/mL
  #slope data from http://www.nature.com/ncomms/2015/150929/ncomms9443/full/ncomms9443.html#supplementary-information
  slopes = c(1.14, 1.56, 1.23, 1.5, 1.15, 1.64, 1.44, 1.22,
             1.53, 1.43, 1.47, 1.34, 1.48, 1.37, 0.71, 1.06)
  
  aveslope = mean(slopes) #here I make the assumption that the HIV envelopes are all equally prevalent
  IIP =log10(1+(conc/IC50)^aveslope) #instantaneous inhibitory potential
  fu = 1/(10^(IIP)) #fraction of infections NOT prevented by NAb
  
  #from Bruel
  h = .5846 #on a graph of %adcc per mL
  EC50 = 2.3 #mcg/mL (for ADCC killing)
  
  IAP = log10(1+(conc/EC50)^h) #instantaneous adcc potential
  fa = 1-1/(10^IAP) #fraction of infected cells killed by ADCC
  datAb = data.frame(times, conc, IIP, fu, IAP, fa)
}
#this function returns a matrix with HIV data
getdatHIV = function(times, timestep, fu, fa) {
  library(deSolve)
  state = c(S=15, I=.0001, V=10^-6)
  ODE = function(time, state, parameters) {
    with(as.list(c(state, parameters)), {
      dS = alpha - deltas*S - fu[(time-starttime)/timestep+1]*beta*S*V
      dI = fu[(time-starttime)/timestep+1]*beta*S*V - (deltai+fa[(time-starttime)/timestep+1])*I 
      dV = pi*I - c*V
      list(c(dS, dI, dV))
    })
  }
  datHIV = as.data.frame(ode(atol=0, y=state, times=times, func=ODE, parms=parameters))
}
getR = function(fu, fa, S) {
  beta = parameters[[3]]
  deltai = parameters[[4]]
  pi = parameters[[5]]
  c = parameters[[6]]
  
  R = fu*beta*pi*S/c/(deltai+fa)
}

