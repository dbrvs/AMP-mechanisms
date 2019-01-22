library(deSolve)
library(minpack.lm)
library(nloptr)
source("fitting Ab curve.R")
source("finding slopes.R")

#geting viral load data
patient20 = make20()
patient21 = make21()
patient22 = make22()
patient23 = make23()
patient24 = make24()
patient25 = make25()
patient26 = make26()
patient27 = make27()
plotall(patient20, patient21, patient22, patient23, 
        patient24, patient25, patient26, patient27)

#getting antibody fits
ab20 = get20()
ab21 = get21()
ab22 = get22()
ab23 = get23()
ab24 = get24()
ab25 = get25()
ab26 = get26()
ab27 = get27()
decay20 = getfit(ab20, 20)
decay21 = getfit(ab21, 21)
decay22 = getfit(ab22, 22)
decay23 = getfit(ab23, 23)
decay24 = getfit(ab24, 24)
decay25 = getfit(ab25, 25)
decay26 = getfit(ab26, 26)
decay27 = getfit(ab27, 27)

#fitting model to each patient
par(mfrow=c(1,1))

completeADCCfit20 = function(patient20, decay20){
  fit20_onlyADCC = fitmodel(patient20, decay20, 
                            parameters = c(alpha = (.8), deltas=(.009), pi=(2500), EC50=(100), h=(20)),
                            V=3.547, low = c(.1, .009, 200, 0, 1), 
                            high = c(.9, .1, 5000, 500, 100)) 
  fit20_onlyADCC2 = fitmodel(patient20, decay20, parameters= fit20_onlyADCC$solution, 
                             V=3.547, low = c(.1, .009, 200, 0, 1), 
                             high = c(.9, .1, 5000, 500, 100)) 
  fit20_onlyADCCfinal = finetune(patient20, decay20, parameters= fit20_onlyADCC2$solution, 
                                 V=3.547, low = c(.1, .009, 200, 0, 1), 
                                 high = c(.9, .1, 5000, 500, 100)) 
  return(fit20_onlyADCCfinal)
}
completeADCCfit21 = function(patient20, decay20){
  fit21_onlyADCC = fitmodel(patient21, decay21, 
                            parameters = c(alpha = (.3), deltas=(.009), pi=(1000), EC50=(1000), h=(50)),
                            V=6.551, low = c(.001, 0.001, 1000, 1000, 40), 
                            high = c(2, .01, 10000, 1500, 100)) 
  fit21_onlyADCC2 = fitmodel(patient21, decay21, parameters= fit21_onlyADCC$solution, 
                             V=6.551, low = c(.001, 0.001, 1000, 1000, 40), 
                             high = c(2, .01, 10000, 1500, 100))
  fit21_onlyADCCfinal = finetune(patient21, decay21, parameters= fit21_onlyADCC2$solution, 
                                 V=6.551, low = c(.001, 0.001, 1000, 1000, 40), 
                                 high = c(2, .01, 10000, 1500, 100)) 
  return(fit21_onlyADCCfinal)
}
completeADCCfit22 = function(patient22, decay22){
  fit22_onlyADCC = fitmodel(patient22, decay22, 
                            parameters = c(alpha = (.15), deltas=(.01), pi=(1200), EC50=(100), h=(25)),
                            V=.745, low = c(.1, .0001, 500, 20, 15), 
                            high = c(3, 1, 5000, 400, 40)) 
  fit22_onlyADCC2 = fitmodel(patient22, decay22, parameters= fit22_onlyADCC$solution, 
                             V=.745, low = c(.1, .0001, 500, 20, 15), 
                             high = c(3, 1, 5000, 400, 40)) 
  fit22_onlyADCCfinal = finetune(patient22, decay22, parameters= fit22_onlyADCC2$solution, 
                                 V=.745, low = c(.1, .0001, 500, 20, 15), 
                                 high = c(3, 1, 5000, 400, 40))  
  return(fit22_onlyADCCfinal)
}
completeADCCfit23 = function(patient20, decay20){
  fit23_onlyADCC = fitmodel(patient23, decay23, 
                            parameters = c(alpha = (.15), deltas=(.001), pi=(330), EC50=(120), h=(12)),
                            V=27.894, low = c(.05, .0005, 100, 20, .1), 
                            high = c(1, 1, 500, 250, 17)) 
  fit23_onlyADCC2 = fitmodel(patient23, decay23, parameters= fit23_onlyADCC$solution, 
                             V=27.894, low = c(.05, .0005, 100, 20, .1), 
                             high = c(1, 1, 500, 250, 17)) 
  fit23_onlyADCCfinal = finetune(patient23, decay23, parameters= fit23_onlyADCC2$solution, 
                                 V=27.894, low = c(.05, .0005, 100, 20, .1), 
                                 high = c(1, 1, 500, 250, 17))  
  return(fit23_onlyADCCfinal)
}
completeADCCfit24 = function(patient24, decay24){
  fit24_onlyADCC = fitmodel(patient24, decay24, 
                            parameters = c(alpha = (.6), deltas=(.04), pi=(6000), EC50=(300), h=(15)),
                            V=5.019, low = c(.05, .001, 1000, 200, 12), 
                            high = c(3, .1, 10000, 355, 55)) 
  fit24_onlyADCC2 = fitmodel(patient24, decay24, parameters= fit24_onlyADCC$solution, 
                             V=5.019, low = c(.001, .001, 300, 50, .1), 
                             high = c(4, 1, 9000, 700, 100)) 
  fit24_onlyADCCfinal = finetune(patient24, decay24, parameters= fit24_onlyADCC2$solution, 
                                 V=5.019, low = c(.001, .001, 300, 50, .1), 
                                 high = c(4, 1, 9000, 700, 100)) 
  return(fit24_onlyADCCfinal)
}
completeADCCfit25 = function(patient25, decay25){
  fit25_onlyADCC = fitmodel(patient25, decay25, 
                            parameters = c(alpha = (.15), deltas=(.01), pi=(500), EC50=(200), h=(.1)),
                            V=27.090, low = c(.05, .001, 100, 20, .05), 
                            high = c(1, 2, 10000, 2000, 12)) 
  fit25_onlyADCC2 = fitmodel(patient25, decay25, parameters= fit25_onlyADCC$solution, 
                             V=27.090, low = c(.1, 10^-5, 100, 20, .05), 
                             high = c(3, 2, 10000, 2000, 12)) #this function call uses the previous parameters
  fit25_onlyADCCfinal = finetune(patient25, decay25, parameters= fit25_onlyADCC2$solution, 
                                 V=27.090, low = c(.1, 10^-5, 100, 20, .05), 
                                 high = c(3, 2, 10000, 2000, 12))
  return(fit25_onlyADCCfinal)
} #this one isnt working in the integrator
completeADCCfit26 = function(patient26, decay26){
  fit26_onlyADCC = fitmodel(patient26, decay26, 
                            parameters = c(alpha = (.15), deltas=(.01), pi=(230), EC50=(800), h=(8)),
                            V=5.141, low = c(.05, .001, 100, 100, 1), 
                            high = c(1, 1, 600, 1200, 12)) 
  fit26_onlyADCC2 = fitmodel(patient26, decay26, parameters= fit26_onlyADCC$solution, 
                             V=5.141, low = c(.05, .001, 100, 100, 1), 
                             high = c(1, 1, 600, 1200, 12)) 
  fit26_onlyADCCfinal = finetune(patient26, decay26, parameters= fit26_onlyADCC2$solution, 
                                 V=5.141, low = c(.05, .001, 100, 100, 1), 
                                 high = c(1, 1, 600, 1200, 12)) 
  return(fit26_onlyADCCfinal)
}
completeADCCfit27 = function(patient27, decay27){
  fit27_onlyADCC = fitmodel(patient27, decay27, 
                            parameters = c(alpha = (1), deltas=(.8), pi=(500), EC50=(350), h=(40)),
                            V=.237, low = c(.01, .001, 50, 100, 1), 
                            high = c(1.5, 1, 2000, 500, 60)) 
  fit27_onlyADCC2 = fitmodel(patient27, decay27, parameters= fit27_onlyADCC$solution, 
                             V=.237, low = c(.01, .001, 50, 100, .1), 
                             high = c(1.5, 1, 2000, 500, 20))
  fit27_onlyADCCfinal = finetune(patient27, decay27, parameters= fit27_onlyADCC2$solution, 
                                 V=.237, low = c(.01, .001, 50, 100, .1), 
                                 high = c(1.5, 1, 2000, 500, 20))
  return(fit27_onlyADCCfinal)
}

#this fits the model to the data with a combination of a global algorithm and a local algorithm
fitmodel = function(patient, decay, parameters, V, low, high) {
  beta=.001; deltai=1;
  known=c(V, deltai, beta)
  times = patient$studytime
  fit=nloptr(x0=parameters, eval_f=cost, lb=low, ub=high,
             patient=patient, decay=decay, times=times, known=known,
             opts=list(algorithm="NLOPT_GN_MLSL_LDS", maxeval=10000, maxtime=120,
              local_opts=list(algorithm="NLOPT_LN_COBYLA", xtol_rel=1e-10)))
  modtime = seq(times[1], 200, 1)#times[length(times)], 1)
  modparam = fit$solution
  modelopt = getmodel(modparam, patient, decay, modtime, known)
  #browser()
  plot(times, log10(patient$viral), ylim=c(1, 5), xlim=c(-20,200))
  lines(modtime, log10(1000*modelopt$V), type='l')
  conc=NULL
  for(i in 1:length(modtime)){
    if(modtime[i]<0){
      conc[i]=0
    }
    else{
      r1 = coef(decay)[1]
      r2 = coef(decay)[2]
      a = coef(decay)[3]
      b = coef(decay)[4]
      conc[i] = a*exp(r1*modtime[i])+a*b*exp(r2*modtime[i])
    }
  }
  #lines(modtime, log10(conc), col=2)
  return(fit)
}

#this fits the model to the data with only a local algorithm
finetune = function(patient, decay, parameters, V, low, high) {
  beta=.001; deltai=1;
  known=c(V, deltai, beta)
  times = patient$studytime
  fit=nloptr(x0=parameters, eval_f=cost, lb=low, ub=high,
             patient=patient, decay=decay, times=times, known=known,
             opts=list(algorithm="NLOPT_LN_COBYLA", xtol_rel=1e-90, xtol_abs=1e-90, maxeval=500))
  modtime = seq(times[1], 200, .5)
  modparam = fit$solution
  modelopt = getmodel(modparam, patient, decay, modtime, known)
  plot(times, log10(patient$viral), ylim=c(1, 5), xlim=c(-20,200))
  lines(modtime, log10(1000*modelopt$V), type='l')
  conc=NULL
  for(i in 1:length(modtime)){
    if(modtime[i]<0){
      conc[i]=0
    }
    else{
      r1 = coef(decay)[1]
      r2 = coef(decay)[2]
      a = coef(decay)[3]
      b = coef(decay)[4]
      conc[i] = a*exp(r1*modtime[i])+a*b*exp(r2*modtime[i])
    }
  }
  return(fit)
}

#function to be optimitzed
cost = function(parameters, patient, decay, times, known){
  model = getmodel(parameters, patient, decay, times, known)
  for(i in 1:length(times)){
    if(model$S[i]<1*10^-8 || is.na(model$S[i])){
      model$S[i]=10^-8
    }
    if(model$I[i]<1*10^-8 || is.na(model$I[i])){
      model$I[i]=10^-8
    }
    if(model$V[i]<1*10^-8 || is.na(model$V[i])){
      model$V[i]=10^-8
    }
  }
  vm = log10(1000*model$V)
  vd = log10(patient$viral)
  #browser()
  diff = 0
  for (i in 1:length(times)){
    diff = diff + abs(vm[i] - vd[i])
  }
  return(diff)
}

#gets data given any set of parameters
getmodel = function(parameters, patient, decay, times, known){
  alpha = unname(parameters[1])
  deltas = unname(parameters[2])
  pi = unname(parameters[3])
  EC50 = unname(parameters[4])
  h = unname(parameters[5])
  V = unname(known[1])
  deltai = unname(known[2])
  beta = unname(known[3])
  S=(alpha)/((deltas)+(beta)*unname(V))
  I=(beta)*unname(S)*unname(V)/(deltai)
  c=pi*I/V
  #browser()
  state = c(S=S, I=I, V=V)
  model = as.data.frame(ode(y=state, times=times, func=ODE, parms=c(alpha=alpha,deltas=deltas,beta=beta,
                                                                    deltai=deltai,pi=pi,c=c,EC50=EC50,h=h), 
                            decay=decay))
}

#diff equations
ODE = function(times, state, parms, decay) {
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
    IAP = log10(1+(conc/EC50)^h)
    fa = 1-1/(10^IAP)
    if(V<10^(-1.7) || is.na(V)){
      V=10^(-1.7)
    }
    dS = alpha - deltas*S -beta*S*V
    dI = beta*S*V - (deltai+fa)*I 
    dV = pi*I - c*V
    if(V<10^(-1.7) || is.na(V)){
      dV=0
    }
    list(c(dS, dI, dV))
  })
}

#this function graphs viral load over time given the fit and patient information
graphing = function(fit, patient, decay, known){
  times = seq(-20, 60, .1)
  param = fit$solution
  alpha = param[1]
  deltas = param[2]
  pi = param[3]
  EC50 = param[4]
  h = param[5]
  V = known[1]
  deltai = known[2]
  beta = known[3]
  S = unname(alpha)/(unname(deltas)+unname(beta)*unname(V))
  I = unname(beta)*unname(S)*unname(V)/unname(deltai)
  c = unname(pi)*unname(I)/unname(V)
  model = as.data.frame(ode(y=c(S=S, I=I, V=V ), times=times, func=ODE, parms=c(alpha=alpha,deltas=deltas,beta=beta,
                                                                                deltai=deltai,pi=pi,c=c,EC50=EC50,h=h), 
                            decay=decay, maxsteps=10e8))
  plot(times, log10(1000*model$V), type='l', ylim=c(1,5), xlim=c(-20,60))
  par(new=T)
  plot(patient$studytime, log10(patient$viral), ylim=c(1,5), xlim=c(-20,60))
}


