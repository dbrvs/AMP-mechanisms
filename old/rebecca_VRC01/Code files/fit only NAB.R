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

completeNABfit20 = function(patient20, decay20){
  fit20_onlyNAB = fitmodel(patient20, decay20, 
                           parameters = c(alpha=(.6), deltas=(.02), pi=(3000), IC50=(150), m=(20)),
                           V=3.547, low = c(.5, .001, 2200, 100, 3), 
                           high = c(1.5, .03, 7000, 270, 100)) 
  fit20_onlyNAB2 = fitmodel(patient20, decay20, parameters= fit20_onlyNAB$solution, 
                            V=3.547, low = c(.01, .01, 50, 50, .1), 
                            high = c(5, 1, 50000, 900, 30)) 
  fit20_onlyNABfinal = finetune(patient20, decay20, parameters= fit20_onlyNAB$solution, 
                                V=3.547, low = c(.01, .01, 50, 50, .1), 
                                high = c(5, 1, 50000, 900, 30)) 
  return(fit20_onlyNABfinal)
}
completeNABfit21 = function(patient21, decay21) {
  fit21_onlyNAB = fitmodel(patient21, decay21, 
                           parameters = c(alpha=(.5), deltas=(.01), pi=(1000), IC50=(1000), m=(7)),
                           V=6.551, low = c(.1, .001, 1000, 900, 4), 
                           high = c(2, .1, 8000, 1700, 200))
  fit21_onlyNAB2 = fitmodel(patient21, decay21, parameters= fit21_onlyNAB$solution, 
                            V=6.551, low = c(.3, .01, 500, 100, 4), 
                            high = c(.9, 1, 8000, 1000, 15))
  fit21_onlyNABfinal = finetune(patient21, decay21, parameters= fit21_onlyNAB2$solution, 
                                V=6.551, low = c(.3, .01, 500, 100, 4), 
                                high = c(.9, 1, 8000, 1000, 15))
  return(fit21_onlyNABfinal)
}
completeNABfit22 = function(patient22, decay22) {
  fit22_onlyNAB = fitmodel(patient22, decay22, 
                           parameters = c(alpha=(.39), deltas=(.009), pi=(1000), IC50=(100), m=(10)),
                           V=.745, low = c(.1, .001, 300, 10, .5), 
                           high = c(1.5, 1, 5000, 500, 30))
  fit22_onlyNAB2 = fitmodel(patient22, decay22, parameters= fit22_onlyNAB$solution, 
                            V=.745, low = c(.1, .001, 300, 10, .5), 
                            high = c(1.5, 1, 5000, 500, 30))
  fit22_onlyNABfinal = finetune(patient22, decay22, parameters= fit22_onlyNAB2$solution, 
                                V=.745, low = c(.1, .001, 300, 10, .5), 
                                high = c(1.5, 1, 5000, 500, 30))
  return(fit22_onlyNABfinal)
}
completeNABfit23 = function(patient23, decay23) {
  fit23_onlyNAB = fitmodel(patient23, decay23, 
                           parameters = c(alpha=(.6), deltas=(.009), pi=(3000), IC50=(160), m=(5)),
                           V=27.894, low = c(.5, .001, 2200, 160, 3), 
                           high = c(1.5, .03, 7000, 270, 100)) 
  fit23_onlyNAB2 = fitmodel(patient23, decay23, parameters= fit23_onlyNAB$solution, 
                            V=27.894, low = c(.5, .001, 2200, 160, 3), 
                            high = c(1.5, .03, 7000, 270, 100)) 
  fit23_onlyNABfinal = finetune(patient23, decay23, parameters= fit23_onlyNAB2$solution, 
                                V=27.894, low = c(.5, .001, 2200, 160, 3), 
                                high = c(1.5, .03, 7000, 270, 100)) 
  return(fit23_onlyNABfinal)
}
completeNABfit24 = function(patient24, decay24) {
  fit24_onlyNAB = fitmodel(patient24, decay24, 
                           parameters = c(alpha=(.6), deltas=(.009), pi=(5500), IC50=(150), m=(3)),
                           V=5.019, low = c(.3, .003, 4800, 50, 3), 
                           high = c(.9, .05, 10000, 400, 7))
  fit24_onlyNAB2 = fitmodel(patient24, decay24, parameters= fit24_onlyNAB$solution, 
                            V=5.019, low = c(.3, .003, 2000, 50, 3), 
                            high = c(.9, .05, 10000, 400, 7))
  fit24_onlyNABfinal = finetune(patient24, decay24, parameters= fit24_onlyNAB2$solution, 
                                V=5.019, low = c(.3, .003, 2000, 50, 3), 
                                high = c(.9, .05, 10000, 400, 7))
  return(fit24_onlyNABfinal)
}
completeNABfit25 = function(patient25, decay25) {
  fit25_onlyNAB = fitmodel(patient25, decay25, 
                           parameters = c(alpha=(.39), deltas=(.009), pi=(1000), IC50=(20000), m=(2.2)),
                           V=27.090, low = c(.1, .003, 1000, 50, .5), 
                           high = c(.9, .02, 10000, 2000000, 7))
  fit25_onlyNAB2 = fitmodel(patient25, decay25, parameters= fit25_onlyNAB$solution, 
                            V=27.090, low = c(.1, .003, 1000, 50, .5), 
                            high = c(.9, .02, 10000, 2000000, 7)) 
  fit25_onlyNABfinal = finetune(patient25, decay25, parameters= fit25_onlyNAB2$solution, 
                                V=27.090, low = c(.1, .003, 1000, 50, .5), 
                                high = c(.9, .02, 10000, 2000000, 7))
  return(fit25_onlyNABfinal)
} #NOT WORKING
completeNABfit26 = function(patient26, decay26) {
  fit26_onlyNAB = fitmodel(patient26, decay26, 
                           parameters = c(alpha=(.6), deltas=(.6), pi=(1000), IC50=(400), m=(2.2)),
                           V=5.141, low = c(.1, .01, 500, 50, .4), 
                           high = c(.9, 1, 10000, 1000, 7))
  fit26_onlyNAB2 = fitmodel(patient26, decay26, parameters= fit26_onlyNAB$solution, 
                            V=5.141, low = c(.1, .01, 500, 50, .4), 
                            high = c(.9, 1, 10000, 1000, 7))
  fit26_onlyNABfinal = finetune(patient26, decay26, parameters= fit26_onlyNAB2$solution, 
                                V=5.141, low = c(.1, .01, 500, 50, .4), 
                                high = c(.9, 1, 10000, 1000, 7))
  return(fit26_onlyNABfinal)
}
completeNABfit27 = function(patient27, decay27) {
  fit27_onlyNAB = fitmodel(patient27, decay27, 
                           parameters = c(alpha=(.39), deltas=(.5), pi=(1000), IC50=(15), m=(10)),
                           V=5.019, low = c(.1, .01, 500, 5, .3), 
                           high = c(4, 1, 6000, 200, 16))
  fit27_onlyNAB2 = fitmodel(patient27, decay27, parameters= fit27_onlyNAB$solution, 
                            V=5.019, low = c(.1, .01, 500, 5, .3), 
                            high = c(4, 1, 6000, 200, 16))
  fit27_onlyNABfinal = finetune(patient27, decay27, parameters= fit27_onlyNAB2$solution, 
                                V=5.019, low = c(.1, .01, 500, 5, .3), 
                                high = c(4, 1, 6000, 200, 16))
  return(fit27_onlyNABfinal)
}

#this fits the model to the data with a combination of a global algorithm and a local algorithm
fitmodel = function(patient, decay, parameters, V, low, high) {
  beta = .001; deltai=1;
  known=c(V, deltai, beta)
  times = patient$studytime
  fit=nloptr(x0=parameters, eval_f=cost, lb=low, ub=high,
             patient=patient, decay=decay, times=times, known=known,
             opts=list(algorithm="NLOPT_GN_MLSL_LDS", maxtime=120, maxeval=10000,
                      local_opts=list(algorithm="NLOPT_LN_COBYLA", xtol_rel=1e-5)))
  modtime = seq(times[1], 200, .5) #times[length(times)], 1)
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
             opts=list(algorithm="NLOPT_LN_COBYLA", xtol_rel=1e-30, xtol_abs=1e-30, maxeval=500))
  modtime = seq(times[1], 200, .5)
  modparam = fit$solution
  modelopt = getmodel(modparam, patient, decay, modtime, known)
  plot(times, log10(patient$viral), ylim=c(1, 5), xlim=c(-20,200))
  lines(modtime, log10(1000*modelopt$V), type='l')
  #browser()
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
  for(i in 1:length(model)){
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
  for (i in 1:length(model)){
    diff = diff + abs(vm[i] - vd[i])
  }
  return(diff)
}

#gets data given any set of parameters
getmodel = function(parameters, patient, decay, times, known){
  alpha = unname(parameters[1])
  deltas = unname(parameters[2])
  pi = unname(parameters[3])
  IC50 = unname(parameters[4])
  m = unname(parameters[5])
  V = unname(known[1])
  deltai = unname(known[2])
  beta = unname(known[3])
  S = (alpha)/((deltas)+(beta)*unname(V))
  I = (beta)*unname(S)*unname(V)/(deltai)
  c = pi*I/V
  #browser()
  state = c(S=S, I=I, V=V)
  model = as.data.frame(ode(y=state, times=times, func=ODE, parms=c(alpha=alpha,deltas=deltas,beta=beta,
                                                                    deltai=deltai,pi=pi,c=c,IC50=IC50,m=m), 
                            decay=decay, maxsteps=10e8))
}

#diff equations
ODE = function(times, state, parms, decay) {
  with(as.list(c(state, parms)), {
    conc = NULL
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
    if(V<10^(-1.7) || is.na(V)){
      V=10^(-1.7)
    }
    IIP = log10(1+(conc/IC50)^m)
    fu = 1/(10^IIP)
    dS = alpha - deltas*S - fu*beta*S*V
    dI = fu*beta*S*V - (deltai)*I 
    dV = pi*I - c*V
    if(V<10^(-1.7) || is.na(V)){
      dV=0
    }
    #browser()
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
  IC50 = param[4]
  m = param[5]
  V = known[1]
  deltai = known[2]
  beta = known[3]
  S = unname(alpha)/(unname(deltas)+unname(beta)*unname(V))
  I = unname(beta)*unname(S)*unname(V)/unname(deltai)
  c = unname(pi)*unname(I)/unname(V)
  model = as.data.frame(ode(y=c(S=S, I=I, V=V ), times=times, func=ODE, parms=c(alpha=alpha,deltas=deltas,beta=beta,
                                                                    deltai=deltai,pi=pi,c=c,IC50=IC50,m=m), 
                            decay=decay, maxsteps=10e8))
  plot(times, log10(1000*model$V), type='l', ylim=c(1,5), xlim=c(-20,60))
  par(new=T)
  plot(patient$studytime, log10(patient$viral), ylim=c(1,5), xlim=c(-20,60))
}

