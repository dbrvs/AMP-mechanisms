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

completebothfit20 = function(patient20, decay20){
  fit20_both = fitmodel(patient20, decay20, 
                        parameters = c(alpha=1, deltas=.2, pi=900, IC50=130, EC50=60, m=20, h=20),
                        V=3.547, low=c(.2, .01, 500, 30, 30, 1, 1), 
                        high=c(2, 1, 10000, 600, 600, 60, 60)) 
  fit20_both2 = fitmodel(patient20, decay20, parameters=fit20_both$solution,
                         V=3.547, low=c(.2, .01, 500, 30, 30, 0, 0), 
                         high=c(3, 1, 1700, 10000, 10000, 10000, 10000)) 
  fit20_bothfinal = finetune(patient20, decay20, parameters=fit20_both2$solution,
                             V=3.547, low=c(.5, .01, 500, 30, 30, 15, 15), 
                             high=c(3, 1, 1700, 400, 400, 35, 35))
  return(fit20_bothfinal)
}
completebothfit21 = function(patient21, decay21){
  fit21_both = fitmodel(patient21, decay21, 
                        parameters = c(alpha=.2, deltas=.8, pi=1000, IC50=400, EC50=400, m=10, h=10),
                        V=6.551, low=c(.1, .3, 900, 100, 100, 5, 10), 
                        high=c(1, 1, 1500, 1000, 1500, 50, 50)) 
  fit21_both2 = fitmodel(patient21, decay21, parameters=fit21_both$solution,
                         V=6.551, low=c(.1, .3, 900, 100, 100, 5, 10), 
                         high=c(1, 1, 1500, 1000, 1500, 50, 50)) 
  fit21_bothfinal = finetune(patient21, decay21, parameters=fit21_both2$solution,
                             V=6.551, low=c(.1, .3, 900, 100, 100, 5, 10), 
                             high=c(1, 1, 1500, 1000, 1500, 50, 50)) 
  return(fit21_bothfinal)
}
completebothfit22 = function(patient22, decay22){
  fit22_both = fitmodel(patient22, decay22, 
                        parameters = c(alpha=.15, deltas=.01, pi=700, IC50=75, EC50=150, m=1, h=1),
                        V=.745, low=c(.01, .01, 500, 50, 50, 1, 1), 
                        high=c(1, .5, 4000, 200, 200, 25, 32)) 
  fit22_both2 = fitmodel(patient22, decay22, parameters=fit22_both$solution,
                         V=.745, low=c(.01, .01, 500, 50, 50, 1, 1), 
                         high=c(1, .5, 4000, 200, 200, 25, 32)) 
  fit22_bothfinal = finetune(patient22, decay22, parameters=fit22_both2$solution,
                             V=.745, low=c(.01, .01, 500, 50, 50, 1, 1), 
                             high=c(1, .5, 4000, 200, 200, 25, 32)) 
  return(fit22_bothfinal)
}
completebothfit23 = function(patient23, decay23){
  fit23_both = fitmodel(patient23, decay23, 
                        parameters = c(alpha=.15, deltas=.4, pi=1000, IC50=100, EC50=100, m=20, h=20),
                        V=27.894, low=c(.05, .001, 500, 50, 50, 10, 10), 
                        high=c(1.5, 1, 3000, 400, 400, 60, 60)) 
  fit23_both2 = fitmodel(patient23, decay23, parameters=fit23_both$solution,
                         V=27.894, low=c(.05, .001, 500, 50, 50, 10, 10), 
                         high=c(1.5, 1, 3000, 400, 400, 60, 60)) 
  fit23_bothfinal = finetune(patient23, decay23, parameters=fit23_both2$solution,
                             V=27.894, low=c(.05, .001, 500, 50, 50, 10, 10), 
                             high=c(1.5, 1, 3000, 400, 400, 60, 60)) 
  return(fit23_bothfinal)
}
completebothfit24 = function(patient24, decay24){
  fit24_both = fitmodel(patient24, decay24, 
                        parameters = c(alpha=.7, deltas=.009, pi=2000, IC50=400, EC50=500, m=55, h=55),
                        V=5.019, low=c(.1, .005, 800, 100, 100, 1, 1), 
                        high=c(3, .1, 9000, 600, 600, 100, 100))  
  fit24_both2 = fitmodel(patient24, decay24, parameters=fit24_both$solution,
                         V=5.019, low=c(.05, .001, 500, 50, 50, 15, 15), 
                         high=c(1, 1, 1700, 100, 100, 50, 50)) 
  fit24_bothfinal = finetune(patient24, decay24, parameters=fit24_both2$solution,
                             V=5.019, low=c(.05, .001, 500, 50, 50, 15, 15), 
                             high=c(1, 1, 1700, 100, 100, 50, 50)) 
  return(fit24_bothfinal)
}
completebothfit25 = function(patient25, decay25){
  fit25_both = fitmodel(patient25, decay25, 
                        parameters = c(alpha=.15, deltas=.01, pi=700, IC50=75, EC50=654, m=20, h=30),
                        V=27.090, low=c(.05, .001, 500, 50, 50, 5, 5), 
                        high=c(1, .5, 1700, 800, 800, 25, 32)) 
  fit25_both2 = fitmodel(patient25, decay25, parameters=fit25_both$solution,
                         V=27.090, low=c(.1, .01, 600, 200, 500, .5, .5), 
                         high=c(1, .5, 1000, 800, 600, 50, 50))
  fit25_bothfinal = finetune(patient25, decay25, parameters=fit25_both2$solution,
                             V=27.090, low=c(.05, .001, 500, 50, 50, 5, 5), 
                             high=c(1, .5, 1700, 800, 850, 40, 50))
  return(fit25_bothfinal)
} #DOESNT WORK IN INTEGRATOR
completebothfit26 = function(patient26, decay26){
  fit26_both = fitmodel(patient26, decay26, 
                        parameters = c(alpha=.5, deltas=.01, pi=700, IC50=350, EC50=350, m=10, h=10),
                        V=5.141, low=c(.1, .001, 500, 300, 300, 5, 5), 
                        high=c(1, 1, 1700, 700, 700, 30, 30)) 
  fit26_both2 = fitmodel(patient26, decay26, parameters=fit26_both$solution,
                         V=5.141, low=c(.1, .001, 500, 300, 300, 5, 5), 
                         high=c(1, 1, 1700, 700, 700, 30, 30)) 
  fit26_bothfinal = finetune(patient26, decay26, parameters=fit26_both2$solution,
                             V=5.141, low=c(.1, .001, 500, 300, 300, 5, 5), 
                             high=c(1, 1, 1700, 700, 700, 30, 30)) 
  return(fit26_bothfinal)
}
completebothfit27 = function(patient27, decay27){
  fit27_both = fitmodel(patient27, decay27, 
                        parameters = c(alpha=1, deltas=.5, pi=700, IC50=400, EC50=700, m=20, h=20),
                        V=.237, low=c(1, .1, 500, 50, 50, 20, 20), 
                        high=c(3, 1, 10000, 800, 800, 50, 50)) 
  fit27_both2 = fitmodel(patient27, decay27, parameters=fit27_both$solution,
                         V=.237, low=c(1, .1, 500, 50, 50, 20, 20), 
                         high=c(3, 1, 10000, 800, 800, 50, 50)) 
  fit27_bothfinal = finetune(patient27, decay27, parameters=fit27_both2$solution,
                             V=.237, low=c(1, .1, 500, 50, 50, 20, 20), 
                             high=c(3, 1, 10000, 800, 800, 50, 50)) 
  return(fit27_bothfinal)
}

#this fits the model to the data with a combination of a global algorithm and a local algorithm
fitmodel = function(patient, decay, parameters, V, low, high) {
  beta=.001; deltai=1;
  known=c(V, deltai, beta) #these are the fixed values
  times = patient$studytime
  fit=nloptr(x0=parameters, eval_f=cost, lb=low, ub=high,
             patient=patient, decay=decay, times=times, known=known,
             opts=list(algorithm="NLOPT_GN_MLSL_LDS", maxeval=10000,# maxtime=120,
                       local_opts=list(algorithm="NLOPT_LN_COBYLA", xtol_rel=1e-10)))
  modtime = seq(times[1], 200, 1)
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
  known=c(V, deltai, beta) #known values
  times = patient$studytime
  fit=nloptr(x0=parameters, eval_f=cost, lb=low, ub=high,
             patient=patient, decay=decay, times=times, known=known,
             opts=list(algorithm="NLOPT_LN_COBYLA", xtol_rel=1e-50, xtol_abs=1e-50, maxeval=500))
  modtime = seq(times[1], 200, 1)
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

#function to be optimized
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
  #browser()
  state = c(S=S, I=I, V=V)
  model = as.data.frame(ode(y=state, times=times, func=ODE, 
                            parms=c(alpha=alpha,deltas=deltas,beta=beta,deltai=deltai,
                                    pi=pi,c=c,IC50=IC50,EC50=EC50,m=m,h=h), 
                            decay=decay, maxsteps=1e6))
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
    IIP = log10(1+(conc/IC50)^m)
    fu = 1/(10^IIP)
    IAP = log10(1+(conc/EC50)^h)
    fa = 1-1/(10^IAP)
    if(V<10^(-1.7) || is.na(V)){
      V=10^(-1.7)
    }
    dS = alpha - deltas*S -fu*beta*S*V
    dI = fu*beta*S*V - (deltai+fa)*I 
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
  EC50 = param[5]
  m = param[6]
  h = param[7]
  V = known[1]
  deltai = known[2]
  beta = known[3]
  S = unname(alpha)/(unname(deltas)+unname(beta)*unname(V))
  I = unname(beta)*unname(S)*unname(V)/unname(deltai)
  c = unname(pi)*unname(I)/unname(V)
  model = as.data.frame(ode(y=c(S=S, I=I, V=V ), times=times, func=ODE, parms=c(alpha=alpha,deltas=deltas,beta=beta,
                                                                                deltai=deltai,pi=pi,c=c,IC50=IC50,EC50=EC50,m=m,h=h), 
                            decay=decay, maxsteps=10e8))
  plot(times, log10(1000*model$V), type='l', ylim=c(1,5), xlim=c(-20,60))
  par(new=T)
  plot(patient$studytime, log10(patient$viral), ylim=c(1,5), xlim=c(-20,60))
}






