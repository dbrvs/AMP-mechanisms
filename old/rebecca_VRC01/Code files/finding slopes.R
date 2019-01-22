#this file fits exponential decays to the patients' viral load data
#in order to determine their rate of viral decay

#data from Lynch study

library(ggplot2)

#the "make" functions create the data sets for each patient
make20 = function(){
  studytime = c(-19.94236311239193, 0, 0.17291066282420786, 0.9798270893371708,
                2.0172910662824144, 3.0547550432276616, 5.014409221902017, 6.974063400576366, 
                8.933717579250715, 14.005763688760805, 15.96541786743516, 20.979827089337174, 
                34.92795389048992, 41.90201729106627, 48.933717579250725, 56, 
                84)
  viral = c(10^3.65891472868217, 10^3.410852713178294, 10^2.9534883720930227,  10^3.3720930232558137,
            10^3.3100775193798446, 10^3.1162790697674416, 10^2.837209302325581, 10^2.581395348837209,
            10^2.410852713178294, 10^2.627906976744186, 10^2.8449612403100772, 10^2.976744186046511,
            10^3.527131782945736, 10^3.457364341085271, 10^3.705426356589147, 10^2.8600550964187312,
            10^3.619283746556472)
  patient20 = data.frame(studytime, viral)
  plot(patient20$studytime, log10(patient20$viral), type='o', ylim=c(1,5))
  return(patient20)
}
make21 = function(){
  studytime = c(-25, -10.024096385542165, -0.048192771084334396, 1.0120481927710898, 2.0240963855421725, 
                2.9879518072289173, 5.01204819277109, 6.987951807228917, 9.01204819277109, 
                12, 13.975903614457835, 16, 20.963855421686745, 
                28, 34.98795180722892, 42.02409638554216,49.01204819277108, 
                56, 84)
  viral = c(10^3.8416666666666655, 10^3.8156000625880147, 10^3.8188389923329686, 10^3.871131278360194, 10^3.845485839461743,
            10^3.748396182131122, 10^3.56723517446409, 10^3.5808637145986544, 10^3.601001408230324,
            10^3.7643091847911125, 10^3.8104052573932092, 10^3.7526208731028006, 10^3.9944922547332182,
            10^3.7370364575183848, 10^4.031513065248005, 10^3.9363949303708337, 10^3.8542481614770767,
            10^3.719559228650136, 10^3.958126721763083)
  patient21 = data.frame(studytime, viral)
  plot(patient21$studytime, log10(patient21$viral), type='o', ylim=c(1,5))
  return(patient21)
}
make22 = function(){
  studytime = c(-28, -18.65036231884058, 0.25724637681159024, 0.365036231884055, 1.192934782608699, 
                2.2346014492753703, 3.1847826086956594, 5.133152173913043, 7.099637681159418, 
                9.020833333333332, 12.03532608695652, 14.006340579710145, 15.97735507246377, 
                20.962862318840582, 27.977355072463762, 34.93387681159422, 48.97463768115942,
                56, 84)
  viral = c(10^2.55, 10^3.1328125000000004, 10^2.7109375000000004, 10^2.640625, 10^2.7812500000000004,
            10^2.7656250000000004, 10^2.4609375000000004, 10^2.2656250000000004, 10^1.7265625,
            10^1.2968750000000009, 10^1.2968750000000009, 10^1.2968750000000009, 10^1.2968750000000009,
            10^1.2968750000000009, 10^1.2968750000000009, 10^1.2968750000000009, 10^1.8984375000000009,
            10^3.215426997245177, 10^3.015977961432505)
  patient22 = data.frame(studytime, viral)
  plot(patient22$studytime, log10(patient22$viral), type='o', ylim=c(1,5))
  return(patient22)
}
make23 = function(){
  studytime = c(-11.018867924528301, 0.15094339622641328, 1.0566037735849, 1.9999999999999964, 
                3.0188679245282977, 4.981132075471695, 6.981132075471695, 8.943396226415093, 
                11.999999999999996, 13.962264150943394, 15.999999999999996, 20.981132075471695, 
                27.924528301886784, 34.9433962264151, 41.92452830188678, 48.90566037735849,
                56, 84)
  viral = c(10^4.5786802030456855, 10^4.269035532994924, 10^4.2639593908629445, 10^4.203045685279188,
            10^4.00507614213198, 10^3.2182741116751274, 10^3.106598984771574, 10^3.106598984771574,
            10^3.6142131979695438, 10^3.8020304568527923, 10^3.9644670050761426, 10^3.9289340101522847,
            10^4.781725888324873, 10^4.8883248730964475, 10^4.4822335025380715, 10^4.324873096446701,
            10^4.256749311294764, 10^4.858953168044075)
  patient23 = data.frame(studytime, viral)
  plot(patient23$studytime, log10(patient23$viral), type='o', ylim=c(1,5))
  return(patient23)
}
make24 = function(){
  studytime = c(-12.798262071689539, 0.24725980053323227, 0.37404957045522025, 1.1975905993877767,
                2.2391626345413265, 3.2151673743458105, 5.135973141107932, 7.129455909943708, 
                9.131233336624867, 12.131924558111983, 14.14160165893157, 16.187617260787988, 
                21.171916658437837, 28.210131332082554, 35.186728547447416, 42.19571442677989,
                49.216944801026955, 56, 84)
  viral = c(10^3.4822158585958323, 10^3.726276291102992, 10^3.3808827885849704,  10^3.6280932161548334,
              10^3.7206082749086598, 10^3.457470129357164, 10^2.523985385602844, 10^2.5389355189098453,
              10^2.6621309370988446, 10^2.8211513775056787, 10^3.047437543201343, 10^3.2479411474276687,
              10^3.293048286758172, 10^3.641749777821665, 10^3.6863434383331697, 10^3.653609163622001,
              10^3.7806655475461644, 10^3.7443526170798878, 10^3.338292011019282)
  patient24 = data.frame(studytime, viral)
  plot(patient24$studytime, log10(patient24$viral), type='o', ylim=c(1,5))
  return(patient24)
}
make25 = function(){
  studytime = c(0.3016236045379337, 0.14265981295712749, 1.1833971151015525, 2.1629944974072117, 
                3.131269672335307, 5.082312447634791, 7.070944951427727, 9.028328162858621, 
                12.075813500600077, 14.127850365707296, 16.128257965172885, 21.166187359887687,
                28.197278141346438,35.13507393401417, 42.16073005593171, 49)
  viral = c(10^4.178800298906274, 10^4.4366748941373, 10^4.38939335612871, 10^4.287448200901247,
            10^3.9901949684110405, 10^2.645682842326939, 10^2.9495935327551446, 10^2.7144539299381796,
            10^3.2835760059781247, 10^3.681211929077692, 10^4.188243019859151, 10^4.1916666666666655, 
            10^4.099999999999999, 10^4.383333333333332, 10^4.0666666666666655, 10^4.249999999999999)
  patient25 = data.frame(studytime, viral)
  plot(patient25$studytime, log10(patient25$viral), type='o', ylim=c(1,5))
  return(patient25)
}
make26 = function(){
  studytime = c(-20.068992706485318, -0.024837374334712337, 0.89256849990144, 4.874039030159668,
                6.9016361127537955, 11.909323871476449, 13.8963138182535, 20.820224719101123,
                27.887246205401144, 34.85314409619555, 41.84525921545438, 48.84249950719497,
                56, 84)
  viral = c(10^3.7876700177409828, 10^3.6330080820027613, 10^3.6428444707273817, 10^3.527390104474671,
            10^3.505726394638282, 10^3.554691504041003, 10^3.5949043958210143, 10^3.534634338655629,
            10^3.732022471910114, 10^3.5737926276365086, 10^3.729982259018334, 10^3.752148630001973,
            10^4.240220385674929, 10^3.759779614325067) 
  patient26 = data.frame(studytime, viral)
  plot(patient26$studytime, log10(patient26$viral), type='o', ylim=c(1,5))
  return(patient26)
}
make27 = function(){
  studytime = c(-63.125, 0.058997050147496566, 0.9439528023598847, 2.005899705014752, 3.0088495575221295,
                4.955752212389381, 6.961651917404129, 8.967551622418878, 11.976401179941004,
                13.982300884955752, 15.988200589970507, 21.002949852507374, 28.023598820059004,
                35.044247787610615, 42.005899705014755,49.02654867256638, 56,
                84)
  viral = c(10^2.4833333333333334, 10^2.2449538940375824, 10^2.442152695514831, 10^2.57642905256312, 10^1.962650686363319,
            10^1.6327085220542124, 10^1.3185376164262657, 10^1.3193273407195782, 10^1.320511927159548,
            10^1.3213016514528606, 10^1.3220913757461732, 10^1.316191670731425, 10^1.3268297215060514,
            10^1.5185701344854028, 10^2.2693424383898924, 10^2.791791512786567, 10^2.438567493112946,
            10^3.015977961432505)
  patient27 = data.frame(studytime, viral)
  plot(patient27$studytime, log10(patient27$viral), type='o', ylim=c(1,5))
  return(patient27)
}
plotall = function(patient20, patient21, patient22, patient23,
                   patient24, patient25, patient26, patient27){
  ggplot() +
    scale_x_continuous(limits = c(-65,90)) +
    scale_y_continuous(limits = c(1,5)) +
    geom_point(data = patient20, aes(studytime, log10(viral)), color="blue") + 
    geom_line(data = patient20, aes(studytime, log10(viral)), color="blue") +
    geom_point(data = patient21, aes(studytime, log10(viral)), color="red") + 
    geom_line(data = patient21, aes(studytime, log10(viral)), color="red") +
    geom_point(data = patient22, aes(studytime, log10(viral)), color="green") + 
    geom_line(data = patient22, aes(studytime, log10(viral)), color="green") +
    geom_point(data = patient23, aes(studytime, log10(viral)), color="purple") + 
    geom_line(data = patient23, aes(studytime, log10(viral)), color="purple") +
    geom_point(data = patient24, aes(studytime, log10(viral)), color="orange") + 
    geom_line(data = patient24, aes(studytime, log10(viral)), color="orange") +
    geom_point(data = patient25, aes(studytime, log10(viral)), color="brown") + 
    geom_line(data = patient25, aes(studytime, log10(viral)), color="brown") +
    geom_point(data = patient26, aes(studytime, log10(viral)), color="pink") + 
    geom_line(data = patient26, aes(studytime, log10(viral)), color="pink") +
    geom_point(data = patient27, aes(studytime, log10(viral)), color="yellow") + 
    geom_line(data = patient27, aes(studytime, log10(viral)), color="yellow") +
    labs(x="days after VRC01 infusion of 40mg/kg", y="log 10 viral load") +
    labs(title = "Patient viral loads") +
    theme(plot.title = element_text(size=27))
    
}
#make the patients
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

#this function takes in a patient, the index of the first data point to fit to, 
#the index of the last data point to fit to, n (the patient number), an initial guess for a (close to V0)
#and the x and y coordinates for where to print the rate value on the graph
getslope = function(patient, start, end, n, aguess, xtext, ytext){
  x = patient$studytime[start:end]
  pred = seq(x[1], x[length(x)], .01) #time series for plotting prediction
  y = patient$viral[start:end]
  Vinit = patient$viral[start]
  Tinit = patient$studytime[start]
  fit = nls(y ~ a*(exp(r*(x-Tinit))), start=c(r=.5,a=aguess))
  plot(x, y, xlab='days after VRC01 infusion', ylab='viral load')
  r = coef(fit)[1]
  a = coef(fit)[2]
  lines(pred, a*exp(r*(pred-Tinit)), col=2)
  title(paste("#",n))
  text(paste("rate = ",round(r, digits=2)), x=xtext,y=ytext, cex=1)
  return(fit)
}

#Fitting and graphing results
par(mfrow=c(2,3))
fit20 = getslope(patient20, 4, 9, 20, 1000, 6, 1750)
#fit21 = getslope(patient21, 5, 7, 21, 1000, 4, 6000) #unresponsive
fit22 = getslope(patient22, 5, 10, 22, 1000, 7, 500)
fit23 = getslope(patient23, 4, 6, 23, 14000, 4, 11000)
fit24 = getslope(patient24, 5, 7, 24, 1000, 4.2, 4000)
fit25 = getslope(patient25, 3, 5, 25, 16000, 2.7, 17000)
#fit26 = getslope(patient26, 3, 5, 26, 1000, 5, 4000) #unresponsive
fit27 = getslope(patient27, 4, 7, 27, 1000, 5, 300)
rs = c(coef(fit20)[1], coef(fit22)[1], coef(fit23)[1],
       coef(fit24)[1], coef(fit25)[1], coef(fit27)[1])
mean(rs) #not including patients #21 and #26 (unresponsive)


