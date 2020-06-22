#library(TSDpredict)
#demo(simulatePredictionExample)
#Simulating examples to predict age
#rm(list=ls())
ntrain = 18
ntest = 1 #repeat test
time = 3
TrueTimes = rep(time,ntest)

#paste0(round(beta,1),collapse=",")
beta = c(5,4.9,5.7,6,4.3,5.6,5.9,5.1,-0.9,-0.8,-1.2,-1,-0.6,-0.9,-0.9,-0.9)
beta = t(matrix(beta,ncol=2))
nsites = ncol(beta)
SIGMA = c(1.1,0.2,0.4,0.8,0.5,0.2,0.1,0.5,0.2,0.8,0.4,0.5,0.5,0.4,0.3,0.4,0.4,0.4,0.6,0.8,0.2,0.2,0.2,0.5,0.8,0.5,0.8,1.7,0.4,0.1,0.4,0.7,0.5,0.5,0.2,0.4,1.2,0.7,0.4,0.2,0.2,0.4,0.2,0.1,0.7,0.9,0.3,0.2,0.1,0.3,0.2,0.4,0.4,0.3,0.5,0.2,0.5,0.4,0.5,0.7,0.2,0.2,0.2,0.5)
SIGMA = matrix(SIGMA,ncol=nsites)

dat = genData(ntrain,nsites,beta,SIGMA,maxTime = 5,seed=1)

#Show train data:
times = dat$Data_time
sites = colnames(dat$Data_Y)
  
if(0) {
  par(mfrow=c(nsites/2,2))
  for(j in 1:nsites) {
    plot( times, dat$Data_Y[,j],main=sites[j],xlab="Time",ylab="Response (log (Rfu+1)")
  }
}

test = genData(ntest,nsites,beta,SIGMA,times = TrueTimes,seed=5)
New_Y = test$Data_Y #responses to test with

pred_NaiveMLE = posteriorTime_NaiveMLE(New_Y, dat$Data_Y, dat$Data_time) 
plotPredictions(pred_NaiveMLE,TrueTimes)
#plotPredictions(pred_NaiveMLE,TrueTimes,pdfFileName="PredWith_NaiveMLE")

#Specify hyperparameters of priors:
Beta0 = c(2,0) 
Lambda0= 10
Sigma0 = 1
pred_FullBayesian = posteriorTime_FullBayesian(New_Y, dat$Data_Y, dat$Data_time,Lambda0=Lambda0, Sigma0=Sigma0, Beta0=Beta0) 
plotPredictions(pred_FullBayesian,TrueTimes)
#plotPredictions(pred_FullBayesian,TrueTimes,pdfFileName="PredWith_FullBayesian")

#Calibrate hyperparameters based on marginal evidence (emprical bayes), Used fixed Beta
pred_EB1 = posteriorTime_EmpiricalBayesian(New_Y, dat$Data_Y, dat$Data_time, Beta0=Beta0) 
plotPredictions(pred_EB1,TrueTimes)

#Calibrate hyperparameters based on marginal evidence (emprical bayes)
pred_EB2 = posteriorTime_EmpiricalBayesian(New_Y, dat$Data_Y, dat$Data_time) 
plotPredictions(pred_EB2,TrueTimes)











