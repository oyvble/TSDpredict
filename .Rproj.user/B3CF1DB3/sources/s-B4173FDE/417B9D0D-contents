#library(TSDpredict)
#demo(simulatePredictionExample)
#Simulating examples to predict age

ntrain = 100
ntest = 100 #repeat test
time = 6
TrueTimes = rep(time,ntest)

nsites = 8
#paste0(round(beta,1),collapse=",")
beta = c(5,4.9,5.7,6,4.3,5.6,5.9,5.1,-0.9,-0.8,-1.2,-1,-0.6,-0.9,-0.9,-0.9)
beta = t(matrix(beta,ncol=2))
SIGMA = c(1.1,0.2,0.4,0.8,0.5,0.2,0.1,0.5,0.2,0.8,0.4,0.5,0.5,0.4,0.3,0.4,0.4,0.4,0.6,0.8,0.2,0.2,0.2,0.5,0.8,0.5,0.8,1.7,0.4,0.1,0.4,0.7,0.5,0.5,0.2,0.4,1.2,0.7,0.4,0.2,0.2,0.4,0.2,0.1,0.7,0.9,0.3,0.2,0.1,0.3,0.2,0.4,0.4,0.3,0.5,0.2,0.5,0.4,0.5,0.7,0.2,0.2,0.2,0.5)
SIGMA = matrix(SIGMA,ncol=nsites)

dat = genData(ntrain,nsites,beta,SIGMA,maxTime = 12,seed=1)

#Show train data:
times = dat$Data_time
sites = colnames(dat$Data_Y)
  
par(mfrow=c(nsites/2,2))
for(j in 1:nsites) {
  plot( times, dat$Data_Y[,j],main=sites[j],xlab="Time",ylab="Response (log (Rfu+1)")
}

test = genData(ntest,nsites,beta,SIGMA,times = TrueTimes,seed=2)
New_Y = test$Data_Y #responses to test with

pred_NaiveMLE = posteriorTime_NaiveMLE(New_Y, dat$Data_Y, dat$Data_time,pdfFileName="PredWith_NaiveMLE") 
plotPredictions(pred_NaiveMLE,TrueTimes)

pred_FullBayesian = posteriorTime_FullBayesian(New_Y, dat$Data_Y, dat$Data_time,pdfFileName="PredWith_FullBayesian") 
plotPredictions(pred_FullBayesian,TrueTimes)
