
#' @title posteriorTime_FullBayesian
#' @author Oyvind Bleka
#' @description Calculating the posterior distribution of time for new individuals
#'
#' @details This function provides posterior distribution of time, based on a multivariate model (Naive MLE approach), where time as a continuous exploratory variable from zero to maxTime
#'
#' Prediction: Applying Bayes Theorem: p(time|Y) = konstant x p(Y|time) x p(time)
#' We assume an uniform prior for p(log(time))=unif, gives p(time)=1/t
#'
#'  Model:
#'  p(Y|time,theta)=MVN(b0 + b1*time, SIGMA), SIGMA= COVARIANCE MATRIX
#'  where theta is estimated using maximum likelihood estimation (MLE), provided as theta_hat
#'
#'  Predictor p(Ynew|time)=p(Ynew|time,theta=theta_hat)
#'
#' @param New_Y The observed response for new individuals (nsamples x nsites) matrix (also possible with vector)
#' @param Data_Y This is observed multi-responses (nsamples x nsites)
#' @param Data_time This is observed times (nsamples vector)
#' @param maxTime maximum of time (upper limit)
#' @param delta Grid size for the time variable
#' @param priorTime A prior distribution for the time variable (NULL means uniform distribution)
#' @importFrom mvtnorm dmvnorm
#' @importFrom stats coef cov dnorm lm
#' @return ret list with posterior distribution for each new individuals (separate list for univariate prediction and all combined)
#' @export
#' @examples 
#' \dontrun{ 
#' ntrain = 100
#' ntest = 100
#' dat = genData(ntrain+ ntest,seed=1)
#' Data_Y = dat$Data_Y[1:ntrain,]
#' Data_time = dat$Data_time[1:ntrain]
#' New_Y = dat$Data_Y[-(1:ntrain),]
#' New_time = dat$Data_time[-(1:ntrain)]
#' predObj = posteriorTime_NaiveMLE(New_Y, Data_Y, Data_time)
#' }
posteriorTime_NaiveMLE = function(New_Y, Data_Y, Data_time, maxTime=12,delta=0.01, priorTime=NULL) {
    
  #Define grid for time variable
  time_grid = seq(0,maxTime,delta)
  
  sites = colnames( Data_Y) #obtain names of sites
  m = length(sites) #number of variables
  n = nrow( Data_Y) #number of samples
  nNew = nrow(New_Y) #number of new data to predict
  if(is.null(nNew)) {
    nNew = 1
    New_Y = rbind(New_Y) #convert to matrix
  }
  
  betahat = matrix(nrow=2,ncol=m) #fitted coefficients
  res <- matrix(nrow=n,ncol=m,dimnames=list(NULL,sites)) #calculate residuals
  #resRE Residuals based on Random effect model (not considered)
  
  for(j in 1:m) { #traverse each sites marginally
    fit = lm(Data_Y[,j]~Data_time) #fit ordinary linear regression model (OLS)
    betahat[,j] = coef(fit)  #estimated coefficients based on fixed model
    res[,j] = fit$res #insert residuals
  }
  #Calculate Sample-Covariance matrix of residuals:
  SIGMA = cov(res) #t(res)%res/(n-2)
  
  
  #CALCULATE MULTIVARIATE PREDICTION (possibly nNew samples):
  pdf_MULTI = matrix(NA,nrow=nNew,ncol=length(time_grid),dimnames = list(rownames(New_Y),time_grid))
  for( t in 1:length(time_grid)) {
    Exp = betahat[1,] +  betahat[2,]*time_grid[t] #Obtain location param
    pdf_MULTI[,t] = dmvnorm(New_Y, Exp, SIGMA) #insert each predictors
    if(!is.null(priorTime)) pdf_MULTI[,t] = pdf_MULTI[,t]*priorTime(time_grid[t]) #scale pdf with prior of time
  }
  pdf_MULTI = pdf_MULTI/(delta*rowSums(pdf_MULTI)) #THIS IS POSTERIOR OF MULTIVARIATE PREDICTION
  #all(sapply(rowSums(pdf_MULTI*delta),function(x) all.equal(x,1)))
     
  
  #CALCULATE MARGINAL PREDICTIONS (still T-distribution):
  pdf_UNIsite = list() #obtain univariate predictions for each site
  for(j in 1:m) { #for each site
    pdf_UNI = matrix(NA,nrow=nNew,ncol=length(time_grid),dimnames = list(rownames(New_Y),time_grid))
    Exp = betahat[1,j] +  betahat[2,j]*time_grid #Obtain location param
    sigma = SIGMA[j,j]
    for(i in 1:nNew) { #for each new predictions
      pdf_UNI[i,] =  dnorm(New_Y[i,j], Exp,sigma) #insert each predictors
      if(!is.null(priorTime)) pdf_UNI[i,] = pdf_UNI[i,]*priorTime(time_grid) #scale pdf with prior of time
    }
    pdf_UNIsite[[sites[j]]] =  pdf_UNI/(delta*rowSums(pdf_UNI))
  }
  return(list(pdf_MULTI= pdf_MULTI,pdf_UNIsite=pdf_UNIsite))
}

