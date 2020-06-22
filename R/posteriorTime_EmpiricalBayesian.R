
#' @title posteriorTime_EmpiricalBayesian
#' @author Oyvind Bleka
#' @description Calculating the posterior distribution of time for new individuals
#'
#' @details This function provides posterior distribution of time, based on a multivariate model (Full Bayesian approach), where time as a continuous exploratory variable from zero to maxTime
#'
#' Prediction: Applying Bayes Theorem: p(time|Y) = konstant x p(Y|time) x p(time)
#' We assume an uniform prior for p(log(time))=unif, gives p(time)=1/t
#'
#'  Model: A full Bayesian model (takes into account the uncertainty of the parameters)
#'  p(Y|time) = int_theta p(Y|time,theta)p(theta) d_theta
#'  Where
#'  p(Y|time,theta)=MVN(b0 + b1*time, SIGMA), SIGMA= COVARIANCE MATRIX
#'  with latent priors 'normal-invWishart'
#'  p(b0,b1|SIGMA)=MVN(0,0, diag(sig0 SIGMA,sig1 SIGMA))
#'  p(SIGMA) = invWishart(Lambda0,nu0), where diag(Lambda0) = lam0
#'
#' Posterior predictor p(Ynew|time) is a m-Multivariate T-distribution(M,A,nu):
#' p(Y) = gamma((nu+m)/2)/gamma(nu/2) (pi*nu)^(-m/2) |A|^-0.5 (1 + (Y-M)'A^{-1}(Y-M)/nu)^-(nu+m)/2
#'
#' @param New_Y The observed response for new individuals (nsamples x nsites) matrix (also possible with vector)
#' @param Data_Y This is observed multi-responses (nsamples x nsites)
#' @param Data_time This is observed times (nsamples vector)
#' @param nu0  degrees of freedom hyperparameter of invWishart prior distribution (should be 0=non-informative or 1=informative)
#' @param Lambda0  Scale matrix hyperparameter of invWishart prior distribution of Covariance (can be numeric or matrix)
#' @param Sigma0 Covariance hyperparameter of Normal prior distribution of beta (can be numeric, 2 or 3 long vector)
#' @param Beta0  Mean hyperparameter of Normal prior distribution of beta (can be 2 long vector or matrix)
#' @param maxTime maximum of time (upper limit)
#' @param delta Grid size for the time variable
#' @param priorTime A prior distribution for the time variable (NULL means uniform distribution)
#' @param optsettings A list of settings for optimizing the marginal evidence wrt hyperparam
#' @return ret list with posterior distribution for each new individuals (separate list for univariate prediction and all combined)
#' @importFrom stats nlm rnorm
#' @export
#' @examples 
#' \dontrun{ 
#' ntrain = 100
#' ntest = 1
#' dat = genData(ntrain+ ntest,seed=1)
#' Data_Y = dat$Data_Y[1:ntrain,]
#' Data_time = dat$Data_time[1:ntrain]
#' New_Y = dat$Data_Y[-(1:ntrain),]
#' New_time = dat$Data_time[-(1:ntrain)]
#' predObj = posteriorTime_EmpiricalBayesian(New_Y, Data_Y, Data_time)
#' plotPredictions(predObj)
#' }
posteriorTime_EmpiricalBayesian = function(New_Y, Data_Y, Data_time , nu0 = 1, Lambda0=NULL, Sigma0 = NULL, Beta0 = NULL, maxTime=12,delta=0.01, priorTime=NULL, optsettings = list(noptim=4,optimsd=1) ) {
    

  sites = colnames( Data_Y) #obtain names of sites
  m = length(sites) #number of variables
  n = nrow( Data_Y) #number of samples

  #Boolean of wheter to calibrate
  calibrate = rep(F,3)
  names(calibrate) = c("lambda0","sigma0","beta0")
  calibrate[1] = ifelse(is.null(Lambda0),T,F)
  calibrate[2] = ifelse(is.null(Sigma0),T,F)
  calibrate[3] = ifelse(is.null(Beta0),T,F)
  
  #CALIBRATE HYPERPARAMETERS BASED ON MarginalEvid
  lam0 = Lambda0 #store as temporary variable (if not provided)
  sig0 = Sigma0 #store as temporary variable (if not provided)
  beta0 = Beta0 #store as temporary variable (if not provided)

  negMargEvid = function(eta)  {
    indCounter = 1 #index counter
    if(calibrate[1]) {
      lam0 = exp(eta[indCounter])
      indCounter = indCounter + 1
    }
    if(calibrate[2]) {
      sig0 = exp(eta[indCounter])
      indCounter = indCounter + 1
    }
    if(calibrate[3]) {
      beta0 = eta[c(indCounter,indCounter+1)]
    }
    val = MarginalEvid_FullBayesian(Data_Y, Data_time,nu0=nu0, Lambda0=lam0, Sigma0=sig0, Beta0=beta0)
    return(-val)
  }
  
  #Prepare optimizing;
  eta_dim = sum(calibrate)
  if(calibrate[3]) eta_dim = eta_dim + 1
  
  bestList = list(min=Inf)
  for(i in 1:optsettings$noptim) {
    foo = nlm(negMargEvid,rnorm(eta_dim,0,optsettings$optimsd))
    if(bestList$min > foo$min) bestList = foo
  }
  eta = bestList$estimate #obtain estimated

  #extract params  
  indCounter = 1 #index counter
  if(calibrate[1]) {
    lam0 = exp(eta[indCounter])
    indCounter = indCounter + 1
  }
  if(calibrate[2]) {
    sig0 = exp(eta[indCounter])
    indCounter = indCounter + 1
  }
  if(calibrate[3]) {
    beta0 = eta[c(indCounter,indCounter+1)]
  }  
  
  return( posteriorTime_FullBayesian(New_Y, Data_Y, Data_time , nu0 = nu0, Lambda0=lam0, Sigma0 = sig0, Beta0=beta0, maxTime=maxTime,delta=delta, priorTime=priorTime) )
}

