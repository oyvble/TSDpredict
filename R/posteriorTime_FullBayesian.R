
#' @title posteriorTime_FullBayesian
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
#' @importFrom stats rnorm runif
#' @return ret list with posterior distribution for each new individuals (separate list for univariate prediction and all combined)
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
#' predObj = posteriorTime_FullBayesian(New_Y, Data_Y, Data_time)
#' plotPredictions(predObj)
#' }
posteriorTime_FullBayesian = function(New_Y, Data_Y, Data_time , nu0 = 1, Lambda0=1, Sigma0 = 1, Beta0=NULL, maxTime=12,delta=0.01, priorTime=NULL) {
    
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
  
  #Hyperparameter settings for parameter priors:
  if(length(Lambda0)==1) {
    Lambda0 = diag(Lambda0,m) #hyperparameter for wishart prior
  } else {
    dim = dim(Lambda0)
    if(is.null(dim) || dim[1]!=m || dim[1]!=dim[2]) stop("Lambda0 parameter must be a [nSites x nSites] square matrix")
  }
  #Lambda0 = matrix(lam0/100,m,m) #hyperparameter for wishart prior
  #diag(Lambda0) = lam0
  if(length(Sigma0)==1) { #if only a numerical given
    Sigma0 = c(rep(Sigma0,2),0)
  } else if(length(Sigma0)==2) { #if only a numerical given
    Sigma0 = c(Sigma0,0)
  } else if(length(Sigma0)!=3) { #otherwise the matrix must be given (3 long vector)
    stop("Sigma0 parameter must be a 3 long vector")
  } 
  Sigma0 = matrix(c(Sigma0[1],Sigma0[3],Sigma0[3],Sigma0[2]),2,2)

  if(is.null(Beta0)) { #if prior not given
    Beta0 = matrix(0,nrow=2,ncol=m)
  } else if(length(Beta0)==2) { #if a two long vector
    Beta0 = replicate(m,Beta0) #assume same mean for each sites
  }  else {
    dim = dim(Beta0)
    if(is.null(dim) || dim[1]!=2 || dim[2]!=m) stop("Beta parameter must be a [nSites x nSites] square matrix")
  }
  
  #CALCULATING NECESSARY VARIABLES FOR MODEL
  X = cbind(1, Data_time)
  XtX = t(X)%*%X #empirical precision matrix
  projY = solve(XtX)%*%(t(X)) #projector matrix: Projecting data Y into space of Xbeta
  betahat = projY%*% Data_Y #obtain OLS for coefficients based on data
  XtX_betahat = XtX%*%betahat
  
  invSigma0 = solve(Sigma0) #inverse of Sigma0 (prior)
  invSigma0_Beta0 = invSigma0%*%Beta0 #taking into account prior mean of beta
  Beta0t_invSigma0_Beta0 = t(Beta0)%*%invSigma0%*%Beta0 #taking into account prior mean of beta
  XtX_betahat_invSigma0_Beta0 = XtX_betahat + invSigma0_Beta0 #temporary variable taking into account prior mean of beta
  
  XtX_invSigma0 = XtX + invSigma0 #precalculate
  inv_XtX_invSigma0 = solve(XtX_invSigma0) #take the inverse
  betaStar =  inv_XtX_invSigma0%*%XtX_betahat_invSigma0_Beta0 #posterior betahat [2 x m] matrix
  Astar = t(Data_Y)%*% Data_Y + Beta0t_invSigma0_Beta0 - t(XtX_betahat_invSigma0_Beta0)%*%inv_XtX_invSigma0%*%XtX_betahat_invSigma0_Beta0 #m x m matrix (residual sum square matrix)
  AstarPrior = Astar + Lambda0 #add Prior
  nu = n + nu0 #degrees of freedom
  konst = lgamma( (nu+m)/2 ) - lgamma( nu/2 ) - m/2*log(pi)#*nu) #m dimension constant
  konst1 = lgamma( (nu+1)/2 ) - lgamma( nu/2 ) - 1/2*log(pi)#*nu) #1 dimension
  
  #Calculate C scalar for each time value:
  Cvec = rep(NA,length(time_grid)) #this is inverse of the predictive variance
  for(t in 1:length(Cvec)) {
    time = time_grid[t]
    x0 = c(1,time)
    x0tx0 = x0%*%t(x0)
    invV = solve( x0tx0 + XtX_invSigma0 ) #obtain invert of V matrix (Varaince of the predictor)
    Cvec[t] = 1 - x0%*%invV%*%x0
  }

  ##############################  
  #CALCULATING TIME PREDICTIONS#
  ##############################  
  
  #CALCULATE MULTIVARIATE PREDICTION (possibly nNew samples):
  inv_AstarPrior = solve(AstarPrior) #take inverse of this as pre-step
  logdet_AstarPrior = determinant(AstarPrior)$mod[1] #log determinant of AstarPrior scale
  
  pdf_MULTI = matrix(NA,nrow=nNew,ncol=length(time_grid),dimnames = list(rownames(New_Y),time_grid))
  for( t in 1:length(time_grid)) {
    #Scale = AstarPrior/Cvec[t] #would be kroenicker product in high dimension (this is matrix multiplication A x C)
    #logdetScale = determinant(Scale)$mod[1] #log determinant of scale
    invScale = Cvec[t]*inv_AstarPrior # #calculating inverse scale
    logdetScale = logdet_AstarPrior - m*log(abs(Cvec[t])) #remember to scale with dimension
    
    Exp = betaStar[1,] +  betaStar[2,]*time_grid[t] #Obtain location param
    res = t(New_Y) - Exp #this is residual (obs - exp)

    logval = konst - 0.5*logdetScale #obtaining common logval part (of predictions)
    for(i in 1:nNew) { #for each new predictions
      pdf_MULTI[i,t] = exp(logval - ((nu + m)/2)*log(1 + res[,i]%*%invScale%*%res[,i])) #transform
    }
    if(!is.null(priorTime)) pdf_MULTI[,t] = pdf_MULTI[,t]*priorTime(time_grid[t]) #scale pdf with prior of time
  }
  pdf_MULTI = pdf_MULTI/(delta*rowSums(pdf_MULTI)) #THIS IS POSTERIOR OF MULTIVARIATE PREDICTION
  #sapply(rowSums(pdf_MULTI*delta),function(x) all.equal(x,1))
 # plot(time_grid,pdf_MULTI,ty="l")
  
  #CALCULATE MARGINAL PREDICTIONS (still T-distribution):
  pdf_UNIsite = list() #obtain univariate predictions for each site
  for(j in 1:m) { #for each site
    pdf_UNI = matrix(NA,nrow=nNew,ncol=length(time_grid),dimnames = list(rownames(New_Y),time_grid))
    Exp = betaStar[1,j] +  betaStar[2,j]*time_grid #Obtain location param
    
    Scale = AstarPrior[j,j]/Cvec #kroenicker product in high dimension (this is the shape matrix)
    logdetScale = log(Scale) #log determinant of scale
    
    logval = konst1 - 0.5*logdetScale #obtaining common logval part (of predictions)
    for(i in 1:nNew) { #for each new predictions
      res = New_Y[i,j] - Exp #this is residual (obs - exp)
      pdf_UNI[i,] = exp( logval - ((nu + 1)/2)* log(1 + res^2/Scale) ) #transform
      if(!is.null(priorTime)) pdf_UNI[i,] = pdf_UNI[i,]*priorTime(time_grid) #scale pdf with prior of time
    }
    pdf_UNIsite[[sites[j]]] =  pdf_UNI/(delta*rowSums(pdf_UNI))
  # lines(time_grid, pdf_UNIsite[[sites[j]]][i,],ty="l")
  }
  return(list(pdf_MULTI= pdf_MULTI,pdf_UNIsite=pdf_UNIsite))
}

