
#' @title MarginalEvid
#' @author Oyvind Bleka
#' @description Calculating the marginal evidence for the data
#'
#' @details This function provides the marginal evidence value based on all the data
#'
#'  Model: A full Bayesian model (takes into account the uncertainty of the parameters)
#'  p(Y|time) = int_theta p(Y|time,theta)p(theta) d_theta
#'  Where
#'  p(Y|time,theta)=MVN(b0 + b1*time, SIGMA), SIGMA= COVARIANCE MATRIX
#'  with latent priors 'normal-invWishart'
#'  p(b0,b1|SIGMA)=MVN(0,0, diag(sig0 SIGMA,sig1 SIGMA))
#'  p(SIGMA) = invWishart(Lambda0,nu0), where diag(Lambda0) = lam0
#'
#' Marginal evidence p(Y|eta) = int_theta p(Y|theta)p(theta|eta)d theta
#'
#' @param Data_Y This is observed multi-responses (nsamples x nsites)
#' @param Data_time This is observed times (nsamples vector)
#' @param nu0  degrees of freedom hyperparameter of invWishart prior distribution (should be 0=non-informative or 1=informative)
#' @param Lambda0  Scale matrix hyperparameter of invWishart prior distribution of Covariance (can be numeric or matrix)
#' @param Sigma0 Covariance hyperparameter of Normal prior distribution of beta (can be numeric, 2 or 3 long vector)
#' @param Beta0  Mean hyperparameter of Normal prior distribution of beta (can be 2 long vector or matrix)
#' @return The numerical marginal evidence
#' @export
#' @examples 
#' \dontrun{ 
#' n = 100
#' dat = genData(n,seed=1)
#' Data_Y = dat$Data_Y
#' Data_time = dat$Data_time
#' modelEvid = MarginalEvid_FullBayesian(Data_Y, Data_time)
#' }
MarginalEvid_FullBayesian = function(Data_Y, Data_time , nu0 = 1, Lambda0=1, Sigma0 = 1, Beta0=NULL) {
    
  sites = colnames( Data_Y) #obtain names of sites
  m = length(sites) #number of variables
  n = nrow( Data_Y) #number of samples
  
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
  Astar = t(Data_Y)%*% Data_Y + Beta0t_invSigma0_Beta0 - t(XtX_betahat_invSigma0_Beta0)%*%inv_XtX_invSigma0%*%XtX_betahat_invSigma0_Beta0 #m x m matrix (residual sum square matrix)
  AstarPrior = Astar + Lambda0 #add Prior
  
  nu = n + nu0 #degrees of freedom
  logdet_AstarPrior = determinant(AstarPrior)$mod[1] #log determinant of AstarPrior scale
  logdet_Lambda0 = determinant(Lambda0)$mod[1] #log determinant of Lambda0
  logdet_Sigma0 = determinant(Sigma0)$mod[1]  #log determinant of Sigma0
  logdet_XtX_invSigma0 = determinant(XtX_invSigma0)$mod[1]  #log determinant of XtX + invSigma0
  
  ###################################
  #Calculate model evidence (logged)#
  ###################################
  
  multilogGamma = function(x) { #multivariate gamma function
    m*(m-1)/4*log(pi) + sum(lgamma(x+(1-1:m)/2))
  }
  
  val1 = -0.5*m*n*log(2*pi)
  val2 = 0.5*m*(2*nu0 + 2*m + n - 1)*log(2)
  val3 = multilogGamma( (nu+m-1)/2 ) - multilogGamma( (nu0+m-1)/2 )
  val4 = 0.5*( (nu0 + m -1)*logdet_Lambda0 - (nu + m -1)*logdet_AstarPrior )
  val5 = -0.5*m*(  logdet_Sigma0 + logdet_XtX_invSigma0 )
    
  val = val1 + val2 + val3 + val4 + val5
  return(val)
}

