#' @title genData
#' @description Generating a example dataset (TSD)
#' @details The function generates a dataset based on multivariate normal distribution where time (uniformly distributed) is a underlying explanatory variable
#' 
#' @param ntot total number of observations (independent)
#' @param nsites number of sites (multivariate response)
#' @param beta Slope parameters (effect on time) for each sites 
#' @param SIGMA Covariance matrix 
#' @param times Time variable for each individuals  
#' @param maxTime maximum of time (upper limit)
#' @param seed The user can set seed if wanted
#' @return A list with variables Data_Y (response), Data_Time 
#' @importFrom stats rnorm runif
#' @export
#' @examples
#' \dontrun{ 
#' dat = genData(seed=1)
#' }
genData = function(ntot = 50, nsites = 10,  beta=NULL, SIGMA = NULL, times=NULL, maxTime = 12, seed=NULL) {
  if(!is.null(seed)) set.seed(seed) #set seed if provided

   #GENERATE Time stamps
  if(is.null(times)) {
    times = runif(ntot, 0, maxTime) #obtain times
  } else {
    if(length(times)!=ntot) stop("times argument must be a ntot long vector")
  }
  

  if(is.null(beta)) {
    beta0 = rnorm(nsites,5,1) #simulate default vector if not provided
    beta1 = rnorm(nsites,0,1) #simulate default vector if not provided
    beta = rbind(beta0,beta1)
  } else {
    dim = dim(beta)
    if(is.null(dim) || dim[1]!=2 || dim[2]!=nsites) stop("beta argument must be a [2 x nsites] matrix")
  }

  if(is.null(SIGMA)) {
    SIGMA = diag(2,nsites) # default matrix (independent sites)
  } else {
    dim = dim(SIGMA)
    if(is.null(dim) || dim[1]!=dim[2] || dim[1]!=nsites) stop("SIGMA argument must be a [nsites x nsites] matrix")
  }
  
  #Simulate responses
  Chol = t(chol(SIGMA)) #obtain cholesky decomposition of SIGMA
  EPS = matrix(rnorm(ntot),ncol=ntot,nrow=nsites) #create random noise
  Y = beta[1,] + outer(beta[2,],times) + Chol%*%EPS #generate response (nsamples x nsites)
  rownames(Y) = paste0("Site",1:nrow(Y))
  
  dat = list(Data_Y=t(Y), Data_time=times, beta=beta,SIGMA=SIGMA) #store data in dataframe
  #dat = data.frame(age=agej,batch=as.factor(batchInd),beta=beta) #store data in dataframe

  return(dat)
}
  