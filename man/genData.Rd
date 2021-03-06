% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genData.R
\name{genData}
\alias{genData}
\title{genData}
\usage{
genData(
  ntot = 50,
  nsites = 10,
  beta = NULL,
  SIGMA = NULL,
  times = NULL,
  maxTime = 12,
  seed = NULL
)
}
\arguments{
\item{ntot}{total number of observations (independent)}

\item{nsites}{number of sites (multivariate response)}

\item{beta}{Slope parameters (effect on time) for each sites}

\item{SIGMA}{Covariance matrix}

\item{times}{Time variable for each individuals}

\item{maxTime}{maximum of time (upper limit)}

\item{seed}{The user can set seed if wanted}
}
\value{
A list with variables Data_Y (response), Data_Time
}
\description{
Generating a example dataset (TSD)
}
\details{
The function generates a dataset based on multivariate normal distribution where time (uniformly distributed) is a underlying explanatory variable
}
\examples{
\dontrun{ 
dat = genData(seed=1)
}
}
