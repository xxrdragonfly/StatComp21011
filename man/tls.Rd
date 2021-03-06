% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TLS.R
\name{tls}
\alias{tls}
\title{A method to estimate paramaters for spatial threshold autoregressive model}
\usage{
tls(X, y, W, q)
}
\arguments{
\item{X}{a matrix of exogenous regressors}

\item{y}{a vector of observation on the dependent variable}

\item{W}{the row normalized spatial weight matrix}

\item{q}{threshold variable which is assumed exogenous with a continuous distribution}
}
\value{
estimation of parameters
}
\description{
A method to estimate paramaters for spatial threshold autoregressive model
}
\examples{
\dontrun{
  n<-nrow(iris)
  W<-matrix(0,n,n)
  for (i in 1:n) {
    for (j in 1:n) {
     if(abs(j-i)<=3)
     W[i,j]<-1
    }
  }
  W <- sweep(W, 1, apply(W, 1, function(x) sqrt(sum(x^2))), "/")
  est <- tls(iris[,2:3],iris[,1],W,iris[,4])
  print(est)
}
}
