% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RMBCaux.R
\name{RMBCaux}
\alias{RMBCaux}
\title{RMBCaux}
\usage{
RMBCaux(
  Y,
  K,
  thetaOld.alpha,
  thetaOld.mu,
  thetaOld.sigma,
  max_iter,
  niterFixedPoint,
  tolerance,
  cutoff = 1 - 0.001
)
}
\arguments{
\item{Y}{A matrix of size n x p.}

\item{K}{The number of clusters.}

\item{thetaOld.alpha}{The initial alpha: An array of K positive real numbers
they must verify the condition sum(thetaOld.mu)== 1.}

\item{thetaOld.mu}{The initial centers: A list with K elements, each of them 
is an array of length p.}

\item{thetaOld.sigma}{The initial stcatter matrix: A list with K matrix, each of them 
has dimension p x p}

\item{max_iter}{a maximum number of iterations used for the
algorithm stopping rule}

\item{niterFixedPoint}{the maximum number of iteration in the internal loop which 
computes sigma an mu separately. The default value is niterFixedPoint=1}

\item{tolerance}{tolerance parameter used for the algorithm stopping
rule}

\item{cutoff}{optional argument for outliers detection - quantiles of chi-square 
to be used as a threshold for outliers detection, defaults to 0.999}
}
\value{
A list including the estimated K centers and labels for the
    observations

\itemize{
\item{\code{centers}}{: matrix of size K
    x p, with the estimated K centers.}
\item{\code{cluster}}{:
    array of size n x 1 integers labels between 1 and K.}
\item{\code{tauPath}}{: sequence of tau scale values at each
    iterations.}
\item{\code{Wni}}{: numeric array of size n x 1
    indicating the weights associated to each observation.}
\item{\code{emptyClusterFlag}}{: a boolean value. True means
    that in some iteration there were clusters totally empty}
\item{\code{niter}}{: number of iterations until convergence
    is achived or maximum number of iteration is reached}
\item{\code{di}}{distance of each observation to its assigned
    cluster-center} 
    }
}
\description{
Robust Model Base Clustering algorithm based on centers, a robust and
efficient version of EM algorithm.
}
