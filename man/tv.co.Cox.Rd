% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tv.co.Cox.R
\name{tv.co.Cox}
\alias{tv.co.Cox}
\title{Multiplicative hazards model with sparse longitudinal covariates}
\usage{
tv.co.Cox(data, n, l, times, bd, scb)
}
\arguments{
\item{data}{An object of class tibble. The structure of the tibble must be: tibble(id = id, X = failure time, covariates = observation for covariates, obs_times = observation times, delta = censoring indicator).}

\item{n}{An object of class integer. The sample size.}

\item{l}{An object of class vector. The selection vector. For example, for the p dimensional regression coefficient function, if we want to construct simultaneous confidence band for the first regression coefficient function, we can take l=c(1,0,...,0).}

\item{times}{An object of class vector. The interest time.}

\item{bd}{An object of class vector. If use auto bandwidth selection, the structure of the vector must be: bd=c(the maximum bandwidth, the minimum bandwidth, the number of bandwidth divided). If use fixed bandwidth, bd is the chosen bandwidth.}

\item{scb}{An object of class vector. If need to construct the simultaneous confidence band, the structure of the vector must be: c(desirable confidence level, repeat times). Otherwise, scb=0.}
}
\value{
a list with the following elements:
\item{est}{The estimation for the corresponding parameters.}
\item{se}{The estimation for the standard error of the estimated parameters.}
\item{scb}{The quantile used to construct simultaneous confidence band.}
}
\description{
Regression analysis of multiplicative hazards model with sparse longitudinal covariates. The kernel weighting approach is employed to impute the missing value and localize the estimating equation. A wild bootstrap-based simultaneous confidence band for the nonparametric function is also provided.
}
\examples{

library(dplyr)
library(gaussquad)
library(MASS)
library(nleqslv)
n=500
beta<-function(t){
  0.5*(t+0.5)^2
}
lqrule64 <- legendre.quadrature.rules(64)[[64]]
simdata <- function(  beta ) {
cen=1
nstep=20
Sigmat_z <- exp(-abs(outer(1:nstep, 1:nstep, "-")) / nstep)
z  <-c(mvrnorm(  1, rep(0,20), Sigmat_z  ))
left_time_points <- (0:(nstep - 1)) / nstep
z_fun <- stepfun(left_time_points, c(0,z  ))
h_fun <- function(x) { beta(x) * z_fun(x)  }
lam_fun <- function(tt) 2 * exp(h_fun(tt))
u <- runif(1)
fail_time <- nleqslv(0, function(ttt)
legendre.quadrature(lam_fun, lower = 0,upper = ttt, lqrule64) + log(u))$x
X <- min(fail_time, cen)
obs=rpois(1,  5)+1
tt= sort(runif(obs, min = 0, max = 1))
obs_times <- tt[which(tt<=cen)]
if (length(obs_times) == 0)
 obs_times <- cen
 covariates_obscov <-z_fun(obs_times)
 return( tibble(X = X,delta = fail_time < cen,
  covariates = covariates_obscov,obs_times = obs_times   )  ) }


data <- replicate(  n, simdata( beta ), simplify = FALSE ) \%>\% bind_rows(.id = "id")
tv.co.Cox(data,n,l,0.2,bd=c(n^(-0.4),n^(-0.4)),scb=0)
}
\references{
Sun, Z. and Cao, H. (2023)  <arXiv:2310.15877>
}
