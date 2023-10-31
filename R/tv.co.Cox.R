#' @title Multiplicative hazards model with sparse longitudinal covariates
#'
#' @description Regression analysis of multiplicative hazards model with sparse longitudinal covariates. The kernel weighting approach is employed to impute the missing value and localize the estimating equation. A wild bootstrap-based simultaneous confidence band for the nonparametric function is also provided.
#'
#' @param  data  An object of class tibble. The structure of the tibble must be: tibble(id = id, X = failure time, covariates = observation for covariates, obs_times = observation times, delta = censoring indicator).
#'
#' @param  n An object of class integer. The sample size.
#'
#' @param  l An object of class vector. The selection vector. For example, for the p dimensional regression coefficient function, if we want to construct simultaneous confidence band for the first regression coefficient function, we can take l=c(1,0,...,0).
#'
#' @param  times An object of class vector. The interest time.
#'
#' @param  bd An object of class vector. If use auto bandwidth selection, the structure of the vector must be: bd=c(the maximum bandwidth, the minimum bandwidth, the number of bandwidth divided). If use fixed bandwidth, bd is the chosen bandwidth.
#'
#' @param  scb An object of class vector. If need to construct the simultaneous confidence band, the structure of the vector must be: c(desirable confidence level, repeat times). Otherwise, scb=0.
#'
#' @references Sun, Z. and Cao, H. (2023)  <arXiv:2310.15877>
#'
#' @return a list with the following elements:
#' \item{est}{The estimation for the corresponding parameters.}
#' \item{se}{The estimation for the standard error of the estimated parameters.}
#' \item{scb}{The quantile used to construct simultaneous confidence band.}
#'
#' @import  gaussquad
#' @importFrom nloptr nloptr
#' @importFrom nleqslv nleqslv
#' @importFrom dplyr %>% bind_rows group_by pull reframe filter
#' @importFrom stats quantile lm rexp
#' @examples
#'
#' library(dplyr)
#' library(gaussquad)
#' library(MASS)
#' library(nleqslv)
#'n=500
#'beta<-function(t){
#'   0.5*(t+0.5)^2
#'}
#'lqrule64 <- legendre.quadrature.rules(64)[[64]]
#'simdata <- function(  beta ) {
#' cen=1
#' nstep=20
#' Sigmat_z <- exp(-abs(outer(1:nstep, 1:nstep, "-")) / nstep)
#' z  <-c(mvrnorm(  1, rep(0,20), Sigmat_z  ))
#' left_time_points <- (0:(nstep - 1)) / nstep
#' z_fun <- stepfun(left_time_points, c(0,z  ))
#' h_fun <- function(x) { beta(x) * z_fun(x)  }
#' lam_fun <- function(tt) 2 * exp(h_fun(tt))
#' u <- runif(1)
#' fail_time <- nleqslv(0, function(ttt)
#' legendre.quadrature(lam_fun, lower = 0,upper = ttt, lqrule64) + log(u))$x
#' X <- min(fail_time, cen)
#' obs=rpois(1,  5)+1
#' tt= sort(runif(obs, min = 0, max = 1))
#' obs_times <- tt[which(tt<=cen)]
#' if (length(obs_times) == 0)
#'  obs_times <- cen
#'  covariates_obscov <-z_fun(obs_times)
#'  return( tibble(X = X,delta = fail_time < cen,
#'   covariates = covariates_obscov,obs_times = obs_times   )  ) }
#'
#'
#' data <- replicate(  n, simdata( beta ), simplify = FALSE ) %>% bind_rows(.id = "id")
#' tv.co.Cox(data,n,l,0.2,bd=c(n^(-0.4),n^(-0.4)),scb=0)
#' @export
#'
#'

tv.co.Cox <- function(data,n,l,times, bd, scb) {

  nt=length(times)
  X=data$X
  id <- data$id
  covariates <- matrix(c(data$covariates),length(X))
  obs_times <- data$obs_times
  delta <- data$delta
  p_z <- dim(covariates)[2]
  if (is.null(p_z)) p_z <- 1
  outf=function(x){x %o% x}
  outmf=function(x,y){x %o% y}
  kerfun <- function(xx){
    pmax((1-xx^2)*0.75,0)
  }

  estproc_Cox_test <- function(data, n,s, h1,h2) {
    X=data$X
    id <- data$id
    covariates <- data$covariates
    obs_times <- data$obs_times
    delta <- data$delta
    kerval <- kerfun((X - s) / h1) * kerfun((obs_times-s)/h2) / h1 /h2


    # Estimating equation and estimation

    Delta<-outer(X,X,"<=")
    kerval_XX<- kerfun((X-s)/h1) %*%t(kerfun((obs_times-s)/h2)  )/h1/h2

    estequ <- function(beta) {

      expbeta <- as.vector(exp(beta*covariates  ))
      S0_XX <-  t(kerval_XX * Delta)* expbeta
      Zbar_XX <-  t(S0_XX) %*% covariates/  colSums(S0_XX)
      Zbar_XX[is.na(Zbar_XX)] <- 0
      res <- sum(delta*kerval*(covariates-Zbar_XX))
      res
    }

    estres <- nleqslv(0,fn=estequ)

    beta_est <- estres$x


    list(est=beta_est)

  }
    if(length(bd)==2){

      h1=bd[1]
      h2=bd[2]

    }else{

      test_idk <- lapply(split(sample(1:n,n),rep(1:2,n/2)),sort)
      hmax <- bd[1]
      hmin <- bd[2]
      testnn <- bd[3]
      hn <- exp(seq(log(hmin),log(hmax),length.out=testnn))

      hat_Mse=0
      for (s in times){
        hat_var= sbeta= bd_id = NULL
         for (n1 in 1:testnn ){
           hh1=hn[n1]
          for (n2 in 1:testnn ){
              hh2=hn[n2]
          beta1=estproc_Cox_test(data %>% filter(!(as.numeric(id) %in% test_idk$`1`)) ,
                                 n/2,s,hh1,hh2)$est
          beta2=estproc_Cox_test(data %>% filter(!(as.numeric(id) %in% test_idk$`2`)),
                                 n/2,s,hh1,hh2)$est
          sbeta=c(sbeta,estproc_Cox_test(data,
                                 n,s ,hh1,hh2)$est)
          hat_var=c(hat_var,sum((beta1-beta2)^2/4))
          bd_id=rbind(bd_id,c(hh1^2,hh1*hh2,hh2^2))
          }
         }
        if(p_z==1){
          C=lm(sbeta~bd_id)$coefficients[-1]
          hat_Mse=hat_Mse+( ( bd_id%*%C)^2 +hat_var )
        }else{
          C=lm(sbeta~bd_id)$coefficients[-1,]
          hat_Mse=hat_Mse+( rowSums(( bd_id%*%C)^2) +hat_var )
        }

      }

      min_Mse= which.min(hat_Mse)
      h1=sqrt(bd_id[min_Mse,1])
      h2=sqrt(bd_id[min_Mse,3])
    }


    if(length(scb)==2){
      alpha=scb[1]
      M=scb[2]
    }

    est.b  = se.b = tilde_S = NULL
    for(i in 1:nt ){
      s=times[i]

    kerval <- kerfun((X - s) / h1) * kerfun((obs_times-s)/h2) / h1 /h2
    Delta<-outer(X,X,"<=")
    kerval_XX<- kerfun((X-s)/h1) %*%t(kerfun((obs_times-s)/h2)  )/h1/h2


    estequ <- function(beta) {

      expbeta <- as.vector(exp( c(covariates %*% beta)  ))
      S0_XX <-  t(kerval_XX * Delta)* expbeta
      Zbar_XX <-  t(S0_XX) %*% covariates/  colSums(S0_XX)
      Zbar_XX[is.na(Zbar_XX)] <- 0
      res <- colSums(delta*kerval*(covariates-Zbar_XX))
      res
    }

    estres <- nleqslv(rep(0,p_z),fn=estequ)

    beta_est <- estres$x
    expbeta <- as.vector(exp(covariates %*% beta_est))
    S0_XX <-  t(kerval_XX * Delta)* expbeta
    Zbar_XX <-   t(S0_XX) %*% covariates/  colSums(S0_XX)
    Zbar_XX[is.na(Zbar_XX)] <- 0

    if(p_z==1){
      dZbar_XX=   t(S0_XX) %*% (covariates^2)*colSums(S0_XX)/colSums(S0_XX)^2-( t(S0_XX) %*% covariates)^2/colSums(S0_XX)^2
      dZbar_XX[is.na(dZbar_XX)] <- 0
      subid=B_indi_inner=NULL
      B_indi <- tibble(B_indi_inner =  delta*kerval*(covariates-Zbar_XX),subid= as.integer(id)) %>% group_by(subid) %>% reframe(B=t(colSums(B_indi_inner)))
      B<-sum((apply(B_indi$B,1,outf)))/n^2
      A_indi <-  - delta*kerval*dZbar_XX
      A <-  sum(A_indi)/n
    }else{
      dZbar_XX= ( t(S0_XX) %*% t(apply(covariates,1,outf))/colSums(S0_XX)-t(apply(Zbar_XX,1,outf)))
      dZbar_XX[is.na(dZbar_XX)] <- 0
      B_indi <- tibble(B_indi_inner =  delta*kerval*(covariates-Zbar_XX),subid= as.integer(id)) %>% group_by(subid) %>% reframe(B=t(colSums(B_indi_inner)))
      B<-matrix(rowSums((apply(B_indi$B,1,outf)))/n^2,ncol=p_z)
      A_indi <-  - delta*kerval*dZbar_XX
      A <-  matrix(colSums(A_indi)/n,ncol=p_z)
    }

    sigma= sqrt(diag(solve(A) %*% B %*% solve(A)))
    est.b = cbind(est.b, beta_est)
    se.b = cbind(se.b, sigma)

    if(length(scb)==2){
    if(p_z==1)sigma=matrix(sigma)
    xi=matrix(rexp(M*n)-1,n)
    res=  (t(as.matrix(B_indi$B))%*%xi)
    tilde_s=l%*%abs( solve(diag(sigma))%*%solve(A)%*%(res))/n
    tilde_S = cbind(tilde_S,t(tilde_s))
     }

    }

    c_alpha=ifelse(length(scb)==2,quantile(apply(tilde_S,1,max),alpha),0)

    list(est=beta_est,sigma=sigma,tilde_s=c_alpha )

  }
