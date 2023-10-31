#' @title Additive hazards model with sparse longitudinal covariates
#'
#' @description Regression analysis of additive hazards model with sparse longitudinal covariates. Three different weighting schemes are provided to impute the missing values.
#'
#' @param  data  An object of class tibble. The structure of the tibble must be: tibble(id = id, X = failure time, covariates = observation for covariates, obs_times = observation times, delta = censoring indicator).
#'
#' @param  n An object of class integer. The sample size.
#'
#' @param  tau An object of class numeric. The pre-specified time endpoint.
#'
#' @param  h An object of class vector. If use auto bandwidth selection, the structure of the vector must be: h = c(the maximum bandwidth, the minimum bandwidth, the number of bandwidth divided). If use fixed bandwidth, h is the chosen bandwidth.
#'
#' @param  method An object of class integer. If use weighted LVCF, method = 1. If use half kernel, method = 2. If use full kernel, method = 3.
#'
#' @references  Sun, Z. et al. (2022) <doi:10.1007/s10985-022-09548-6>
#'
#' @return a list with the following elements:
#' \item{est}{The estimation for the corresponding parameters.}
#' \item{se}{The estimation for the standard error of the estimated parameters.}
#'
#' @import  gaussquad
#' @importFrom dplyr %>% bind_rows group_by pull reframe filter
#' @importFrom tidyr expand_grid
#' @importFrom purrr array_branch
#' @importFrom tibble tibble
#' @importFrom stats  lm
#'
#' @examples
#'
#' library(gaussquad)
#' library(dplyr)
#' library(nleqslv)
#' library(MASS)
#' n=500
#' lqrule64 <- legendre.quadrature.rules(64)[[64]]
#' simdata <- function(alpha,beta ) {
#' cen=1
#' nstep=20
#' Sigmat_z <- exp(-abs(outer(1:nstep, 1:nstep, "-")) / nstep)
#' z  <-   c(mvrnorm(  1, c(1: nstep)/2, Sigmat_z  ))
#' left_time_points <- (0:(nstep - 1)) / nstep
#' z_fun <- stepfun(left_time_points, c(0,z  ))
#' lam_fun <- function(tt) {  alpha(tt)+beta*z_fun(tt)}
#' u <- runif(1)
#' fail_time <- nleqslv( 0 , function(ttt)
#'  legendre.quadrature(lam_fun,
#'                    lower = 0,
#'                    upper = ttt,
#'                    lqrule64) + log(u))$x
#'
#' X <- min(fail_time, cen)
#' obs=rpois(1,5)+1
#' tt= sort(runif(obs, min = 0, max = 1))
#' obs_times <- tt[which(tt<=cen)]
#' if (length(obs_times) == 0)
#'  obs_times <- cen
#'  covariates_obscov <-z_fun(obs_times)
#'  return( tibble(X = X,delta = fail_time < cen,
#'  covariates = covariates_obscov,obs_times = obs_times, censoring = cen  )  ) }
#' data <- replicate(n, simdata(alpha = function(tt)  tt, 1  ),
#'               simplify = FALSE ) %>% bind_rows(.id = "id")
#'
#' add.haz(data,n,1,n^(-0.5),3)
#' @export
#'
#'
add.haz <- function(data,n,tau, h, method) {
  lqrule64 <- legendre.quadrature.rules(64)[[64]]
  kerfun <- function(xx){
    pmax((1-xx^2)*0.75,0)
  }
  dist=subid=res=A_indi_inner=B_indi_inner=subjid=id=NULL
  f<-function(data,n,tau,h,method){
    X <- data$X/tau
    id <- data$id
    covariates <- data$covariates
    obs_times <- data$obs_times/tau
    delta <- data$delta

    p_z <- length(covariates)/length(X)
    lqrulepoints <- 0.5 * lqrule64$x + 0.5
    outf=function(x){x %o% x}
    outmf=function(x,y){x %o% y}
    dist_tt <- outer(lqrulepoints, obs_times, "-")
    dist_XX <- outer(X, obs_times, "-")


    if(method==1){
    j_tt<-dist_tt
    ismin<-function(a){
      a[which(a<0)]=100
      ifelse((abs(a)==min(abs(a)))&(a<99),1,0) }
    ismin_mat= function(a){
      if(dim(a)[1]==1){matrix(rep(1,dim(a)[2]),1)}else{apply(a, 2,ismin)}
    }
    j_tt=tibble(dist=t(dist_tt),subid= as.integer(id)) %>%
      group_by(subid)%>%reframe(res=ismin_mat(dist))  %>%
      pull(res)
    j_tt=t(j_tt)
    kerval_tt <- (kerfun(dist_tt / h) / h)*j_tt*(dist_tt>=0)

    j_XX<-dist_XX
    j_XX=tibble(dist=t(dist_XX),subid= as.integer(id))%>%
      group_by(subid)%>%reframe(res=ismin_mat(dist))%>%
      pull(res)
    j_XX=t(j_XX)
    kerval_XX <- kerfun(dist_XX / h) / h  * j_XX*(dist_XX >= 0)
    kerval <- kerfun((X - obs_times) / h) / h *diag(j_XX)*((X - obs_times) >= 0)

  }
    if(method==2){

      kerval_tt <- (kerfun(dist_tt / h) / h)*(dist_tt>=0)
      kerval <- kerfun((X - obs_times) / h) / h *((X - obs_times) >=0)
      kerval_XX <- kerfun(dist_XX / h) / h * (dist_XX >= 0)

    }
    if(method==3){

      kerval_tt <- (kerfun(dist_tt / h) / h)
      kerval <- kerfun((X - obs_times) / h) / h
      kerval_XX <- kerfun(dist_XX / h) / h

    }

    S0_tt <- outer(lqrulepoints, X, "<=") * kerval_tt
    Zbar_tt <- S0_tt %*% covariates / rowSums(S0_tt)
    Zbar_tt[is.na(Zbar_tt)] <- 0


    if(p_z==1){
      cov_Zbar_mat <-  -outer(as.vector(Zbar_tt),covariates ,"-")
      A_indi <-tibble(A_indi_inner=c(t(lqrule64$w)%*%(S0_tt*t(t(cov_Zbar_mat)*covariates))*0.5),id=as.integer(id) )   %>%
        group_by(id) %>% reframe(A = sum(A_indi_inner))
      A <- sum(A_indi$A)/n
    }else{
      cov_Zbar_mat <- expand_grid(cov=covariates,Zbar=Zbar_tt)
      cov_m_Zbar <- cov_Zbar_mat$cov - cov_Zbar_mat$Zbar
      A_indi <-
        tibble(
          A_indi_inner = t(mapply(
            outmf,
            array_branch(cov_m_Zbar, 1),
            array_branch(cov_Zbar_mat$cov, 1)
          )) * as.vector(S0_tt) * lqrule64$w,
          subjid = rep(as.integer(id), rep(64, length(id))),
          tid = rep(1:64, length(id))
        ) %>% group_by(subjid) %>% reframe(A = t(colSums(A_indi_inner) * 0.5))
      A <- matrix(colSums(A_indi$A)/n,ncol=p_z)
    }



    S0_XX <- outer(X, X, "<=") * kerval_XX
    Zbar_XX <- S0_XX %*% covariates / rowSums(S0_XX)
    Zbar_XX[is.na(Zbar_XX)] <- 0
    B_indi <- tibble(B_indi_inner = (covariates-Zbar_XX)* kerval *delta,subid= as.integer(id)) %>% group_by(subid) %>% reframe(B=t(colSums(B_indi_inner)))
    B <- colSums(B_indi$B)/n
    beta <- solve(A) %*% B

    betaextmat <- do.call(rbind,lapply(beta,function(bb) diag(x=bb,p_z,p_z)))
    A_indi_beta <- A_indi$A %*% betaextmat

    if(p_z==1){
      Sigma <-sum((B_indi$B-A_indi_beta)^2)/n^2
    }else{
      Sigma <- matrix(rowSums(apply(B_indi$B-A_indi_beta,1,outf))/n^2,nrow=p_z)
    }

    sigma=solve(A) %*% Sigma %*% solve(A)
    se=sqrt(diag(solve(A) %*% Sigma %*% solve(A)))
    list(est=beta,se=se)}

  if(length(h)==3){
    test_idk <- lapply(split(sample(1:n,n),rep(1:2,n/2)),sort)
    nn=h[3]
    hmin<-h[2]
    hmax<-h[1]
    hn <- exp(seq(log(hmin),log(hmax),length.out=nn))

    hat_var=hat_beta=NULL
    covariates <- data$covariates
    p_z <- length(covariates)/dim(data)[1]

    if(p_z > 1){
     for(i in 1:nn) {
        hh=hn[i]
        beta1=f( data %>% filter(!(as.numeric(id) %in% test_idk$`1`)),n/2,tau,h=hh,method)$est
        beta2=f( data %>% filter(!(as.numeric(id) %in% test_idk$`2`)),n/2,tau,h=hh,method)$est
        sbeta=f(data, n ,tau,h=hh,method)$est
        hat_var=c(hat_var,sum((beta1-beta2)^2/4))
        hat_beta=cbind(hat_beta,sbeta)
      }
      C=lm(t(hat_beta)~hn)$coefficients[-1,]
      hat_Mse= colSums((matrix(C)%*%hn )^2) +hat_var
    }else{
        for(i in 1:nn) {
        hh=hn[i]
        beta1=f( data %>% filter(!(as.numeric(id) %in% test_idk$`1`)),n/2,tau,h=hh,method)$est
        beta2=f( data %>% filter(!(as.numeric(id) %in% test_idk$`2`)),n/2,tau,h=hh,method)$est
        sbeta=f(data, n ,tau,h=hh,method)$est
        hat_var=c(hat_var,sum((beta1-beta2)^2/4))
        hat_beta=rbind(hat_beta,sbeta)
        }

      C=lm(hat_beta~hn)$coefficients[-1]
      hat_Mse= (hn*C)^2 +hat_var

    }

    h=hn[which.min(hat_Mse)]

  }


    f(data,n,tau,h,method)

  }
