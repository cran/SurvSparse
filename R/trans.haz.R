#' @title Transformed hazards model with sparse longitudinal covariates
#'
#' @description Statistical inference on transformed hazards model with sparse longitudinal covariates. Kernel-weighted log-likelihood and sieve maximum log-likelihood estimation are combined to conduct statistical inference on the hazards function.
#'
#' @param  data  An object of class tibble. The structure of the tibble must be: tibble(id = id, X = failure time, covariates = observation for covariates, obs_times = observation times, delta = censoring indicator).
#'
#' @param  n An object of class integer. The sample size.
#'
#' @param  nknots An object of class integer. The number of knots for B-spline.
#'
#' @param  norder An object of class integer. The order of B-spline.
#'
#' @param  tau An object of class numeric. The maximum follow-up time.
#'
#' @param  s An object of class numeric. The parameter for Box-Cox transformation.
#'
#' @param  h An object of class vector. If use auto bandwidth selection, the structure of the vector must be: h = c(the maximum bandwidth, the minimum bandwidth, the number of bandwidth divided). If use fixed bandwidth, h is the chosen bandwidth.
#'
#' @references Sun, D. et al. (2023) <arXiv:2308.15549>
#'
#' @return a list with the following elements:
#' \item{est}{The estimation for the corresponding parameters.}
#' \item{se}{The estimation for the standard error of the estimated parameters.}
#'
#' @import gaussquad MASS tidyr foreach
#' @importFrom nloptr nloptr
#' @importFrom nleqslv nleqslv
#' @importFrom splines bs
#' @importFrom dplyr %>% bind_rows group_by pull reframe filter mutate ungroup row_number
#' @importFrom tibble tibble
#'
#' @examples
#'
#'
#' library(dplyr)
#' library(gaussquad)
#' library(nleqslv)
#' library(MASS)
#' n=200
#' lqrule64 <- legendre.quadrature.rules(64)[[64]]
#' simdata <- function(  beta ) {
#' cen=1
#' nstep=20
#' Sigmat_z <- exp(-abs(outer(1:nstep, 1:nstep, "-")) / nstep)
#' z  <-  2*(pnorm(c(mvrnorm(  1, rep(0,20), Sigmat_z  )))-0.5)
#' left_time_points <- (0:(nstep - 1)) / nstep
#' z_fun <- stepfun(left_time_points, c(0,z  ))
#' h_fun <- function(x) { beta  * z_fun(x)  }
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
#'  covariates = covariates_obscov,obs_times = obs_times, censoring = cen  )  )
#'   }
#' beta=1
#' data <- replicate( n, simdata( beta ), simplify = FALSE ) %>% bind_rows(.id = "id")
#' trans.haz(data,n,3,3,1,s=0,n^(-0.35))
#' @export
#'
#'
trans.haz <- function(data, n, nknots, norder,tau, s, h) {

  lqrule64 <- legendre.quadrature.rules(64)[[64]]
  outf=function(x){x %o% x}
  trans_fun <- function(x, s) {
    if (s == 0) {
      res <- exp(x)
    } else {
      res = ifelse(x < -1 / s, .Machine$double.eps, (s * x + 1) ^ (1 / s))
    }
    return(res)
  }

  trans_fun_d <- function(x, s) {
    if (s == 0) {
      res <- exp(x)
    } else {
      res = ifelse(x < -1 / s, .Machine$double.eps, (s * x + 1) ^ (1 / s - 1))
    }
    return(res)
  }

  trans_fun_d1o1 <- function(x, s) {
    if (s == 0) {
      res <- rep(1, length(x))
    } else {
      res = ifelse(x < -1 / s, .Machine$double.eps, (s * x + 1) ^ (-1))
    }
    return(res)
  }

  trans_fun_d12o1 <- function(x, s) {
    if (s == 0) {
      res <- exp(x)
    } else {
      res = ifelse(x < -1 / s, .Machine$double.eps, (s * x + 1) ^ (1 / s - 2))
    }
    return(res)
  }

  trans_fun_d12o2 <- function(x, s) {
    if (s == 0) {
      res <- rep(1, length(x))
    } else {
      res = ifelse(x < -1 / s, .Machine$double.eps, (s * x + 1) ^ (-2))
    }
    return(res)
  }

  kerfun <- function(xx){
    pmax((1-xx^2)*0.75,0)
  }

  logll_val <- function(para,data, n, nknots, norder, s, h, pl = 0) {

    gammap<- nknots+norder-1
    X <- data$X
    id <- data$id
    covariates <- data$covariates
    obs_times <- data$obs_times
    delta <- data$delta
    kerval <- kerfun((X - obs_times) / h) / h * (X > obs_times)

    knots <- (1:nknots) / (nknots + 1)

    bsmat <-
      bs(
        X,
        knots = knots,
         degree = norder,
        intercept = TRUE,
        Boundary.knots = c(0, 1)
      )

    loglik1_1 <- function(beta, gamma) {
      alphaX <- bsmat %*% gamma
      res <-
        sum(log(trans_fun(alphaX +  covariates %*% beta, s)) * delta * kerval)

      res
    }

    loglik1_1_d <- function(beta, gamma) {
      alphaBeta <- bsmat %*% gamma +covariates %*% beta
      temp1 <- trans_fun_d1o1(alphaBeta, s) * delta * kerval
      res <- as.vector(t(temp1) %*% cbind(covariates, bsmat))

      res
    }

    loglik2_inner_1 <- function(tt, beta, gamma) {
      dist <- outer(tt, obs_times, "-")
      kerval_tt <- kerfun(dist / h) / h * (dist > 0)
      alpha_tt <-
        bs(  tt,
             knots = knots,
             degree = norder,
             intercept = TRUE,
             Boundary.knots = c(0, 1)  ) %*% gamma %>% as.vector()
      res <-
        trans_fun(outer(alpha_tt,c(covariates %*% beta), "+"), s) * kerval_tt
      rowSums(outer(tt, X, "<") * res)
    }

    loglik2_inner_1_d <- function(tt, beta, gamma) {

      dist <- outer(tt, obs_times, "-")
      kerval_tt <- kerfun(dist / h) / h * (dist > 0)
      bsmat_tt <-
        bs(
          tt,
          knots = knots,
           degree = order,
          intercept = TRUE,
          Boundary.knots = c(0, 1)
        )

      alpha_tt <-
        bsmat_tt %*% gamma %>% as.vector()
      temp1 <-
        outer(tt, X, "<") * trans_fun_d(outer(alpha_tt,c(covariates %*% beta), "+"), s) * kerval_tt
      res <- cbind(temp1 %*% covariates, rowSums(temp1) * bsmat_tt)

      res
    }

    loglik_1 <- function(beta, gamma) {
      res <-
        loglik1_1(beta, gamma) -
        legendre.quadrature(
          function(tt)
            loglik2_inner_1(tt, beta, gamma),
          lower = 0,
          upper = 1,
          lqrule64
        )

      res / n

    }



    f <- function(xx) {
      -loglik_1(xx[1:ncol(data$covariates)], xx[-(1:ncol(data$covariates))])
    }
    f(para)


  }

  estproc_ori_dh <- function(data, n, nknots, norder, tau,s, h) {
    gammap <- nknots + norder-1
    X <- data$X/tau
    id <- data$id
    covariates <- matrix(data$covariates,length(X))
    obs_times <- data$obs_times
    delta <- data$delta
    kerval <- kerfun((X - obs_times) / h) / h * (X > obs_times)
    knots <- (1:nknots) / (nknots + 1)


    bsmat <-bs(  X,knots = knots,degree = norder, intercept=TRUE,Boundary.knots = c(0, 1)  )

    loglik1_1 <- function(beta, gamma) {
      alphaX <- bsmat %*% gamma
      res <-  sum(log(trans_fun(alphaX +  covariates %*% beta, s)) * delta * kerval)

      res
    }

    loglik1_1_d <- function(beta, gamma) {
      alphaBeta <- bsmat %*% gamma +   covariates %*% beta
      temp1 <- trans_fun_d1o1(alphaBeta, s) * delta * kerval
      as.vector(t(temp1) %*% cbind(covariates, bsmat))
    }

    loglik2_inner_1 <- function(tt, beta, gamma) {

      dist <- outer(tt, obs_times, "-")
      kerval_tt <- kerfun(dist / h) / h * (dist > 0)
      alpha_tt <-  bs(tt, knots = knots,degree = norder, intercept=TRUE,  Boundary.knots = c(0, 1)  ) %*% gamma %>% as.vector()
      res <- c(trans_fun(outer(alpha_tt,  covariates %*% beta, "+"), s)) * kerval_tt
      rowSums(outer(tt, X, "<") * res)
    }

    loglik2_inner_1_d <- function(tt, beta, gamma) {
      dist <- outer(tt, obs_times, "-")
      kerval_tt <- kerfun(dist / h) / h * (dist > 0)
      bsmat_tt <- bs(tt, knots = knots,degree = norder,intercept=TRUE, Boundary.knots = c(0, 1) )
      alpha_tt <-  bsmat_tt %*% gamma %>% as.vector()
      temp1 <- outer(tt, X, "<") * c(trans_fun_d(outer(alpha_tt,   covariates %*% beta, "+"), s)) * kerval_tt
      cbind(temp1 %*% covariates, rowSums(temp1) * bsmat_tt)
    }

    loglik_1 <- function(beta, gamma) {
      res <- loglik1_1(beta, gamma) -
        legendre.quadrature(
          function(tt)
            loglik2_inner_1(tt, beta, gamma),
          lower = 0,
          upper = 1,
          lqrule64
        )

      res / n

    }

    loglik_1_d <- function(beta, gamma) {
      res1 <-loglik1_1_d(beta, gamma)
      res2 <- loglik2_inner_1_d(0.5 * lqrule64$x + 0.5, beta, gamma)
      res2 <- as.vector(0.5 * colSums(lqrule64$w * res2))
      (res1 - res2) / n

    }

    f <- function(xx) {
      -loglik_1(xx[1:nc], xx[-(1:nc)])
    }

    estres <- nloptr(
      x0 = rep(0, nc + gammap),
      eval_f = f,
      eval_grad_f = function(xx)
        - loglik_1_d(xx[1:nc], xx[-(1:nc)]),
      opts = list(
        "algorithm" = "NLOPT_LD_LBFGS",
        "maxeval" = 10000,
        "print_level" = 0
      )
    )

    return(estres$solution)
  }

  if(length(h)==3){
    K=h[3]
    test_idk <- lapply(split(sample(1:n,n),rep(1:K,n/K)),sort)
    lengthh <- h[3]
    hmin <- h[2]
    hmax <- h[1]
    hn <- exp(seq(log(hmin),log(hmax),length.out=lengthh))
    test_idx_one=hh=bd=NULL
    res <- foreach(hh=hn) %do% {
      foreach(test_idx_one = test_idk) %do% {
        foldkpar=estproc_ori_dh(data %>%
                                  filter(!(as.numeric(id) %in% test_idx_one)),n*(1-1/K), 3, 3, s, hh, pl = 0)
        testing <- data %>% filter(as.numeric(id) %in% test_idx_one)
        test_logll <- logll_val(foldkpar,testing,n/K, 3, 3, s, hh, pl = 0)
        tibble(test_logll,bd=hh)
      } %>% bind_rows(.id="fold")
    }

    h=res %>% bind_rows() %>% group_by(bd) %>% reframe(cvloss=mean(test_logll))}

  gammap<- nknots+norder+1
  X <- data$X/tau
  id <- data$id
  covariates <- matrix(data$covariates,length(X))
  obs_times <- data$obs_times
  delta <- data$delta
  kerval <- kerfun((X - obs_times) / h) / h * (X > obs_times)
  knots <- (1:nknots) / (nknots + 1)
  nc=length(data$covariates)/length(X)
  bsmat <-bs( X, knots = knots,degree = norder,  intercept = TRUE,   Boundary.knots = c(0, 1) )

  loglik1_1 <- function(beta, gamma) {
    alphaX <- bsmat %*% gamma
    sum(log(trans_fun(alphaX +  covariates %*% beta, s)) * delta * kerval)
  }

  loglik1_1_d <- function(beta, gamma) {
    alphaBeta <- bsmat %*% gamma +covariates %*% beta
    temp1 <- trans_fun_d1o1(alphaBeta, s) * delta * kerval
    as.vector(t(temp1) %*% cbind(covariates, bsmat))
  }

  loglik2_inner_1 <- function(tt, beta, gamma) {

    dist <- outer(tt, obs_times, "-")
    kerval_tt <- kerfun(dist / h) / h * (dist > 0)
    alpha_tt <- bs(  tt,  knots = knots, degree = norder, intercept = TRUE,Boundary.knots = c(0, 1)  ) %*% gamma %>% as.vector()
    res <-trans_fun(outer(alpha_tt,c(covariates %*% beta), "+"), s)
    res <- ifelse(kerval_tt==0,0,res)* kerval_tt
    res <- ifelse(outer(tt, X, "<"),res,0)
    rowSums(res)
  }

  loglik2_inner_1_d <- function(tt, beta, gamma) {

    dist <- outer(tt, obs_times, "-")
    kerval_tt <- kerfun(dist / h) / h * (dist > 0)
    bsmat_tt <- bs( tt,knots = knots,  degree = norder,intercept = TRUE, Boundary.knots = c(0, 1)   )
    alpha_tt <- bsmat_tt %*% gamma %>% as.vector()
    temp1 <- trans_fun_d(outer(alpha_tt,c(covariates %*% beta), "+"), s)
    temp1 <- ifelse(outer(tt, X, "<"),temp1,0)
    temp1 <- ifelse(kerval_tt==0,0,temp1) * kerval_tt
    cbind(temp1 %*% covariates, rowSums(temp1) * bsmat_tt)
  }

  loglik_1 <- function(beta, gamma) {
    res <-
      loglik1_1(beta, gamma) -
      legendre.quadrature(
        function(tt)
          loglik2_inner_1(tt, beta, gamma),
        lower = 0,
        upper = 1,
        lqrule64
      )

    res / n

  }

  loglik_1_d <- function(beta, gamma) {
    res1 <-
      loglik1_1_d(beta, gamma)
    res2 <- loglik2_inner_1_d(0.5 * lqrule64$x + 0.5, beta, gamma)
    res2 <- as.vector(0.5 * colSums(lqrule64$w * res2))
    (res1 - res2) / n
  }


  f <- function(xx) {
    -loglik_1(xx[1:nc], xx[-(1:nc)])
  }


  A <- function(beta, gamma) {

    dist_XX <- outer(X, obs_times, "-")
    kerval_XX <- kerfun(dist_XX / h) / h * (dist_XX > 0)
    bsmat_XX <- bs(  X,   knots = knots,degree = norder, intercept = TRUE, Boundary.knots = c(0, 1)  )
    alpha_XX <- bsmat_XX %*% gamma %>% as.vector()
    inner <- outer(alpha_XX, covariates %*% beta, "+")
    S0 <- outer(X, X, "<=") * matrix(trans_fun_d12o1(inner, s),length(X)) * kerval_XX
    S1 <- S0 %*% covariates/ rowSums(S0)
    outf=function(x){x %o% x}
    if(nc==1){
      S2 <- S0 %*% apply(covariates,1,outf) / rowSums(S0)
      outerprod <-  (S2 - (apply(S1,1,outf))) * c(trans_fun_d1o1(alpha_XX +  covariates %*% beta, s) ) ^ 2
    }else{
      S2 <- S0 %*% ( t(apply(covariates,1,outf))) / rowSums(S0)
      outerprod <-  (S2 - t(apply(S1,1,outf))) * c(trans_fun_d1o1(alpha_XX +  covariates %*% beta, s) ) ^ 2

    }
    outerprod[is.na(outerprod)] <-0
    matrix(colSums( outerprod *  kerval * delta) / n,nc)
  }

  B <- function(beta, gamma) {

    dist_XX <- outer(X, obs_times, "-")
    kerval_XX <- kerfun(dist_XX / h) / h * (dist_XX > 0)
    bsmat_XX <- bs(X, knots = knots,norder,   intercept = TRUE,  Boundary.knots = c(0, 1)    )
    alpha_XX <- bsmat_XX %*% gamma %>% as.vector()
    inner <- outer(alpha_XX,   covariates %*% beta, "+")
    S0 <- outer(X, X, "<=") * matrix(trans_fun_d12o1(inner, s),length(X)) * kerval_XX
    S1 <- S0 %*% covariates / rowSums(S0)
    outerproducinner1 <- (S1 - covariates) * c(trans_fun_d1o1(alpha_XX +  covariates %*% beta, s))
    outerproducinner1[is.na(outerproducinner1)] <-0
    bb1=outerproducinner1*  kerval * delta * ( (X-obs_times)>0 )
    a=coe=cla=NULL
    bb1=tibble(a=bb1 ,id=as.numeric(id)) %>% group_by(id) %>%  reframe( cla=colSums( a ))
    nc=nc
    bb1 <- bb1 %>% group_by(id) %>% mutate(coe = row_number()) %>% pivot_wider(values_from=`cla`,names_from=coe) %>% ungroup %>% dplyr::select(-id)%>% as.matrix()


    b_inner <- function(tt, beta, gamma) {
      dist <- outer(tt, obs_times, "-")
      kerval_tt <- kerfun(dist / h) / h * (dist > 0)
      bsmat_tt <-bs(tt, knots = knots,degree = norder,   intercept = TRUE,    Boundary.knots = c(0, 1))

      alpha_tt <- bsmat_tt %*% gamma %>% as.vector()
      inner <- outer(alpha_tt,   covariates %*% beta  , "+")
      S0_tt <- outer(tt, X, "<=") * matrix(trans_fun_d12o1(inner, s),length(tt)) * kerval_tt
      S1_tt <- S0_tt %*% covariates / rowSums(S0_tt)
      res=NULL

      for( jj in 1:nc){
        temp1  <-
          outer(tt, X, "<=") * matrix(trans_fun_d(outer(alpha_tt, covariates %*% beta,"+"), s),length(tt))*
          kerval_tt     *outer(as.vector(S1_tt[,jj]),covariates[,jj],"-")
        temp1[is.na(temp1)] <-0
        cla=NULL
          r1=tibble(a=t(temp1) ,id=as.numeric(id)) %>% group_by(id) %>%  reframe( cla=colMeans(a))
          r1 <- r1 %>% group_by(id) %>% mutate(coe = row_number()) %>% pivot_wider(values_from=`cla`,names_from=coe) %>% ungroup %>% dplyr::select(-id)%>% as.matrix()
        res=cbind(res,r1)
      }
      res
    }
    binner <- b_inner (0.5 * lqrule64$x + 0.5, beta, gamma)
    bb2=NULL
    for( jj in 1:nc){
      b1 <- as.vector(0.5 * rowSums(lqrule64$w * binner[,(1+(nc-1)*64):(64*nc)]))
      bb2=cbind(bb2,b1)
    }
    bb=bb1-bb2
    if(nc==1){
       sum( apply(bb,1,outf)) / n
    }else{
      matrix(colSums(t(apply(bb,1,outf))) / n,nc)
    }
  }


  if(s==0) {
    estres <- nloptr(
      x0 = rep(0, nc + gammap),
      eval_f = f,
      eval_grad_f = function(xx)
        - loglik_1_d(xx[1:nc], xx[-(1:nc)]),
      opts = list(
        "algorithm" = "NLOPT_LD_SLSQP",
        "xtol_rel"=1.0e-6,
        "maxeval" = 10000,
        "print_level" = 0
      )
    )
  } else {
    ineqmat <- cbind(covariates,bsmat)
    ineqmat <- ineqmat[(X>obs_times)& (abs(X-obs_times) <= h),]
    estres <- nloptr(
      x0 = rep(0, nc + gammap),
      eval_f = f,
      eval_grad_f = function(xx)
        - loglik_1_d(xx[1:nc], xx[-(1:nc)]),
      eval_g_ineq = function(xx) -ineqmat %*%xx-1/s,
      eval_jac_g_ineq = function(xx) -ineqmat,
      opts = list(
        "algorithm" = "NLOPT_LD_SLSQP",
        "xtol_rel"=1.0e-6,
        "maxeval" = 10000,
        "print_level" = 0
      )
    )
  }

  A_est <- A(estres$solution[1:nc], estres$solution[-(1:nc)])
  B_est <- B(estres$solution[1:nc], estres$solution[-(1:nc)])

  se=sqrt(diag(solve(A_est) %*% B_est %*% solve(A_est)/n))
  list(est=estres$solution,se=se)

}
