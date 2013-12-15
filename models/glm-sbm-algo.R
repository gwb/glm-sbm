require(mc2d)
require(VGAM)


get.B <- function(q, l, m, yij){
  mql <- m[q,l]
  return( mql^yij * (1-mql)^(1-yij) )
}

setup.glm <- function(x, N, K){
  y.mod <- matrix(rep(diag(K), N), nrow=K*N, byrow=T)
  x.mod <- matrix(rep(x, each=K), nrow=K*N)
  return(list(ymod=y.mod, xmod=x.mod))
}

# POSTERIOR P(Z_i | ...) # # # # # # # # # 


update.tau.step <- function(tau, y, m, theta, N, K){
  new.tau <- NULL
  for(i in seq(N)){
    unnorm.tau.i <- NULL
    for(q in seq(K)){  
      # Tau[i,q] update has two components
      # First component
      C.iq <- 0
      for( j in seq(N)[-i] ){
        for(l in seq(K)){
          C.iq <- C.iq + tau[j,l] * log(get.B(q,l,m, y[i,j]))
        }
      }
      # Second component
      R.iq <- log(theta[i,q])
      unnorm.tau.i <- c(unnorm.tau.i, exp(C.iq + R.iq)) # Normalizing
    }    
    # Immediately update tau for following iterations
    tau[i,] <- unnorm.tau.i / sum(unnorm.tau.i) 
  }
  return(tau)
}


update.tau <- function(tau, y, m, theta, N, K, tau.tol=0.0001, tau.maxiter=100){
  old.tau <- tau
  new.tau <- 0 * tau
  err <- NULL
  current.err <- 1
  i <- 0
  while(current.err > tau.tol){
    new.tau <- get.norm.tau(old.tau, y, m, theta, N, K)
    current.err <- max((new.tau - old.tau)^2 / old.tau)
    err <- c(err, current.err)
    old.tau <- new.tau
    i <- i + 1
    if(i >= tau.maxiter){
      break
    }
  }
  return(list(tau=new.tau, err=err))
}


# Maximization of M given posterior

eval.M.factor <- function(m, y, tau, N, K){
  cumul.res <- 0
  for(i in seq(N)){
    for(j in seq(N)){
      if(i < j){
        for(q in seq(K)){
          for(l in seq(K)){
            cumul.res <- cumul.res + tau[i,q] * tau[j,l] * log(get.B(q,l,m,y[i,j]))
          }
        }
      }
    }
  }
  return(cumul.res)
}

rectify.M <- function(m.val){
  if(m.val > 0.999){
    m.val <- 0.99
  }
  if(m.val < 0.001){
    m.val <- 0.001
  }
  return(m.val)
}

opt.M <- function(y, tau, N, K){
  new.M <- matrix(0, nrow=K, ncol=K)
  for(k in seq(K)){
    for(l in seq(K)){
      if(k <= l){
        # Performs the [k,l] update
        num.kl <- 0
        denom.kl <- 0
        for(i in seq(N)){
          for(j in seq(N)){
            if(i < j){
              num.kl <- num.kl + tau[i,k] * tau[j,l] * y[i,j]
              denom.kl <- denom.kl + tau[i,k] * tau[j,l]
            }
          }
        }
        new.M[k,l] <- rectify.M(num.kl / denom.kl)
        new.M[l,k] <- new.M[k,l]
      }
    }
  }
  return(list(m=new.M, log.lik=eval.M.factor(new.M, y, tau, N, K)))
}


# Maximization of J given the posterior

opt.J <- function(y.mod, x.mod, tau, N, K, glm.formula='y.mod~x.mod'){
  w <- as.vector(t(tau)) # transforms taus into regression weights
  glm.formula <- formula(glm.formula)
  fit.res <- vglm(formula=glm.formula, family=multinomial, weights=w)

  # extracts the thetas from the regression results
  theta <- fitted.values(fit.res)[as.logical(rep(c(1,rep(0,K-1))), N),]
  return(list(betas=fit.res@coefficients,  theta=theta, log.lik=logLik(fit.res)))
}


# Full Algo

em.algo <- function(tau.init, m.init, theta.init, y, x, N, K, M.tol=0.0001, tau.tol=0.0001, maxiter=10, err.bound=0.01, glm.formula='y.mod~x.mod', burnin=20, max.up=4, max.stop.iter=3){

  # initialize containers and variables
  res <- list(tau=list(), m=list(), theta=list(), betas=list())
  tau <- tau.init
  m <- m.init
  theta <- theta.init

  # transforms inputs
  inputs.mods <- setup.glm(x, N, K)
  y.mod <- inputs.mods$ymod
  x.mod <- inputs.mods$xmod

  # iterator control variables
  i <- 1
  old.lbound <- 1
  stop.iter = 0
  
  while(stop.iter < max.stop.iter){
    cat("iterations: ", i, "\n")

    # E step
    E.res <- update.tau(tau.init, y, m, theta, N, K, tau.tol)
    tau <- E.res$tau
    err <- E.res$err

    # M step
    res.M <- opt.M(y, tau, N, K)
    res.J <- opt.J(y.mod, x.mod, tau, N, K, glm.formula)
    m <- res.M$m
    theta <- res.J$theta

    # Is it really the lower bound, isn't there some constant I forgot?
    l.bound <- res.J$log.lik + res.M[[2]]
    
    
    # logging
    res$tau <- c(res$tau, list(tau))
    res$m <- c(res$m, list(m))
    res$theta <- c(res$theta, list(theta))
    res$betas <- c(res$betas, list(res.J$betas))


    cat("lower bound: ", l.bound, "\n")
    if(i > 1){
      tau.err <- sum((res$tau[[i]] - res$tau[[i-1]])^2) / length(res$tau[[i-1]])
      m.err <- sum((res$m[[i]] - res$m[[i-1]])^2) / length(res$tau[[i-1]])
      betas.err <- sum((res$betas[[i]] - res$betas[[i-1]])^2) / length(res$tau[[i-1]])
      cat("Tau error:   ", tau.err, "\n")
      cat("M error:     ", m.err, "\n")
      cat("Betas error: ", betas.err, "\n")

      if(abs(l.bound - old.lbound) < err.bound & i > burnin){
        stop.iter <- stop.iter + 1
      }
      if(l.bound < old.lbound & i > burnin){
        stop.iter <- stop.iter + 1
      }
      
      cat("Bound Err: ", abs(l.bound - old.lbound), "\n")
      
    }

    old.lbound <- l.bound
    i <- i + 1
    
    if(i >= maxiter){
      stop.iter <- max.stop.iter
    }
  }

  return(res)
  
}


