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
  return(list(new.tau, err))
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


get.M.kl.update <- function(k,l,tau, y, N){
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
  return(num.kl / denom.kl)
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

      }
    }
  }

}

opt.M <- function(y, tau, N, K){
  new.M <- matrix(0, nrow=K, ncol=K)
  for(k in seq(K)){
    for(l in seq(K)){
      if(k <= l){
        new.M[k,l] <- get.M.kl.update(k,l,tau,y,N)
        new.M[k,l] <- rectify.M(new.M[k,l])
        new.M[l,k] <- new.M[k,l]
      }
    }
  }
  return(list(new.M, get.H(new.M, y, tau, N, K)))
}

# Maximization of J given the posterior

#get.J <- function(theta, tau, N, K){
#  cumul.res <- 0
#  for(i in seq(N)){
#    for(q in seq(K)){
#      cumul.res <- cumul.res + tau[i,q] * log(theta[i,q])
#    }
#  }
#  return(cumul.res)
#}
#
#simplexify <- function(unconstrained.thet){
#  expsum <- sum(exp(unconstrained.thet))
#  return(sapply(unconstrained.thet, function(x) exp(x) / expsum))
#}
#
#get.fit.J <- function(tau, N, K){
#  f <- function(unconstrained.thet){
#    unconstrained.thet.mat <- matrix(unconstrained.thet, nrow=N, byrow=T)
#    theta <- t(apply(unconstrained.thet.mat, 1, simplexify))
#    return(get.J(theta, tau, N, K))
#  }
#  return(f)
#}
#
#
#opt.J <- function(init.par, tau, N, K){
#  fit.J <- get.fit.J(tau, N, K)
#  res <- optim(init.par, fit.J, control=list(fnscale=-1, maxit=5000))
#  unconstrained.thet.mat <- matrix(res$par, nrow=N, byrow=T)
#  theta <- t(apply(unconstrained.thet.mat, 1, simplexify))
#  return(theta)
#}
#

tau.2.w <- function(tau){
  return(as.vector(t(tau)))
}

glmfitted.2.theta <- function(glm.fitted, K, N){
  return(glm.fitted[as.logical(rep(c(1,rep(0,K-1))), N),])
}


opt.J <- function(y.mod, x.mod, tau, N, K, glm.formula='y.mod~x.mod'){
  w <- tau.2.w(tau)
  glm.formula <- formula(glm.formula)
  fit.res <- vglm(formula=glm.formula, family=multinomial, weights=w)
  theta <- glmfitted.2.theta(fitted.values(fit.res), K, N)
  return(list(betas=fit.res@coefficients,  theta=theta, log.lik=logLik(fit.res)))
}


# Full Algo

em.algo <- function(tau.init, m.init, theta.init, y, x, N, K, M.tol=0.0001, tau.tol=0.0001, maxiter=10, err.bound=0.01, glm.formula='y.mod~x.mod', burnin=20, max.up=4){

  res <- list(tau=list(), m=list(), theta=list(), betas=list())
  tau <- tau.init
  m <- m.init
  theta <- theta.init

  # transforms inputs
  inputs.mods <- setup.glm(x, N, K)
  y.mod <- inputs.mods$ymod
  x.mod <- inputs.mods$xmod

  i <- 1
  old.lbound <- 1
  stop.iter = 0
  
  while(stop.iter < 3){

    cat("iterations: ", i, "\n")

    #if(i == 15){
    #  browser()
    #}
    
    # E step
    E.res <- opt.tau(tau.init, y, m, theta, N, K, tau.tol)
    tau <- E.res[[1]]
    err <- E.res[[2]]

    
    # M step
    #res.M <- opt.M(m[upper.tri(m, diag=T)], y, tau, N, K, M.tol)
    res.M <- opt.M(y, tau, N, K)
    m <- res.M[[1]]
    res.J <- opt.J(y.mod, x.mod, tau, N, K, glm.formula)
    theta <- res.J$theta
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
      stop.iter <- 3
    }
  }

  return(res)
  
}


