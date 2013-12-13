require(mc2d)
require(VGAM)

#N <- 15
#K <- 3

# creates a random symetrix matrix
#M <- matrix(runif(K*K), nrow=K, ncol=K)
#M <- M + t(M)
#M <- M / max(M) * 0.9

# the multinomial parameters
#Theta <- rdirichlet(N, rep(1, K))

# creates a random binary symetric matrix
#Y <- matrix(runif(N*N), nrow=N, ncol=N)
#Y <- Y + t(Y)
#Y <- Y / max(Y)
#Y <- round(Y)
#diag(Y) <- 1

#tau.0 <- rdirichlet(N, rep(1,K))



gen.data <- function(N, m, theta){
  ys <- NULL
  zs <- apply(theta, 1, function(x) rmultinom(1, 1, x))
  ks <- apply(apply(zs, 1:2, as.logical), 2, which) # the vector of assignments
  for(i in seq(N)){
    zi <- zs[,i]
    for(j in seq(N)){
      if(i < j){
        zj <- zs[,j]
        ys <- c(ys, rbern(1, t(zi) %*% m %*% zj))
      }
    }
  }
  y <- matrix(0, nrow=N, ncol=N)
  y[lower.tri(y)] <- ys
  y <- t(y)
  y[lower.tri(y)] <- ys
  diag(y) <- 1
  return(list(y, ks))
}

get.true.empirical.m <- function(y, ks, K){
  uks <- unique(ks)
  empirical.m <- matrix(0, nrow=K, ncol=K)
  for(i in seq(K)){
    for(j in seq(K)){
      if(i <= j){
        if(i %in% ks && j %in% ks){
          i.ind <- ks == i
          j.ind <- ks == j
          propensity <- sum(y[i.ind, j.ind]) / (sum(i.ind) * sum(j.ind))
          empirical.m[i, j] <- propensity
          empirical.m[j, i] <- propensity
        } else {
          empirical.m[i, j] <- 0.00001
        }
      }
    }
  }
  return(empirical.m)
}

init.M <- function(K){
  m <- matrix(0, nrow=K, ncol=K)
  for(i in seq(K)){
    for(j in seq(K)){
      if(i <= j){
        m[i,j] <- runif(1)
        m[j,i] <- m[i,j]
      }
    }
  }
  return(m)
}


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

get.C <- function(i, q, tau, y, m, N, K){
  res.cumul = 0
  for(j in seq(N)){
    if( j != i ){
      for(l in seq(K)){
        res.cumul <- res.cumul + tau[j,l] * log(get.B(q, l, m, y[i,j]))
      }
    }   
  }
  return(res.cumul)
}


get.R <- function(i, q, theta){
  return(log(theta[i,q]))
}


get.tau.iq <- function(i, q, tau, y, m, theta, N, K){
  C.iq <- get.C(i, q, tau, y, m, N, K)
  R.iq <- get.R(i, q, theta)
  return(exp(C.iq + R.iq))
}

get.normtau.i <- function(i, tau, y, m, theta, N, K){
  unnorm.tau.i <- sapply(seq(K),
                         function(q) get.tau.iq(i, q, tau, y, m, theta, N, K))
  return(unnorm.tau.i / sum(unnorm.tau.i))
}

#get.norm.tau <- function(tau, y, m, theta, N, K){
#  new.tau <- NULL
#  for(i in seq(N)){
#    new.tau <- rbind(new.tau, get.normtau.i(i, tau, y, m, theta, N, K))
#  }
#  return(new.tau)
#}

get.norm.tau <- function(tau, y, m, theta, N, K){
  new.tau <- NULL
  for(i in seq(N)){
    tau[i,] <- get.normtau.i(i, tau, y, m, theta, N, K)
    new.tau <- rbind(new.tau, tau[i,])
  }
  return(new.tau)
}

opt.tau <- function(tau, y, m, theta, N, K, tau.tol=0.0001){
  old.tau <- tau
  new.tau <- 0 * tau
  err <- NULL
  current.err <- 1
  while(current.err > tau.tol){
    new.tau <- get.norm.tau(old.tau, y, m, theta, N, K)
    current.err <- sum( (new.tau - old.tau)^2 )
    err <- c(err, current.err)
    old.tau <- new.tau
  }
  return(list(new.tau, err))
}


# Maximization of M given posterior

get.H <- function(m, y, tau, N, K){
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

get.fit.H <- function(y, tau, N, K){
  f <- function(m.ls){
    m <- matrix(0, nrow=K, ncol=K)
    m[upper.tri(m, diag=T)] <- m.ls
    m <- t(m)
    m[upper.tri(m, diag=T)] <- m.ls
    return(get.H(m, y, tau, N, K))
  }
  return(f)
}


#opt.M <- function(init.par, y, tau, N, K, tol=0.0001){
#  fit.H <- get.fit.H(y, tau, N, K)
#  res <- optim(init.par, fit.H, control=list(fnscale=-1, maxit=5000,reltol=tol))
#  m.ls <- res$par
#  m <- matrix(0, nrow=K, ncol=K)
#  m[upper.tri(m, diag=T)] <- m.ls
#  m <- t(m)
#  m[upper.tri(m, diag=T)] <- m.ls
#  return(list(m, res$value))
#}

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

