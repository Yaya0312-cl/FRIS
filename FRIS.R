########################
#### load libraries ####
########################
library(MASS)
library(fda)
library(fpca)
library(MAVE) 

########################
#### initial values ####
########################

#### get spline basis 

get.spline.t <- function(time, rangeval.t, norder.t, nk.t){
  breaks.t <- seq(rangeval.t[1], rangeval.t[2], length=nk.t)
  norder <- norder.t
  bt <- list()
  for (i in 1:sample.len) {bt[[i]] <- bsplineS(time[[i]], breaks = breaks.t, norder = norder.t, nderiv = 0)}
  nbisis.t <- norder.t + nk.t - 2
  return(list(bt = bt, nbisis.t = nbisis.t))
}

get.spline.u <- function(u, rangeval.u, norder.u, nk.u){
  breaks.u <- c(rangeval.u[1], quantile(u, seq(0, 1, length=nk.u)[-c(1,nk.u)]), rangeval.u[2])
  bu.0<-bsplineS(u, breaks = breaks.u, norder = norder.u, nderiv = 0)
  bu.1<-bsplineS(u, breaks = breaks.u, norder = norder.u, nderiv = 1)
  bu.2<-bsplineS(u, breaks = breaks.u, norder = norder.u, nderiv = 2)
  nbisis.u <- norder.u + nk.u - 2
  return(list(bu.0 = bu.0, bu.1 = bu.1, bu.2 = bu.2, nbisis.u = nbisis.u))
}

#### get initial value of mean function

get.mu.ini <- function(time, x, y, beta.ini, gamma.ini, max.iter, 
                       threshold, spline.par) {
  Beta <- beta.ini
  Gamma <- gamma.ini

  iter <- 1
  err <- 1
  get.bt <- get.spline.t(time = time, rangeval.t = spline.par$rangeval.t, 
                         norder.t = spline.par$norder.t, nk.t = spline.par$nk.t)
  bt <- get.bt$bt
  nbisis.t <- get.bt$nbisis.t
  
  ## least square interate
  while (err > threshold && iter < max.iter) {
    beta<-Beta; gamma<-Gamma
    get.bu <- get.spline.u(u = x %*% beta, rangeval.u = spline.par$rangeval.u,
                           norder.u = spline.par$norder.u, nk.u = spline.par$nk.u)
    bu.0 <- get.bu$bu.0
    bu.1 <- get.bu$bu.1
    bu.2 <- get.bu$bu.2
    nbisis.u <- get.bu$nbisis.u
    nbisis.mu <- nbisis.u * nbisis.t
    #### updata beta
    Bn <- list(); Bn1 <- list(); Bn2 <- list()
    for (i in 1:sample.len) { 
      Bn[[i]] <- kronecker(t(bt[[i]]), bu.0[i, ]) 
      Bn1[[i]] <- kronecker(t(bt[[i]]), bu.1[i, ]) 
      Bn2[[i]] <- kronecker(t(bt[[i]]), bu.2[i, ]) } 
    
    ## get mean curve
    mean.0 <- list(); mean.1 <- list(); mean.2 <- list()
    for (i in 1:sample.len) { 
      mean.0[[i]] <- t(gamma) %*% Bn[[i]] 
      mean.1[[i]] <- t(gamma) %*% Bn1[[i]]
      mean.2[[i]] <- t(gamma) %*% Bn2[[i]] }
    
    beta.der<-matrix(NA, p, sample.len)
    for (i in 1:sample.len) { 
      beta.der[ ,i]<-x[i,] %*% (y[[i]] - mean.0[[i]]) %*% t(mean.1[[i]])
      }
    beta.hessian<-diag(0, p, p)
    for (i in 1:sample.len) {
      beta.hessian <- beta.hessian - 
        x[i,] %*% mean.1[[i]] %*% t(mean.1[[i]]) %*% t(x[i,]) + 
        x[i,] %*% (y[[i]] - mean.0[[i]]) %*% t(mean.2[[i]]) %*% t(x[i,])
    }
    
    Beta <- beta - solve(beta.hessian) %*% rowSums(beta.der)
    Beta <- Beta/sqrt(sum(Beta^2))
    if(Beta[1] < 0) { Beta = -Beta }
    
    #### updata gamma
    gamma.der<-matrix(NA, nbisis.mu, sample.len)
    for (i in 1:sample.len) {
      gamma.der[,i] <- Bn[[i]] %*% t(y[[i]] - mean.0[[i]])
    }
    gamma.hessian<-diag(0, nbisis.mu, nbisis.mu)
    for (i in 1:sample.len) {
      gamma.hessian <- gamma.hessian - Bn[[i]] %*% t(Bn[[i]])
    }
    
    Gamma <- gamma - solve(gamma.hessian) %*% rowSums(gamma.der)
    
    err <- max(abs(beta-Beta))
    iter <- iter + 1
  }
  
  return(list(beta.ini = beta, gamma.ini = gamma, mean = mean.0))
}

#### get initial value of eigenfunctions

get.eigenfun.ini <- function(K, spline.par){
  rangeval.t1 <- spline.par$rangeval.t1
  nk.t1 <- spline.par$nk.t1
  norder.t1 <- spline.par$norder.t1
  #### get residual
  mean <- mu.ini$mean
  residual <- list()
  for (i in 1:sample.len) {
    residual[[i]] <- y[[i]] - mean[[i]]
  }
  
  #### FPCA of residual
  z <- cbind(rep(1:sample.len, time.len), unlist(residual), unlist(time))
  k <- K
  fpc <- fpca.mle(z, M.set = k, r.set = k)
  eigenfun <- fpc$eigenfunctions
  grids <- fpc$grid
  eigenval <- fpc$eigenvalues
  #### generate bn1
  breaks <- seq(rangeval.t1[1], rangeval.t1[2], length = nk.t1)
  get.bt1 <- bsplineS(grids, breaks, norder.t1)
  nbisis.t1 <- norder.t1 + nk.t1 - 2
  
  #### get eta
  eta <- matrix(NA, nbisis.t1, K)
  for (k in 1:K) {
    eta[ ,k] <- lm(eigenfun[k,] ~ get.bt1-1)$coef
  }
  eta <- qr.Q(qr(eta))

  return(list(eta.ini = eta, residual = residual, eigenfun = eigenfun, 
              eigenval = eigenval, grids = grids, K = K))
}

#### get initial value of score variance functions

get.scorevar.ini <- function(time, x, n.cl, spline.par){
  rangeval.u1 <- spline.par$rangeval.u1
  nk.u1 <- spline.par$nk.u1
  norder.u1 <- spline.par$norder.u1
  K <- eigenfun.ini$K
  eigenfun <- eigenfun.ini$eigenfun
  grids <- eigenfun.ini$grids
  ## cluster
  kmean <- kmeans(x, n.cl)
  coef<-NULL
  for (factor in 1:n.cl) {
    index <- which(kmean$cluster == factor)
    coef0 <- matrix(0, length(index), K)
    n<-1
    for (i in index) {
      x.cl <- NULL
      for (k in 1:K) {
        x.cl <- cbind(x.cl, eigenfun[k,ceiling(length(grids) * time[[i]])])
      }
      coef0[n,] <- lm(t(residual[[i]]) ~ x.cl - 1)$coef
      n<-n+1
    }
    coef0<-apply(coef0, 2, sd)
    coef<-rbind(coef,t(matrix(rep(coef0,length(index)), K, length(index))))
  }
  index<-NULL
  for (factor in 1:n.cl) {
    index<-c(index,which(kmean$cluster == factor))
  }
  
  ## updata alpha
  alpha <- matrix(NA, p, K)
  for (k in 1:K) {
    alpha[ ,k] <- coef(mave(coef[ ,k] ~ x[index, ]),1)
  }

  ## updat lambda
  uu <- x[index, ] %*% alpha
  breaks.u1 <- c(rangeval.u1[1], quantile(uu, seq(0, 1, length = nk.u1)[-c(1, nk.u1)]), 
                 rangeval.u1[2])
  nbisis.u1 <- norder.u1 + nk.u1 - 2
  lambda <- matrix(NA, nbisis.u1, K)
  zeta <- matrix(NA, sample.len, K)
  for (k in 1:K) {
    get.bu1 <- bsplineS(uu[ ,k], breaks.u1, norder.u1)
    lambda[ ,k] <- lm(coef[ ,k] ~ get.bu1 - 1)$coef
    zeta[ ,k] <- (get.bu1 %*% lambda[ ,k])^2
  }
  
  return(list(alpha.ini = alpha, lambda.ini = lambda, zeta.ini = zeta))
}

#### get initial value of sigma

get.sigma.ini <- function(time, grids, eigenfun.ini, zeta.ini, residual){
  x.cov <- NULL
  for (i in 1:sample.len) {
    x.cov.k <- 0
    for (k in 1:K) {
      x.cov.k <- x.cov.k + eigenfun.ini[k,ceiling(length(grids) * time[[i]])] %*% 
        t(eigenfun.ini[k,ceiling(length(grids) * time[[i]])]) * zeta.ini[i,k]
    }
    x.cov <- c(x.cov, mean(diag(x.cov.k - 
                                 t(residual[[i]]) %*% residual[[i]])))
  }
  return(list(sigma.ini = mean(x.cov)))
}

###################
#### FRIS code ####
###################

FRIS <- function(Beta.ini, Gamma.ini, Eta.ini, Alpha.ini, Lambda.ini, 
                 Sigma.ini, Zeta.ini, max.iter, threshold, spline.par, 
                 tuning.par, K){
  Beta <- Beta.ini
  Gamma <- Gamma.ini
  Eta <- Eta.ini
  Alpha <- Alpha.ini
  Lambda <- Lambda.ini
  Sigma <- Sigma.ini
  Zeta <- Zeta.ini
  C <- matrix(0,sample.len,4)
  rangeval.t = spline.par$rangeval.t; norder.t = spline.par$norder.t; nk.t = spline.par$nk.t
  rangeval.u = spline.par$rangeval.u; norder.u = spline.par$norder.u; nk.u = spline.par$nk.u
  rangeval.t1 = spline.par$rangeval.t1; norder.t1 = spline.par$norder.t1; nk.t1 = spline.par$nk.t1
  rangeval.u1 = spline.par$rangeval.u1; norder.u1 = spline.par$norder.u1; nk.u1 = spline.par$nk.u1
  err <- 1
  iter <- 1
  get.bt <- get.spline.t(time, rangeval.t, norder.t, nk.t)
  bt <- get.bt$bt
  nbisis.t <- get.bt$nbisis.t
  v <- tuning.par$v
  h <- tuning.par$h
  scad.l <- tuning.par$scad.l
  scad.a <- tuning.par$scad.a
  step <- tuning.par$step
  while (iter < max.iter && err[iter] > threshold ) {
    beta <- Beta; gamma <- Gamma; eta <- Eta;
    alpha <- Alpha; lambda <- Lambda; sigma <- Sigma;
    zeta <- Zeta; c <- C
    u <- x %*% beta
    uu <- x %*% alpha
    
    ## get spline basis
    get.bu <- get.spline.u(u = x %*% beta, rangeval.u = spline.par$rangeval.u,
                           norder.u = spline.par$norder.u, nk.u = spline.par$nk.u)
    bu.0 <- get.bu$bu.0
    bu.1 <- get.bu$bu.1
    bu.2 <- get.bu$bu.2
    nbisis.u <- get.bu$nbisis.u

    bu1.0 <- list()
    bu1.1 <- list()
    for (k in 1:K) {
      get.bu1 <- get.spline.u(uu[,k], rangeval.u1, norder.u1, nk.u1)
      bu1.0[[k]] <- get.bu1$bu.0
      bu1.1[[k]] <- get.bu1$bu.1
    }
    nbasis.u1 <- get.bu1$nbisis.u
    
    nbisis.mu <- nbisis.u * nbisis.t
    
    Bn <- list(); Bn1 <- list(); Bn2 <- list()
    for (i in 1:sample.len) { 
      Bn[[i]] <- kronecker(t(bt[[i]]), bu.0[i, ]) 
      Bn1[[i]] <- kronecker(t(bt[[i]]), bu.1[i, ]) 
      Bn2[[i]] <- kronecker(t(bt[[i]]), bu.2[i, ]) } 
    
    ## get mean curve
    mean.0 <- list(); mean.1 <- list(); mean.2 <- list()
    for (i in 1:sample.len) { 
      mean.0[[i]] <- t(gamma) %*% Bn[[i]] 
      mean.1[[i]] <- t(gamma) %*% Bn1[[i]]
      mean.2[[i]] <- t(gamma) %*% Bn2[[i]] }
    
    ## get covarince function of each individual
    sigma.ni <- list()
    for (i in 1:sample.len) {
      cov <- 0
      for (k in 1:K) {
        cov <- cov + (zeta[i,k])^2*eta[,k]%*%t(eta[,k])
      }
      sigma.ni[[i]] <- bt[[i]] %*% cov %*% t(bt[[i]]) + diag(sigma, time.len[i], time.len[i]) 
    }
    
    ## update beta
    beta.der<-matrix(NA, p, sample.len)
    for (i in 1:sample.len) { 
      beta.der[ ,i]<-x[i,] %*% (y[[i]] - mean.0[[i]] + t(beta) %*% x[i,] %*% mean.1[[i]]) %*%
        solve(sigma.ni[[i]]) %*% t(mean.1[[i]])
    }
    beta.hessian<-diag(0, p, p)
    for (i in 1:sample.len) {
      beta.hessian <- beta.hessian + 
        x[i, ] %*% mean.1[[i]] %*% solve(sigma.ni[[i]]) %*% t(mean.1[[i]]) %*% t(x[i,])
    }
    Beta <- solve(beta.hessian) %*% rowSums(beta.der)
    Beta <- Beta/sqrt(sum(Beta^2))
    if(Beta[1] < 0) { Beta = -Beta }
    
    
    #####update gamma
    gamma.der<-matrix(NA, nbisis.mu, sample.len)
    for (i in 1:sample.len) {
      gamma.der[,i] <- Bn[[i]] %*% solve(sigma.ni[[i]]) %*% y[[i]]
    }
    gamma.hessian<-diag(0, nbisis.mu, nbisis.mu)
    for (i in 1:sample.len) {
      gamma.hessian <- gamma.hessian + Bn[[i]] %*%solve(sigma.ni[[i]]) %*% t(Bn[[i]])
    }
    
    Gamma <- solve(gamma.hessian) %*% rowSums(gamma.der)
    
    ## update lambda
    lambda.der <- matrix(0, nbasis.u1, K)
    for (k in 1:K) {
      for (i in 1:sample.len) {
        lambda.der[ ,k] <- lambda.der[ ,k] + (zeta[i,k] + c[i,k]/v) %*% bu1.0[[k]][i, ]
      }
    }
    lambda.hessian <- array(0, c(nbasis.u1, nbasis.u1, K))
    for (k in 1:K) {
      for (i in 1:sample.len) {
        lambda.hessian[,,k] <- lambda.hessian[,,k] + bu1.0[[k]][i, ] %*% t(bu1.0[[k]][i, ])
      }
    }
    Lambda <- NULL
    for (k in 1:K) {
      Lambda <- cbind(Lambda, solve(lambda.hessian[,,k]) %*% lambda.der[,k])
    }
    
    ## update alpha
    alpha.der <- matrix(0, p, K)
    for (k in 1:K) {
      for (i in 1:sample.len) {
        alpha.der[ ,k] <- alpha.der[ ,k] + t(lambda[ ,k] %*% bu1.1[[k]][i, ] %*% x[i,]) %*% 
          (zeta[i,k] + c[i,k]/v - lambda[ ,k] %*% bu1.0[[k]][i, ] + 
             lambda[ ,k] %*% bu1.1[[k]][i, ] * (x[i, ] %*% alpha[ ,k]))
      }
    }
    alpha.hessian <- array(0, c(p, p, K))
    for (k in 1:K) {
      for (i in 1:sample.len) {
        alpha.hessian[,,k] <- alpha.hessian[,,k] + x[i, ]%*%(lambda[ ,k]%*%bu1.1[[k]][i, ])^2%*%t(x[i, ])
      }
    }
    Alpha <- matrix(NA, p, K)
    for (k in 1:K) {
      Alpha[ ,k] <- solve(alpha.hessian[,,k]) %*% alpha.der[ ,k]
      Alpha[ ,k] <- Alpha[,k]/sqrt(sum(Alpha[ ,k]^2))
      if(Alpha[1,k] < 0) { Alpha[,k] = -Alpha[,k] }
    }
    
    ## update zeta
    sigma.ni.1 <- list()
    for (i in 1:sample.len) {
      sigma.ni.1[[i]] <- solve(sigma.ni[[i]])/2 - solve(sigma.ni[[i]])%*%
        t(y[[i]]-mean.0[[i]]) %*% (y[[i]]-mean.0[[i]]) %*% solve(sigma.ni[[i]])/2
    }
    
    W <- matrix(NA, sample.len, K)
    for (k in 1:K) {
      for (i in 1:sample.len) {
        W[i,k] <- lambda[,k] %*% bu1.0[[k]][i, ]
      }
    }
    H.der<-matrix(NA, sample.len, K)
    for (k in 1:K) {
      for (i in 1:sample.len) {
        H.der[i,k]<-as.vector(sigma.ni.1[[i]]) %*% as.vector(
          2 * bt[[i]] %*% eta[ ,k] %*% t(eta[ ,k]) %*% t(bt[[i]]) * zeta[i,k])+
          v * (zeta[i,k] - W[i,k] + c[i,k]/v)
      }
    }
    y.vir <- zeta - h * H.der
    Zeta <- ifelse(abs(y.vir) < scad.l * (1 + h), 1, 0) * ifelse(y.vir > 0, 1, -1) * 
      ifelse(abs(y.vir) - scad.l > 0, abs(y.vir) - scad.l, 0)+
      y.vir * ifelse(abs(y.vir) > scad.a * scad.l, 1, 0) + 
      ifelse(abs(y.vir) < scad.a * scad.l, 1, 0) * 
      ifelse(abs(y.vir) > (1 + h) * scad.l, 1, 0) * 
      ((scad.a - 1) * y.vir - ifelse(y.vir > 0, 1, -1) * scad.a * scad.l)/(scad.a - (1 + h))
    Zeta[which(Zeta[,1]<=0),1]<-0
    Zeta[which(Zeta[,2]<=0),2]<-0
    Zeta[which(Zeta[,3]<=0),3]<-0
    Zeta[which(Zeta[,4]<=0),4]<-0
    ## update sigma
    sigma.der <- NULL
    for (i in 1:sample.len) {
      sigma.der <- cbind(sigma.der, sum(diag(sigma.ni.1[[i]])))
    }
    Sigma <- sigma - mean(sigma.der) * step
    
    ## update eta
    eta.der<-array(NA,c(nbasis.u1, sample.len, K))
    for (k in 1:K) {
      for (i in 1:sample.len) {
        veceta1<-list()
        z = NULL
        for (j in 1:nbasis.u1) {
          z <- c(z, as.vector(2 * bt[[i]] %*% diag(nbasis.u1)[ ,j] %*% 
            (bu1.0[[k]][i, ] %*% lambda[,1])^2 %*% t(eta[ ,k]) %*% t(bt[[i]])))
        }
        eta.der[ ,i,k] <- as.vector(sigma.ni.1[[i]]) %*% matrix(z, time.len[i]^2,nbasis.u1)
    }
      Eta[ ,k] <- eta[ ,k] + (rowMeans(eta.der[,,k])) * step
    }
    Eta <- qr.Q(qr(Eta))
    for (k in 1:K) { if(Eta [1,k] < 0) { Eta[ ,k] = -Eta[ ,k] } }
   
    ## update c
    C <- c + v * (Zeta - W)
    err<-c(err, max(max(abs(beta-Beta)), max(abs(Zeta-zeta))))
    iter<-iter+1  
  }
}




