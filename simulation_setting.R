########################
#### load libraries ####
########################


sim.example1 <- function(sample.len, time.len, p, sigma, beta, alpha){
  
  # sample.len: sample sizes 
  # time.len: number of time points for each individual
  # p: dimension of predictors
  # sigma: sd of error (response)
  # beta: index coefficients in mean function with norm 1
  # alpha: index coefficients in covariance function with norm 1
  
  # generate time points
  t <- runif(sum(time.len))
  time.index <- c(0,cumsum(time.len))
  time <- list()
  for (i in 1:sample.len) {time[[i]] <- t[(time.index[i]+1):time.index[i+1]]}

  # eigenfunctions
  eigenfun1 <- function(t) {cos(pi * t) * sqrt(2)}
  eigenfun2 <- function(t) {sin(pi * t) * sqrt(2)}
  eigenfun3 <- function(t) {cos(3 * pi * t) * sqrt(2)}

  
  # generate predictor data
  x <- matrix(runif(p * sample.len, -1, 1), sample.len, p)
  x <- apply(x, 2, scale)
  u <- x %*% beta
  uu <- x %*% alpha
  
  # true mean function 
  mu <- function(t, u) {10 * (u * cos(t) + (1 - u) * sin(t))}
  mean <- list()
  for (i in 1:sample.len) {mean[[i]] <- mu(time[[i]], u[[i]])}
    
  # true covariance function 
  rho1 <- function(x) {(1 - cos(x*pi))*10}
  rho2 <- function(x) {(1 - cos(x*pi))}
  rho3 <- function(x) {(1 - cos(x*pi))/10}
  
  # covariance function of each individual
  Sigma <- function(i) {
      eigenfun1(time[[i]]) %*% t(eigenfun1(time[[i]])) * rho1(uu[i,1])+
      eigenfun2(time[[i]]) %*% t(eigenfun2(time[[i]])) * rho2(uu[i,2])+ 
      eigenfun3(time[[i]]) %*% t(eigenfun3(time[[i]])) * rho3(uu[i,3])+
      diag(sigma, time.len[i], time.len[i])
  }

  # generate response function
  y <- list()
  for (i in 1:sample.len) {y[[i]] <- mvrnorm(1, mean[[i]], Sigma(i))}
  
  return(list(time=time, x=x, y=y, mean=mean, 
              rho1=rho1(uu[,1]), rho2=rho2(uu[,2]), rho3=rho3(uu[,3])))
}


sim.example2 <- function(sample.len, time.len, p, sigma, beta, alpha){
  
  # sample_len: sample sizes
  # time_len: number of time points
  # p: dimension of predictors
  # sigma: sd of error (response)
  # beta: index coefficients in mean function 
  # alpha: index coefficients in covariance function
  
  # generate time points
  t <- runif(sum(time.len))
  time.index <- c(0,cumsum(time.len))
  time <- list()
  for (i in 1:sample.len) {time[[i]] <- t[(time.index[i]+1):time.index[i+1]]}
  
  # eigenfunctions
  eigenfun1 <- function(t) {cos(pi * t) * sqrt(2)}
  eigenfun2 <- function(t) {sin(pi * t) * sqrt(2)}
  eigenfun3 <- function(t) {cos(3 * pi * t) * sqrt(2)}
  
  
  # generate predictor data
  x <- matrix(runif(p * sample.len, -1, 1), sample.len, p)
  x <- apply(x, 2, scale)
  u <- x %*% beta
  uu <- x %*% alpha
  
  # true mean function 
  mu <- function(t, u) {10 * (u * cos(t) + (1 - u) * sin(t))}
  mean <- list()
  for (i in 1:sample.len) {mean[[i]] <- mu(time[[i]], u[[i]])}
  
  # true eigenvalue
  rho<-t(matrix(c(5,1,0.5), p, sample.len))
  
  # covariance function of each individual
  Sigma <- function(i) {
    eigenfun1(time[[i]]) %*% t(eigenfun1(time[[i]])) * rho[i,1]+
      eigenfun2(time[[i]]) %*% t(eigenfun2(time[[i]])) * rho[i,2]+ 
      eigenfun3(time[[i]]) %*% t(eigenfun3(time[[i]])) * rho[i,3]+
      diag(sigma, time.len[[i]], time.len[[i]])
  }
  
  # generate response function
  y <- list()
  for (i in 1:sample_len) {y[[i]] <- mvrnorm(1, mean[[i]], Sigma(i))}
  
  return(list(time=time, x=x, y=y, mean=mean, 
              rho1=rho[,1], rho2=rho[,2], rho3=rho[,3]))
}


sim.example3 <- function(sample.len, time.len, p, beta, alpha){
  
  # sample_len: sample sizes
  # time_len: number of time points
  # p: dimension of predictors
  # sigma: sd of error (response)
  # beta: index coefficients of score function 
  # alpha: index coefficients in covariance function
  
  # generate time points
  t <- runif(sum(time.len))
  time.index <- c(0,cumsum(time.len))
  time <- list()
  for (i in 1:sample_len) {time[[i]] <- t[(time.index[i]+1):time.index[i+1]]}
  
  # eigenfunctions
  eigenfun1 <- function(t) {cos(pi * t) * sqrt(2)}
  eigenfun2 <- function(t) {sin(pi * t) * sqrt(2)}
  eigenfun3 <- function(t) {cos(3 * pi * t) * sqrt(2)}
  
  
  # generate predictor data
  x <- matrix(runif(p * sample.len, -1, 1), sample.len, p)
  x <- apply(x, 2, scale)
  u <- x %*% beta
  uu <- x %*% alpha
  
  # true score function 
  score1 <- function(x) {(1 - cos(x*pi))}
  score2 <- function(x) {(1 - cos(x*pi))/5}
  score3 <- function(x) {(1 - cos(x*pi))/10}
  
  # true mean function 
  mu <- function(t) {t^2 + 1} 
  mean <- list()
  for (i in 1:sample_len) {mean[[i]] <- mu(time[[i]])+
    eigenfun1(time[[i]])*score1(u[i,1]) + 
    eigenfun2(time[[i]])*score2(u[i,2]) + 
    eigenfun3(time[[i]])*score3(u[i,3])
  }
  

  # true eigenvalue function 
  eigenvalue1 <- function(x) {sqrt((1 - cos(x*pi)))}
  eigenvalue2 <- function(x) {sqrt((1 - cos(x*pi))/5)}
  eigenvalue3 <- function(x) {sqrt((1 - cos(x*pi))/10)}
  
  # covariance function of each individual
  Sigma <- function(i) {
    eigenfun1(time[[i]]) %*% t(eigenfun1(time[[i]])) * eigenvalue1(uu[i,1])+
      eigenfun2(time[[i]]) %*% t(eigenfun2(time[[i]])) * eigenvalue1(uu[i,2])+ 
      eigenfun3(time[[i]]) %*% t(eigenfun3(time[[i]])) * eigenvalue1(uu[i,3])
  }
  
  # generate response function
  y <- list()
  for (i in 1:sample.len) {y[[i]] <- mvrnorm(1, mean[[i]], Sigma(i))}
  
  return(list(time=time, x=x, mean=mean, y=y))
}

