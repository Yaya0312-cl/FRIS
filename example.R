###################
#### load codes ###
###################
source("FRIS.R")
source("simulation_setting.R")

###################################
#### generation of the example ####
###################################
sample.len <- 100  # training sample size
time.len <- rep(10,100)  # number of time points
p <- 3  # dimension of predictors
sigma <- 1  # sd of error (response)
K <- 3  # number of eigenfunction
beta <- c(0.2, 0.8, 0.6)  # index coefficients in mean function with norm 1
alpha <- matrix(c(0.9, 0.1, 0.4, 0.2, 0.6, 0.8, 0.5, 0.8, 0.3), p, K) 
  # index coefficients in covariance function with norm 1
sim.example <-sim.example1  # choice of example

dat <- sim.example(sample.len=sample.len, time.len=time.len, 
                    p=p, sigma=sigma, beta=beta, alpha=alpha)
time <- dat$time; x <- dat$x; mean <- dat$mean; y <- dat$y; 
rho1 <- dat$rho1; rho2 <- dat$rho2; rho3 <- dat$rho3; 

####get initial values
## beta, gamma
threshold <- 0.001
max.iter <- 100
spline.par <- list(rangeval.t = c(0, 1), norder.t = 3, nk.t = 6,
                   rangeval.u = c(-5, 5), norder.u = 3, nk.u = 6,
                   rangeval.t1 = c(0, 1), norder.t1 = 3, nk.t1 = 6,
                   rangeval.u1 = c(-5, 5), norder.u1 = 3, nk.u1 = 6)
nbisis.t <- spline.par$norder.t + spline.par$nk.t -2
nbisis.u <- spline.par$norder.u + spline.par$nk.u -2
beta.ini<-matrix(rep(0.6, p),3,1)
gamma.ini<-rnorm(nbisis.u * nbisis.t, 0, 1)
mu.ini <- get.mu.ini(time, x, y, beta.ini, gamma.ini, max.iter, threshold, spline.par)
Beta.ini <- mu.ini$beta.ini
Gamma.ini <- mu.ini$gamma.ini

## eta
K <- 4
eigenfun.ini <- get.eigenfun.ini(K, spline.par)
Eta.ini <- eigenfun.ini$eta.ini
Eigenfun.ini <- eigenfun.ini$eigenfun
grids <- eigenfun.ini$grids
residual <- eigenfun.ini$residual


## lambda, alpha
n.cl <- 20
scorevar.ini <- get.scorevar.ini(time, x, n.cl, spline.par)
Alpha.ini <- scorevar.ini$alpha.ini
Lambda.ini <- scorevar.ini$lambda.ini
Zeta.ini <- scorevar.ini$zeta.ini

## sigma
sigma.ini <- get.sigma.ini(time, grids, Eigenfun.ini, Zeta.ini, residual)
Sigma.ini <- sigma.ini$sigma.ini

#### run code
tuning.par <- list(v = 0.001, h = 0.0001, scad.l = 0.1, scad.a = 3.7, step = 0.001)
max.iter <- 50
threshold <- 1e-3
K <- 4
runcode <- FRIS(Beta.ini, Gamma.ini, Eta.ini, Alpha.ini, Lambda.ini, 
                 Sigma.ini, Zeta.ini, max.iter, threshold, spline.par, tuning.par, K)

