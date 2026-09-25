library(mvtnorm)

################################################################################
### Functions for simulating high-dimensional logistic regression data,
### used by generate_datasets.R.
################################################################################

### Maps a (design, p.act) combination to (p, p_active).
### design: 1 -> p = 250, 2 -> p = 500, 3 -> p = 1000
### p.act:  0 -> 0 active covariates, 1 -> 4 active, 2 -> 8 active
gen_p <- function(design, p.act) {
  if (design == 1) p <- 250
  if (design == 2) p <- 500
  if (design == 3) p <- 1000
  if (p.act == 0) {
    p_active <- 0
  } else if (p.act == 1) {
    p_active <- 4
  } else if (p.act == 2) {
    p_active <- 8
  } else {
    stop("gen_p(): unsupported p.act value '", p.act,
         "'. Only 0 (0 active), 1 (4 active), or 2 (8 active) are defined.")
  }
  return(c(p, p_active))
}

### Generating X
gen_X <- function(n, design, p.act, case, corr) {
  p <- gen_p(design, p.act)[1]
  ar1 <- function(n, rho) rho^(abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - (1:n - 1)))
  # Case 1: Isotropic design, where Sigma = Ip.
  if (case == 1) sigm <- diag(p)
  # Case 2: AR(1)/CS-type design with correlation `corr`.
  if (case == 2) sigm <- ar1(p, corr)
  X <- rmvnorm(n = n, mean = rep(0, p), sigma = sigm)
  return(X)
}

### Generating Beta
gen_Beta <- function(design, p.act, set) {
  set.seed(2024)
  p <- gen_p(design, p.act)[1]
  p_active <- gen_p(design, p.act)[2]
  beta <- rep(0, p)
  if (p_active > 0) {
    # Set 1: active betas ~ Unif(0.5, 1.5)
    if (set == 1) beta[1:p_active] <- runif(n = p_active, min = 0.5, max = 1.5)
    # Set 2: active betas ~ Unif(1.5, 3)
    if (set == 2) beta[1:p_active] <- runif(n = p_active, min = 1.5, max = 3)
  }

  return(beta)
}

################################################################################
### Main simulation data generator.
### design: 1, 2, 3   (p = 250, 500, 1000)
### p.act:  0, 1, 2    (0, 4, 8 active covariates)
### set:    1 (weak signal), 2 (strong signal)
### case:   1 (independent covariates), 2 (AR(1)-correlated, needs `corr`)
################################################################################
simdata <- function(n, design, p.act, set, case, corr = NULL) {
  p <- gen_p(design, p.act)[1]
  pact <- gen_p(design, p.act)[2]
  if (case == 2 && is.null(corr) == TRUE) message("corr cannot be null under case 2")
  X0 <- gen_X(n, design, p.act, case, corr)

  Bc <- gen_Beta(design, p.act, set)
  ## Logit probabilities
  Y <- plogis(X0 %*% Bc)
  ## Binary response variable
  E <- rbinom(n, 1, Y)
  Xunscaled <- as.matrix(X0)
  X0 <- as.matrix(scale(X0)) ### scaling
  X <- as.matrix(X0)
  return(list(X = X0, Y = Y, E = E, B = Bc, p_act = pact, p = p, Xunscaled = Xunscaled,
              case = case, set = set, simtype = "main"))
}


##################### Supplementary data simulation ######################

gen_Beta_sup <- function(design, p.act) {
  set.seed(2024)
  p <- gen_p(design, p.act)[1]
  pact <- gen_p(design, p.act)[2]
  b0 <- rep(0, p)
  if (pact == 0) { b0 <- b0 }
  if (pact == 4) { b0[1:pact] <- c(-1.5, 2.0, -2.5, 3.0) }
  if (pact == 8) { b0[1:pact] <- c(-2, 2.143, -2.286, 2.429, -2.571, 2.714, -2.857, 3) }
  return(b0)
}

logit_data <- function(n, design, p.act, rho_11, rho_12, rho_33) {
  p <- gen_p(design, p.act)[1]
  pact <- gen_p(design, p.act)[2]
  p_active <- pact
  Bc <- gen_Beta_sup(design, p.act)
  cov.1 <- (1 - rho_11) * diag(pact) + array(rho_11, c(pact, pact))       ## cov matrix of active covariates
  cov.3 <- (1 - rho_33) * diag(p - pact) + array(rho_33, c(p - pact, p - pact)) ## cov matrix of inactive covariates
  cov.2 <- array(rho_12, c(pact, p - pact))                               ## cov between active and inactive
  cov.f <- rbind(cbind(cov.1, cov.2), cbind(t(cov.2), cov.3))             ### full covariance matrix
  ## Eigen decomposition of the covariance matrix to get the square-root matrix
  cov.E <- eigen(cov.f)
  cov.sq <- cov.E$vectors %*% diag(sqrt(cov.E$values)) %*% t(cov.E$vectors)
  X00 <- matrix(rnorm(n * p), nrow = n)
  X0 <- cov.sq %*% t(X00) # covariates with the induced correlation structure
  X0 <- t(X0)
  Y <- plogis(X0 %*% Bc)
  E <- rbinom(n, 1, Y)
  Xunscaled <- as.matrix(X0)
  X0 <- as.matrix(scale(X0)) ### scaling
  X <- as.matrix(X0)
  return(list(X = X, Y = Y, E = E, B = Bc, p_act = p_active, p = p, Xunscaled = Xunscaled,
              simtype = "supp", rho_act = rho_11, rho_inact = rho_33, rho_bn = rho_12))
}
