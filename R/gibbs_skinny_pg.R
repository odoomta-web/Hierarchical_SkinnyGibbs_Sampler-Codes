################################################################################
### Hierarchical SKINNY Gibbs sampler for logistic regression using
### Polya-Gamma latent variables (H-SGPG).
################################################################################

library(MASS)
library(invgamma)
library(BayesLogit)

hsgpg_func <- function(X, E, Z, gammaZ, q, r, s, niter = 4000, nburn = 2000, cutoff = 0.5,
                        intercept = c(TRUE, FALSE),
                        standardize = c(TRUE, FALSE)) {
  Xnames <- colnames(X)
  n <- nrow(X)
  if (standardize == TRUE) { X <- as.matrix(scale(X)) }

  if (intercept == TRUE) {
    X <- as.matrix(cbind(rep(1, n), X))
    gammaZ <- c(1, gammaZ)
    Xnames <- c(intercept, Xnames)
  }
  p <- ncol(X)
  outbeta <- numeric()
  outgammaZ <- numeric()
  outtau1_2 <- numeric()
  gammaSize <- numeric()

  tau0_2 <- 1 / n
  inv_tau0_2 <- 1 / tau0_2
  tau1_2 <- rinvgamma(1, shape = r, rate = s)
  E_minus_0_5 <- E - 0.5
  w <- rpg.sp(n, 1, 0)

  if (intercept == TRUE) {
    gammaZ[1] <- 1
    startj <- 2
  } else {
    startj <- 1
  }

  for (itr in seq_len(nburn + niter)) {
    if (itr %% 1000 == 0) cat("finish iteration", itr, "\n")

    Y <- E_minus_0_5 / w
    idx_active <- which(gammaZ == 1)
    idx_inactive <- which(gammaZ == 0)
    n_active <- length(idx_active)
    inv_tau1_2 <- 1 / tau1_2
    const <- (q * sqrt(tau0_2)) / ((1 - q) * sqrt(tau1_2))

    ##### Update beta
    beta <- numeric(p)
    if (n_active == 0) {
      beta <- rnorm(p, mean = 0, sd = 1 / sqrt(n - 1 + inv_tau0_2))
    } else {
      X_active <- X[, idx_active, drop = FALSE]

      XtWX <- sweep(t(X_active), 2, w, FUN = "*") %*% X_active + inv_tau1_2 * diag(n_active)
      sigma1 <- solve(XtWX)

      mean1 <- sweep(sigma1 %*% t(X_active), 2, w, FUN = "*") %*% Y
      beta[idx_active] <- t(unlist(rmvn(1, mu = mean1, sigma = sigma1)))
      beta[idx_inactive] <- rnorm(p - n_active, mean = 0, sd = 1 / sqrt(n - 1 + inv_tau0_2))
    }

    ##### Update gammaZ
    if (n_active > 0) {
      temp_2_12 <- Y - X_active %*% beta[idx_active]
    } else {
      temp_2_12 <- Y - numeric(n)
    }

    beta_active <- beta[idx_active]
    for (j in startj:p) {
      temp_0 <- 0.5 * beta[j]^2 * (inv_tau0_2 - inv_tau1_2)
      temp_2_1 <- beta[j] * X[, j] * w
      if (j %in% idx_active) {
        ### residual term excluding covariate j from the active set
        temp_2_11 <- (Y - X_active[, -which(idx_active == j), drop = FALSE] %*% beta_active[-which(idx_active == j)])
        temp_2 <- sum(temp_2_1 * temp_2_11)
      } else {
        temp_2 <- sum(temp_2_1 * temp_2_12)
      }

      # Equivalent to, but far cheaper than, sum(X[,j] * diag(1 - w) * X[,j])
      temp_3 <- 0.5 * (beta[j]^2) * sum(X[, j]^2 * (1 - w))
      log_d_j <- log(const) + temp_0 + temp_2 + temp_3

      if (is.infinite(exp(log_d_j))) {
        gammaZ[j] <- 1
      } else {
        gammaZ[j] <- rbinom(1, 1, exp(log_d_j) / (1 + exp(log_d_j)))
      }
    }

    # Update tau1^2
    tau1_2 <- rinvgamma(1, shape = (sum(gammaZ) / 2 + r), rate = (s + sum(gammaZ * beta * beta) / 2))

    # Update Polya-Gamma latent variables (omega/w)
    if (n_active > 0) {
      w <- rpg.sp(n, 1, X_active %*% beta_active)
    }

    if (itr > nburn) {
      outbeta <- cbind(outbeta, unlist(beta))
      outgammaZ <- rbind(outgammaZ, gammaZ)
      gammaSize[(itr - nburn)] <- sum(gammaZ)
      outtau1_2[(itr - nburn)] <- tau1_2
    }
  }

  gammaEst <- ifelse(apply(outgammaZ, 2, mean) >= cutoff, 1, 0)
  betaEst <- apply(outbeta, 1, mean)
  betaEst <- betaEst * gammaEst
  tau1_sq <- mean(outtau1_2)
  return(list("betahat" = betaEst, "gammahat" = gammaEst, "tau1_sqhat" = tau1_sq,
              "outbeta" = outbeta, "outgamma" = outgammaZ, "outtau1_sq" = outtau1_2,
              "size" = gammaSize))
}
