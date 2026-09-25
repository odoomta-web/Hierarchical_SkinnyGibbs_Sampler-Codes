################################################################################
### Hierarchical EXACT Gibbs sampler for logistic regression using
### Polya-Gamma latent variables (H-EGPG).
################################################################################

library(Matrix)
library(MASS)
library(invgamma)
library(BayesLogit)
library(mvnfast)

hegpg_func <- function(X, E, gammaZ, q, r, s, niter = 4000, nburn = 2000, cutoff = 0.5,
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
  eps_diag <- .Machine$double.eps * diag(p)

  if (intercept == TRUE) {
    gammaZ[1] <- 1
    startj <- 2
  } else {
    startj <- 1
  }

  for (itr in 1:(nburn + niter)) {
    if (itr %% 1000 == 0) cat("Completed iteration", itr, "\n")

    Y <- E_minus_0_5 / w
    inv_tau1_2 <- 1 / tau1_2
    const <- (q * sqrt(tau0_2)) / ((1 - q) * sqrt(tau1_2))

    # Update beta
    D_z <- diag((gammaZ * inv_tau1_2) + (1 - gammaZ) * inv_tau0_2)
    V <- sweep(t(X), 2, w, FUN = "*") %*% X + D_z + eps_diag  # t(X) %*% W %*% X + D_z + eps_diag
    sigma1 <- solve(V)
    mean1 <- sweep(sigma1 %*% t(X), 2, w, FUN = "*") %*% Y     # sigma1 %*% t(X) %*% W %*% Y
    beta <- t(unlist(rmvn(1, mu = mean1, sigma = sigma1)))

    # Update gammaZ
    for (j in startj:p) {
      temp1 <- beta[j]^2 * 0.5 * (inv_tau0_2 - inv_tau1_2)
      log_d_j <- log(const) + temp1
      if (is.infinite(exp(log_d_j)) == TRUE) {
        gammaZ[j] <- 1
      } else {
        gammaZ[j] <- rbinom(1, 1, exp(log_d_j) / (1 + exp(log_d_j)))
      }
    }

    # Update tau1^2
    tau1_2 <- rinvgamma(1, shape = (sum(gammaZ) / 2 + r), rate = (s + sum(gammaZ * beta * beta) / 2))

    # Update Polya-Gamma latent variables (omega/w)
    w <- rpg.sp(n, 1, X %*% beta)

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
