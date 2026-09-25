rm(list = ls())
library(MASS)
library(truncnorm)
library(BayesLogit)
## Download and install the skinnybasad package from the supplementary materials of the paper
## "Skinny Gibbs: A Consistent and Scalable Gibbs Sampler for Model Selection"
## link: https://www.tandfonline.com/doi/suppl/10.1080/01621459.2018.1482754?scroll=top
library(skinnybasad)
library(glmnet)
library(ncvreg)

source("R/gibbs_skinny_pg.R")
source("R/gibbs_exact_pg.R")
source("R/evaluation_metrics.R")

################################################################################
### Driver script: loads simulated datasets (produced by generate_datasets.R),
### fits each competing method, evaluates variable-selection / prediction
### performance, and aggregates results across replicates.
################################################################################

nsim <- 20

HSGPG_results <- matrix(rep(0, nsim * 9), ncol = 9) ## Hierarchical Skinny Gibbs with Polya-Gamma
HEGPG_results <- matrix(rep(0, nsim * 9), ncol = 9) ## Hierarchical Exact Gibbs with Polya-Gamma
SGT_results <- matrix(rep(0, nsim * 9), ncol = 9)   ## Skinny Gibbs with t-approximation
EGT_results <- matrix(rep(0, nsim * 9), ncol = 9)   ## Exact Gibbs with t-approximation
lasso_results <- matrix(rep(0, nsim * 9), ncol = 9)
mcp_results <- matrix(rep(0, nsim * 9), ncol = 9)
scad_results <- matrix(rep(0, nsim * 9), ncol = 9)

for (i in 1:nsim) {

  # NOTE: this path corresponds to the p = 500, 8-active-covariate, set 2,
  # case 2, corr = 0.5 scenario. Update it to match whichever scenario you
  # generated with generate_datasets.R.
  load(file = paste("data/p500_8set2case2_0.5/dflist_", i, ".RData", sep = ""))

  p.all <- dflist$p     # total no. of covariates
  p.act <- dflist$p_act # no. of active covariates
  K <- dflist$K
  ind <- sample(1:p.all, K)
  gammaZ0 <- rep(0, p.all); gammaZ0[ind] <- 1 ### initial value for gammaZ
  n0 <- length(dflist$E)
  ################################ Models
  hexact_pg <- hegpg_func(X = as.matrix(dflist$X), E = unlist(dflist$E), gammaZ = gammaZ0,
                           q = dflist$qn, r = 2, s = 1, nburn = 2000, niter = 4000, cutoff = 0.5,
                           intercept = FALSE, standardize = FALSE)

  hskinny_pg <- hsgpg_func(X = as.matrix(dflist$X), E = unlist(dflist$E), gammaZ = gammaZ0,
                            q = dflist$qn, r = 2, s = 1, nburn = 2000, niter = 4000, cutoff = 0.5,
                            intercept = FALSE, standardize = FALSE)

  skinny_t <- skinnybasad(X = as.matrix(dflist$X), E = dflist$E, pr = dflist$qn, B0 = dflist$B0, Z0 = gammaZ0, nsplit = 10, a0 = 0.01, b0 = 1,
                           modif = 1, nburn = 2000, niter = 4000, printitrsep = 2000, maxsize = max(K, sqrt(n0)))

  exact_t <- skinnybasad(X = as.matrix(dflist$X), E = unlist(dflist$E), pr = dflist$qn, B0 = dflist$B0, Z0 = gammaZ0, nsplit = 10, a0 = 0.01, b0 = 1,
                          modif = 0, nburn = 2000, niter = 4000, printitrsep = 2000, maxsize = max(K, sqrt(n0)))

  lasso_f <- cv.glmnet(x = as.matrix(dflist$X), y = unlist(dflist$E), alpha = 1, intercept = TRUE,
                       family = "binomial", standardize = FALSE)

  mcp_f <- ncvreg(X = as.matrix(dflist$Xunscaled), y = unlist(dflist$E), family = "binomial", penalty = "MCP", gamma = 3,
                   nlambda = 50, warn = FALSE)

  scad_f <- ncvreg(X = as.matrix(dflist$Xunscaled), y = unlist(dflist$E), family = "binomial", penalty = "SCAD", gamma = 3.7,
                    nlambda = 50, warn = FALSE)

  z_act <- c(rep(1, p.act), rep(0, p.all - p.act))

  if (length(hskinny_pg$gammahat) > p.all) {
    hskinny_pg_pred <- hskinny_pg$gammahat[-1]
    hexact_pg_pred <- hexact_pg$gammahat[-1]
  } else {
    hskinny_pg_pred <- hskinny_pg$gammahat
    hexact_pg_pred <- hexact_pg$gammahat
  }

  skinny_t_pred <- ifelse(skinny_t$marZ >= 0.5, 1, 0)
  exact_t_pred <- ifelse(exact_t$marZ >= 0.5, 1, 0)
  if (length(coef(lasso_f)) > p.all) {
    lasso_pred <- ifelse(coef(lasso_f, s = lasso_f$lambda.min) != 0, 1, 0)[-1]
  } else {
    lasso_pred <- ifelse(coef(lasso_f, s = lasso_f$lambda.min) != 0, 1, 0)
  }
  mcp_pred <- ifelse(mcp_f$beta[, mcp_f$convex.min][-1] != 0, 1, 0)
  scad_pred <- ifelse(scad_f$beta[, scad_f$convex.min][-1] != 0, 1, 0)

  if (length(hskinny_pg$betahat) > p.all) {
    eval.hskinnypg <- evaluation(actual = z_act, predicted = hskinny_pg_pred, beta = hskinny_pg$betahat[-1],
                                  X = as.matrix(dflist$X), Y = dflist$E, method = "logit")
    eval.hexactpg <- evaluation(actual = z_act, predicted = hexact_pg_pred, beta = hexact_pg$betahat[-1],
                                 X = as.matrix(dflist$X), Y = dflist$E, method = "logit")
  } else {
    eval.hskinnypg <- evaluation(actual = z_act, predicted = hskinny_pg_pred, beta = hskinny_pg$betahat,
                                  X = as.matrix(dflist$X), Y = dflist$E, method = "logit")
    eval.hexactpg <- evaluation(actual = z_act, predicted = hexact_pg_pred, beta = hexact_pg$betahat,
                                 X = as.matrix(dflist$X), Y = dflist$E, method = "logit")
  }

  eval.skinnyt <- evaluation(actual = z_act, predicted = skinny_t_pred, beta = NULL, X = as.matrix(dflist$X),
                              Y = dflist$E, method = "logit")
  eval.exactt <- evaluation(actual = z_act, predicted = exact_t_pred, beta = NULL, X = as.matrix(dflist$X),
                             Y = dflist$E, method = "logit")
  eval.lasso <- evaluation(actual = z_act, predicted = lasso_pred, beta = as.numeric(coef(lasso_f, s = lasso_f$lambda.min))[-1],
                            X = as.matrix(dflist$X), Y = dflist$E, method = "logit")
  eval.mcp <- evaluation(actual = z_act, predicted = mcp_pred, beta = mcp_f$beta[, mcp_f$convex.min][-1],
                          X = as.matrix(dflist$X), Y = dflist$E, method = "logit")
  eval.scad <- evaluation(actual = z_act, predicted = scad_pred, beta = scad_f$beta[, scad_f$convex.min][-1],
                           X = as.matrix(dflist$X), Y = dflist$E, method = "logit")

  HSGPG_results[i, ] <- c(unlist(eval.hskinnypg), length(dflist$E[dflist$E == 1]) / length(dflist$E))
  HEGPG_results[i, ] <- c(unlist(eval.hexactpg), length(dflist$E[dflist$E == 1]) / length(dflist$E))
  SGT_results[i, ] <- c(unlist(eval.skinnyt), length(dflist$E[dflist$E == 1]) / length(dflist$E))
  EGT_results[i, ] <- c(unlist(eval.exactt), length(dflist$E[dflist$E == 1]) / length(dflist$E))
  lasso_results[i, ] <- c(unlist(eval.lasso), length(dflist$E[dflist$E == 1]) / length(dflist$E))
  mcp_results[i, ] <- c(unlist(eval.mcp), length(dflist$E[dflist$E == 1]) / length(dflist$E))
  scad_results[i, ] <- c(unlist(eval.scad), length(dflist$E[dflist$E == 1]) / length(dflist$E))

  if (i %% 10 == 0) print(paste("finish simulation replicate", i))

}

f_results <- data.frame(HSGPG = apply(HSGPG_results, 2, mean),
                         HEGPG = apply(HEGPG_results, 2, mean),
                         SGT = apply(SGT_results, 2, mean),
                         EGT = apply(EGT_results, 2, mean),
                         AL = apply(lasso_results, 2, mean),
                         MCP = apply(mcp_results, 2, mean),
                         SCAD = apply(scad_results, 2, mean),
                         row.names = c("SEN", "SPE", "MCC", "TP", "FP", "TN", "FN", "MSPE", "Ratio"))

# NOTE: previously this saved the last loop iteration's raw `dflist` instead of
# the aggregated results table -- fixed to save `f_results`.
save(f_results, file = paste("results/result_p_", p.all, "_", p.act, ".RData", sep = ""))
