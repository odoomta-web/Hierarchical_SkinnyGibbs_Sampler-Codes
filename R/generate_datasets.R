rm(list = ls())

source("R/simulate_data.R")

################################################################################
### Driver script: generates and saves the simulated datasets used by
### run_model_comparison.R. Each replicate is written to its own .RData file
### under data/<scenario>/dflist_<i>.RData -- create that directory before
### running (see README.md).
################################################################################

sim <- c("main", "supp") ##
nsim <- 20  # number of simulation replicates
n0 <- 200   # number of observations
d0 <- c(1, 2, 3)       ### design index: maps to total no. of covariates (250, 500, 1000) via gen_p()
pact0 <- c(0, 1, 2)    ### no. of active-covariate index: 0 -> 0 active, 1 -> 4 active, 2 -> 8 active (see gen_p())
set0 <- c(1, 2)        # 1 = weak signal; 2 = strong signal
case0 <- c(1, 2)       # 1 = independent; 2 = AR(1)-correlated
corr0 <- c(0.5, 0.75, 0.90)

rho_110 <- c(0, 0.25, 0.5, 0.75, 0.9) ## Correlation between active covariates
rho_120 <- c(0, 0.25, 0.5, 0.75, 0.9) # Correlation between active and inactive covariates
rho_330 <- c(0, 0.25, 0.5, 0.75, 0.9) # Correlation between inactive covariates

K <- max(30, log(n0))

## c0: index 1 = main simulation, 2 = supplementary simulation
## c1: design index into d0 (total number of covariates)
## c2: index into pact0 (number of active covariates)
##
## NOTE: the indices below (c1 = 1, c2 = 3) reproduce a p = 250, 8-active-covariate
## scenario. run_model_comparison.R's example load path is
## "data/p500_8set2case2_0.5/dflist_<i>.RData" (p = 500, 8 active, set 2, case 2,
## corr 0.5) -- set c1 = 2 below (design 2 -> p = 500) if you want the datasets
## to match that example exactly, and update the load path in
## run_model_comparison.R if you choose a different scenario.
for (c0 in c(2)) {          ### main or supplementary simulation
  for (c1 in c(1)) {        # design index: number of total covariates
    for (c2 in c(3)) {      # index of number of active covariates
      d <- d0[c1]
      pact <- pact0[c2]
      if (c0 == 1) { # Main simulation
        for (c3 in c(2)) {   # setting index: weak vs strong signal
          for (c4 in c(1:2)) { # case index: independent or correlated
            for (c5 in c(1)) { # correlation value for AR(1) design
              s <- set0[c3]
              cc <- case0[c4]
              corval <- corr0[c5]
              for (i in 1:(nsim)) {
                dflist <- simdata(n = n0, design = d, p.act = pact, set = s, case = cc, corr = corval)
                choicep <- function(x) { return(x - K + qnorm(0.90) * sqrt((x / dflist$p) * (1 - x / dflist$p))) }
                cp <- uniroot(choicep, c(1, K))$root
                qn <- cp / dflist$p
                dflist$index <- i
                dflist$qn <- qn
                dflist$K <- K
                dflist$B0 <- rep(0, dflist$p)
                print(i)
                save(dflist, file = paste("data/p", dflist$p, "_", dflist$p_act, "set", s, "case", cc, "_", corval, "/dflist_", i, ".RData", sep = ""))
              }
            } # End c5
          } # End c4
        } # End c3

      } else if (c0 == 2) { # Supplementary simulation

        for (c6 in c(2)) {   # correlation index for active covariates (and inactive, reused below)
          for (c8 in c(2)) { # correlation index between active and inactive covariates
            rho110 <- rho_110[c6]
            rho330 <- rho_330[c6]
            rho120 <- rho_110[c8]

            for (i in 1:(nsim)) {
              dflist <- logit_data(n = n0, design = d, p.act = pact, rho_11 = rho110, rho_12 = rho120, rho_33 = rho330)
              choicep <- function(x) { return(x - K + qnorm(0.90) * sqrt((x / dflist$p) * (1 - x / dflist$p))) }
              cp <- uniroot(choicep, c(1, K))$root
              qn <- cp / dflist$p
              dflist$index <- i
              dflist$qn <- qn
              dflist$K <- K
              dflist$B0 <- rep(0, dflist$p)
              print(paste("supp", i))
              save(dflist, file = paste("data/p", dflist$p, "_", dflist$p_act, "_", rho110, "_sup/dflist_", i, ".RData", sep = ""))
            }
          } # End c8
        } # End c6

      } # End if c0 statement
    } # End c2
  } # End c1

} ## end c0
