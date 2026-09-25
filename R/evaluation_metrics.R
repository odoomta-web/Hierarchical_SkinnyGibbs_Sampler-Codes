################################################################################
### Evaluation metrics used to compare variable-selection / prediction
### performance across methods (Sensitivity, Specificity, MCC, MSPE).
################################################################################

evaluation <- function(actual, predicted, beta, X, Y, method = c("logit", "probit")) {

  true.idx <- which(actual == 1)
  false.idx <- which(actual == 0)
  positive.idx <- which(predicted == 1)
  negative.idx <- which(predicted == 0)

  TP <- length(intersect(true.idx, positive.idx))
  FP <- length(intersect(false.idx, positive.idx))
  FN <- length(intersect(true.idx, negative.idx))
  TN <- length(intersect(false.idx, negative.idx))

  Sensitivity <- TP / (TP + FN)
  if ((TP + FN) == 0) Sensitivity <- 1
  Specific <- TN / (TN + FP)
  if ((TN + FP) == 0) Specific <- 1

  MCC.denom <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  if (MCC.denom == 0) MCC.denom <- 1
  MCC <- (TP * TN - FP * FN) / MCC.denom
  if ((TN + FP) == 0) MCC <- 1

  n <- nrow(X)
  if (is.null(beta) == TRUE) {
    MSPE <- NA
  } else {
    if (method == "logit") MSPE <- mean((plogis(X %*% beta) - Y)^2)
    if (method == "probit") MSPE <- mean((pnorm(X %*% beta) - Y)^2)
  }

  return(list(Sensitivity = Sensitivity, Specific = Specific, MCC = MCC,
              TP = TP, FP = FP, TN = TN, FN = FN, MSPE = MSPE))
}
