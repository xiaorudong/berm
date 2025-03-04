#' Evaluate Model Results
#'
#' This function calculates the R-squared value. The R^2 value provides an indication
#' of goodness of fit and therefore a measure of how well unseen samples are likely
#' to be predicted by the model, with values ranging from 0 to 1.
#'
#' @param true A numeric vector of observed true values of the dependent variable.
#' @param predicted A numeric vector of predicted values from the model, corresponding
#'                  to the true values.
#'
#' @return A single numeric value representing the R-squared of the model.


eval_results <- function(true, predicted) {
  SSE <- sum((predicted - true)^2)
  SST <- sum((true - mean(true))^2)
  R_square <- 1 - SSE / SST

  return(R_square)
}
