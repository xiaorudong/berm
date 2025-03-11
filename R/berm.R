
#' Bootstrap-Enhanced Regularization Method (BERM)
#' This function applies a bootstrap-enhanced regularization technique
#' integrating bootstrapped confidence intervals with elastic net models to
#' obtain robust variable selection and accurate coefficient estimation.
#'
#' @param x A matrix or data frame of predictor variables.
#' @param y A vector of the response variable.
#' @param lambda_grid A numeric vector of lambda values to be used in regularization.
#' @param alpha_grid A numeric vector specifying alpha values to be used in regularization.
#'                   This parameter only affects the results when 'unrestricted' is TRUE.
#' @param standardize Logical, indicates whether to standardize predictor variables before fitting the model.
#'                    Default is standardize=TRUE.
#' @param nfold The number of folds for cross-validation, used to optimize model parameters. Default is nfold=10.
#' @param K An integer specifying the number of bootstrap samples to draw. Default is K=100.
#' @param unrestricted Logical, if TRUE, allows 'alpha_grid' to vary; if FALSE, 'alpha_grid' is fixed at 0.5. Default is FALSE.
#' @return A list containing several components:
#'         - 'step1res': A data frame of coefficients across all bootstrap samples.
#'         - 'step1select': A data frame indicating which variables were consistently selected across samples.
#'         - 'step2pars': Model parameters from the final elastic net model.
#'         - 'coef': The estimated coefficients.
#'         - 'intercept': The estimated intercept.
#'         - 'finalmodel': The final model object.
#' @import glmnet
#' @import caret
#' @importFrom dplyr full_join
#' @importFrom purrr reduce
#' @importFrom tibble column_to_rownames
#' @importFrom stats quantile coef
#' @importFrom magrittr %>%

berm <- function(x, y,
                       lambda_grid = 10^seq(2, -3, by = -.1),
                       alpha_grid = seq(1, 0, length.out=20),
                       standardize = TRUE, nfold=10, K=100,
                       unrestricted=FALSE) {
  
  x <- data.matrix(x)
  y <- data.matrix(y)
  
  res_ls <- list()
  # pars_ls <- list()
  for (i in 1:K) {
    #Creating a resampled dataset from the sample data
    train_rows <- sample(1:nrow(x), nrow(x), replace = TRUE)
    x.train.sub <- x[train_rows, ]
    y.train.sub <- y[train_rows,]

    #Running the regression on these data
    if(unrestricted) test_alpha <- alpha_grid
    if(!unrestricted) test_alpha <- 0.5

    model_bootstrap <- ElasticNet.model(x=x.train.sub, y=y.train.sub,
                                        lambda_grid=lambda_grid, alpha_grid=test_alpha,
                                        standardize = standardize, nfold=nfold)
    res_ls[[i]] <- model_bootstrap$coef
    # pars_ls[[i]] <- model_bootstrap$pars
  }

  res_df <- res_ls %>% purrr::reduce(full_join, by = "Predictor") %>% column_to_rownames(var = "Predictor") %>% as.data.frame()
  colnames(res_df) <- 1:K
  ci_df <- data.frame(VAR=rownames(res_df),
                      lw=apply(res_df, 1, function(x) quantile(x, probs=0.025)),
                      up=apply(res_df, 1, function(x) quantile(x, probs=0.975)))

  ci_df$cover0 <- ifelse(ci_df$up>0 & ci_df$lw<0, T, F)
  selectedvars <- ci_df[which(!ci_df$cover0),]$VAR
  penalty.factor <- ifelse(colnames(x) %in% selectedvars, 1, Inf)

  allresLasso <- ElasticNet_weight.model(x, y,
                                         lambda_grid=lambda_grid, alpha_grid=alpha_grid,
                                         standardize = standardize, nfold=nfold,
                                         penalty.factor=penalty.factor)
  return(list(step1res=res_df, step1select=ci_df,
              step2pars = allresLasso$pars,
              coef=allresLasso$coef, intercept=allresLasso$intercept,
              finalmodel=allresLasso$finalmodel))
}


#' Elastic Net Model
#'
#' This function fits an elastic net model to the provided data using cross-validation
#' to select the best lambda and alpha values.
#'
#' @param x A matrix or data frame of predictor variables where rows are observations
#'          and columns are variables.
#' @param y A vector of the response variable.
#' @param lambda_grid A numeric vector of lambda values to be tested.
#' @param alpha_grid A numeric vector of alpha values to be tested.
#' @param standardize A logical value indicating whether the predictor variables should
#'                    be standardized before fitting the model.
#' @param nfold The number of folds to be used in cross-validation. This parameter
#'              determines how the data is split during the cross-validation process.
#' @return A list containing several components:
#'         - 'cvres': A data frame containing the cross-validation results.
#'         - 'pars': The best tuning parameters (lambda and alpha) selected based on
#'                   cross-validation.
#'         - 'coef': A data frame of coefficient estimates at the best lambda and alpha.
#'         - 'intercept': The intercept term from the model at the best lambda.
#'         - 'cvmodels': The full cross-validated `train` object including all models
#'                       fitted during the tuning process.
#'         - 'finalmodel': The final model object with parameters set at the best lambda
#'                         and alpha.
#' @import glmnet
#' @import caret
#' @importFrom stats coef

ElasticNet.model <- function(x, y,
                             lambda_grid = 10^seq(2, -3, by = -.1),
                             alpha_grid = seq(1, 0, length.out=20),
                             standardize = TRUE, nfold=10) {

  srchGrid <- expand.grid(alpha = alpha_grid, lambda = lambda_grid)

  x.train.sub <- data.matrix(x)
  y.train.sub <- data.matrix(y)

  model <- train(
    x=x.train.sub, y = y.train.sub[,1], method = "glmnet", standardize = standardize,
    trControl = trainControl("cv", number = nfold),
    tuneGrid = srchGrid
  )


  model.coef <- as.matrix(coef(model$finalModel, model$bestTune$lambda, exact = T,
                               x=x.train.sub, y = y.train.sub[,1]))
  coefs <- data.frame(Predictor = names(model.coef[,1]), Value= model.coef[,1])[-1,]

  return(list(cvres=model$results, pars = model$bestTune,
              coef=coefs, intercept=model.coef[1],
              cvmodels=model, finalmodel=model$finalModel))
}


#' Weighted Elastic Net Model
#'
#' This function extends the elastic net model by allowing for differential
#' penalization of variables through the `penalty.factor` parameter. It fits
#' the model using cross-validation to select the best lambda and alpha values,
#' incorporating these weights.
#'
#' @param x A matrix or data frame of predictor variables where rows are observations
#'          and columns are variables.
#' @param y A vector of the response variable.
#' @param lambda_grid A numeric vector of lambda values to be tested.
#' @param alpha_grid A numeric vector of alpha values to be tested.
#' @param standardize Logical, indicating whether the predictor variables should
#'                    be standardized before fitting the model.
#' @param nfold The number of folds to be used in cross-validation. This parameter
#'              determines how the data is split during the cross-validation process.
#' @param penalty.factor A numeric vector or a single number that specifies the penalty
#'                       multiplier for each predictor variable.
#' @return A list containing several components:
#'         - 'cvres': A data frame containing the cross-validation results.
#'         - 'pars': The best tuning parameters (lambda and alpha) selected based on
#'                   cross-validation.
#'         - 'coef': A data frame of coefficient estimates at the best lambda and alpha,
#'                   adjusted for the penalty factor.
#'         - 'intercept': The intercept term from the model at the best lambda.
#'         - 'cvmodels': The full cross-validated `train` object including all models
#'                       fitted during the tuning process.
#'         - 'finalmodel': The final model object with parameters set at the best lambda
#'                         and alpha.
#' @import glmnet
#' @import caret
#' @importFrom stats coef

ElasticNet_weight.model <- function(x, y,
                                    lambda_grid = 10^seq(2, -3, by = -.1),
                                    alpha_grid = seq(1, 0, length.out=20),
                                    standardize = TRUE, nfold=10,
                                    penalty.factor) {

  srchGrid <- expand.grid(alpha = alpha_grid, lambda = lambda_grid)

  x.train.sub <- data.matrix(x)
  y.train.sub <- data.matrix(y)

  model <- train(
    x=x.train.sub, y = y.train.sub[,1], method = "glmnet", standardize = standardize,
    trControl = trainControl("cv", number = nfold), penalty.factor=penalty.factor,
    tuneGrid = srchGrid
  )

  model.coef <- as.matrix(coef(model$finalModel, model$bestTune$lambda, exact = T,
                               x=x.train.sub, y = y.train.sub[,1], penalty.factor=penalty.factor))
  coefs <- data.frame(Predictor = names(model.coef[,1]), Value= model.coef[,1])[-1,]

  return(list(cvres=model$results, pars = model$bestTune,
              coef=coefs, intercept=model.coef[1],
              cvmodels=model, finalmodel=model$finalModel))
}


