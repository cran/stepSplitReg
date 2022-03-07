#' 
#' @title Coefficients for stepSplitReg Object
#' 
#' @description \code{coef.stepSplitReg} returns the coefficients for a stepSplitReg object.
#' 
#' @method coef stepSplitReg
#' 
#' @param object An object of class stepSplitReg
#' @param group_index Groups included in the ensemble. Default setting includes all the groups.
#' @param ... Additional arguments for compatibility.
#' 
#' @return The coefficients for the stepSplitReg object.
#' 
#' @export
#' 
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#' 
#' @seealso \code{\link{stepSplitReg}}
#' 
#' @examples 
#' # Required Libraries
#' library(mvnfast)
#' 
#' # Setting the parameters
#' p <- 100
#' n <- 30
#' n.test <- 1000
#' sparsity <- 0.2
#' rho <- 0.5
#' SNR <- 3
#' 
#' # Generating the coefficient
#' p.active <- floor(p*sparsity)
#' a <- 4*log(n)/sqrt(n)
#' neg.prob <- 0.2
#' nonzero.betas <- (-1)^(rbinom(p.active, 1, neg.prob))*(a + abs(rnorm(p.active)))
#' 
#' # Correlation structure
#' Sigma <- matrix(0, p, p)
#' Sigma[1:p.active, 1:p.active] <- rho
#' diag(Sigma) <- 1
#' true.beta <- c(nonzero.betas, rep(0 , p - p.active))
#' 
#' # Computing the noise parameter for target SNR
#' sigma.epsilon <- as.numeric(sqrt((t(true.beta) %*% Sigma %*% true.beta)/SNR))
#' 
#' # Simulate some data
#' set.seed(1)
#' x.train <- mvnfast::rmvn(n, mu=rep(0,p), sigma=Sigma)
#' y.train <- 1 + x.train %*% true.beta + rnorm(n=n, mean=0, sd=sigma.epsilon)
#' x.test <- mvnfast::rmvn(n.test, mu=rep(0,p), sigma=Sigma)
#' y.test <- 1 + x.test %*% true.beta + rnorm(n.test, sd=sigma.epsilon)
#' 
#' # Stepwise Split Regularized Regression
#' step.out <- stepSplitReg(x.train, y.train, n_models = 3, max_variables = NULL, keep = 4/4,
#'                          model_criterion = c("F-test", "RSS")[1],
#'                          stop_criterion = c("F-test", "pR2", "aR2", "R2", "Fixed")[1], 
#'                          stop_parameter = 0.05, 
#'                          shrinkage = TRUE, alpha = 4/4, include_intercept = TRUE, 
#'                          n_lambda = 50, tolerance = 1e-2, max_iter = 1e5, n_folds = 5, 
#'                          model_weights = c("Equal", "Proportional", "Stacking")[1])
#' step.coefficients <- coef(step.out, group_index = 1:step.out$n_models)
#' step.predictions <- predict(step.out, x.test, group_index = 1:step.out$n_models)
#' mspe.step <- mean((step.predictions-y.test)^2)/sigma.epsilon^2
#' 
coef.stepSplitReg <- function(object, group_index = NULL, ...){
   
   if(is.null(group_index)){

      final_coef <- numeric(length(object$betas[[1]]) + 1)
      for(model.ind in 1:object$n_models)
         final_coef <- final_coef + c(object$intercepts[[model.ind]], object$betas[[model.ind]]) * object$models_weights[[model.ind]]
      return(final_coef)
   } else{
      
      if(any(!(group_index %in% 1:object$n_models)))
         stop("The group index is invalid.")
      
      final_coef <- numeric(length(object$betas[[1]]) + 1)
      for(model.ind in group_index)
         final_coef <- final_coef + c(object$intercepts[[model.ind]], object$betas[[model.ind]]) * object$models_weights[[model.ind]]
      return(final_coef)
   }
}
#' 
#' @title Coefficients for cv.stepSplitReg Object
#' 
#' @description \code{coef.cv.stepSplitReg} returns the coefficients for a cv.stepSplitReg object.
#' 
#' @method coef cv.stepSplitReg
#' 
#' @param object An object of class cv.stepSplitReg
#' @param group_index Groups included in the ensemble. Default setting includes all the groups.
#' @param ... Additional arguments for compatibility.
#' 
#' @return The coefficients for the cv.stepSplitReg object.
#' 
#' @export
#' 
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#' 
#' @seealso \code{\link{cv.stepSplitReg}}
#' 
#' @examples 
#' # Required Libraries
#' library(mvnfast)
#' 
#' # Setting the parameters
#' p <- 100
#' n <- 30
#' n.test <- 500
#' sparsity <- 0.2
#' rho <- 0.5
#' SNR <- 3
#' 
#' # Generating the coefficient
#' p.active <- floor(p*sparsity)
#' a <- 4*log(n)/sqrt(n)
#' neg.prob <- 0.2
#' nonzero.betas <- (-1)^(rbinom(p.active, 1, neg.prob))*(a + abs(rnorm(p.active)))
#' 
#' # Correlation structure
#' Sigma <- matrix(0, p, p)
#' Sigma[1:p.active, 1:p.active] <- rho
#' diag(Sigma) <- 1
#' true.beta <- c(nonzero.betas, rep(0 , p - p.active))
#' 
#' # Computing the noise parameter for target SNR
#' sigma.epsilon <- as.numeric(sqrt((t(true.beta) %*% Sigma %*% true.beta)/SNR))
#' 
#' # Simulate some data
#' set.seed(1)
#' x.train <- mvnfast::rmvn(n, mu=rep(0,p), sigma=Sigma)
#' y.train <- 1 + x.train %*% true.beta + rnorm(n=n, mean=0, sd=sigma.epsilon)
#' x.test <- mvnfast::rmvn(n.test, mu=rep(0,p), sigma=Sigma)
#' y.test <- 1 + x.test %*% true.beta + rnorm(n.test, sd=sigma.epsilon)
#' 
#' # Stepwise Split Regularized Regression
#' step.out <- cv.stepSplitReg(x.train, y.train, n_models = c(2, 3), max_variables = NULL, keep = 4/4,
#'                             model_criterion = c("F-test", "RSS")[1],
#'                             stop_criterion = c("F-test", "pR2", "aR2", "R2", "Fixed")[1], 
#'                             stop_parameter = 0.05, 
#'                             shrinkage = TRUE, alpha = 4/4, include_intercept = TRUE, 
#'                             n_lambda = 50, tolerance = 1e-2, max_iter = 1e5, n_folds = 5, 
#'                             model_weights = c("Equal", "Proportional", "Stacking")[1], 
#'                             n_threads = 1)
#' step.coefficients <- coef(step.out, group_index = 1:step.out$n_models_optimal)
#' step.predictions <- predict(step.out, x.test, group_index = 1:step.out$n_models_optimal)
#' mspe.step <- mean((step.predictions-y.test)^2)/sigma.epsilon^2
#' 
coef.cv.stepSplitReg <- function(object, group_index = NULL, ...){
   
   if(is.null(group_index)){
      
      final_coef <- numeric(length(object$betas[[1]]) + 1)
      for(model.ind in 1:object$n_models)
         final_coef <- final_coef + c(object$intercepts[[model.ind]], object$betas[[model.ind]]) * object$models_weights[[model.ind]]
      return(final_coef)
   } else{
      
      if(any(!(group_index %in% 1:object$n_models_optimal)))
         stop("The group index is invalid.")
      
      final_coef <- numeric(length(object$betas[[1]]) + 1)
      for(model.ind in group_index)
         final_coef <- final_coef + c(object$intercepts[[model.ind]], object$betas[[model.ind]]) * object$models_weights[[model.ind]]
      return(final_coef)
   }
}
