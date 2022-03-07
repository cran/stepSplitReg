#' 
#' @title Cross Validation - Stepwise Split Regularized Regression
#' 
#' @description \code{cv.stepSplitReg} performs the CV procedure for stepwise split regularized regression.
#' 
#' @param x Design matrix.
#' @param y Response vector.
#' @param n_models Number of models into which the variables are split.
#' @param max_variables Maximum number of variables that a model can contain.
#' @param keep Proportion of models to keep based on their individual cross-validated errors. Default is 1.
#' @param model_criterion Criterion for adding a variable to a model. Must be one of c("F-test", "RSS"). Default is "F-test".
#' @param stop_criterion Criterion for determining when a model is saturated. Must be one of c("F-test", "pR2", "aR2", "R2", "Fixed"). Default is "F-test".
#' @param stop_parameter Parameter value for the stopping criterion. Default is 0.05 for "F-test".
#' @param shrinkage TRUE or FALSE parameter for shrinkage of the final models. Default is TRUE.
#' @param alpha Elastic net mixing parmeter for model shrinkage. Default is 3/4.
#' @param include_intercept TRUE or FALSE parameter for the inclusion of an intercept term.
#' @param n_lambda Number of candidates for the sparsity penalty parameter. Default is 100.
#' @param tolerance Convergence criteria for the coefficients. Default is 1e-3.
#' @param max_iter Maximum number of iterations in the algorithm. Default is 1e5.
#' @param n_folds Number of cross-validation folds. Default is 10.
#' @param model_weights Criterion to determine the weights of the model for prediciton. Must be one of c("Equal", "Proportional", "Stacking"). Default is "Equal".
#' @param n_threads Number of threads. Default is 1.
#' 
#' @return An object of class cv.stepSplitReg.
#' 
#' @export
#' 
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#' 
#' @seealso \code{\link{coef.cv.stepSplitReg}}, \code{\link{predict.cv.stepSplitReg}}
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
cv.stepSplitReg <- function(x, y, n_models = NULL, max_variables = NULL, keep = 1,
                            model_criterion = c("F-test", "RSS")[1],
                            stop_criterion = c("F-test", "pR2", "aR2", "R2", "Fixed")[1], stop_parameter = 0.05, 
                            shrinkage = TRUE, alpha = 3/4, include_intercept = TRUE, 
                            n_lambda = 100, tolerance = 1e-3, max_iter = 1e5, n_folds = 10, 
                            model_weights = c("Equal", "Proportional", "Stacking")[1], 
                            n_threads = 1){
  
  # Check function input
  Data_Check_CV(x, y, n_models, max_variables, keep,
                model_criterion,
                stop_criterion, stop_parameter, 
                shrinkage, alpha, include_intercept, 
                n_lambda, tolerance, max_iter, n_folds, 
                model_weights,
                n_threads)
  
  # Setting the numerical index for the stop criterion
  if(stop_criterion=="F-test")
    stop_criterion <- 4 else if(stop_criterion=="pR2")
      stop_criterion <- 3 else if(stop_criterion=="aR2")
        stop_criterion <- 2 else if(stop_criterion=="R2")
          stop_criterion <- 1 else if(stop_criterion=="Fixed")
            stop_criterion <- 0
          
  # Setting the numerical index for the stop criterion
  if(model_criterion=="RSS")
    model_criterion <- 1 else if(model_criterion=="F-test")
      model_criterion <- 2 
  
  # Shuffle the data
  n <- nrow(x)
  p <- ncol(x)
  random.permutation <- sample(1:n, n)
  x.permutation <- x[random.permutation, ]
  y.permutation <- y[random.permutation]
  
  # Setting the default number of models
  if(is.null(n_models))
    n_models <- floor(2*sqrt(p))
  # Setting the default maximum number of variables per model
  if(is.null(max_variables))
    max_variables <- n
  
  # CPP input formatting
  shrinkage.cpp <- sum(shrinkage)
  include_intercept.cpp <- sum(include_intercept)
  
  # Stepwise split algorithm - CV output
  output <-  CV_Stepwise_Split(x = x.permutation, y = y.permutation,
                               n_models = n_models, max_variables_per_model = max_variables,
                               model_criterion = model_criterion,
                               stop_criterion = stop_criterion, stop_parameter = stop_parameter,
                               shrinkage = shrinkage.cpp, alpha = alpha, include_intercept = include_intercept.cpp,
                               n_lambda = n_lambda, tolerance = tolerance, max_iter = max_iter,
                               n_folds = n_folds,
                               n_threads = n_threads)
  output$n_models <- n_models
  n_models_optimal <- length(output$variables)
  
  # Removing incomplete models
  incomplete.models <- which(sapply(output$variables, length, simplify=TRUE)<1)
  if(length(incomplete.models)>=1){
    output$variables[incomplete.models] <- rep(NULL, length(incomplete.models))
    output$betas[incomplete.models] <- rep(NULL, length(incomplete.models))
    # Adjusting the number of models
    n_models_optimal <- n_models_optimal - length(incomplete.models)
  }
  
  # Eliminating models
  models.keep <- order(unlist(output$cv_error), decreasing=FALSE)[1:round(keep*n_models_optimal,0)]
  output$variables <- output$variables[models.keep]
  output$betas <- output$betas[models.keep]
  output$cv_error <- output$cv_error[models.keep]
  output$n_models_optimal <- length(models.keep)
  
  # Computing the final betas output
  for(k in 1:output$n_models_optimal){
    beta.temp <- rep(0, p)
    beta.temp[output$variables[[k]]] <- as.numeric(output$betas[[k]])
    output$betas[[k]] <- beta.temp
  }
  
  # Model weights methods
  if (model_weights == "Stacking"){ # Model Stacking
    
    prediction.matrix <- models.predictions(x.permutation, y.permutation, output,
                                            shrinkage, alpha, n_lambda, tolerance, max_iter, n_folds)
    output$models_weights <- model.stacking.matrix(y.permutation, prediction.matrix)
    
  } else if (model_weights == "Proportional"){ # Proportional weights
    
    output$models_weights <- (1/unlist(output$cv_error))/sum(1/unlist(output$cv_error))
    
  } else if (model_weights == "Equal"){ # Equal weights
    
    output$models_weights <- rep(1/output$n_models_optimal, output$n_models_optimal)
    
  }
  
  # Create the object of class "cv.stepSplitReg"
  class(output) <- append("cv.stepSplitReg", class(output))
  
  # returning the output
  return(output)
}

