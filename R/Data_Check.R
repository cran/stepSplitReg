# ----------------------------------------------
# Checking Input Data for stepSplitReg Function
# ----------------------------------------------
Data_Check <- function(x, y, n_models, max_variables, keep,
                       model_criterion,
                       stop_criterion, stop_parameter, 
                       shrinkage, alpha, include_intercept, 
                       n_lambda, tolerance, max_iter, n_folds, 
                       model_weights){
  
  # Checking the input for the design matrix (x) and the response vector (y)
  if (all(!inherits(x, "matrix"), !inherits(x, "data.frame"))) {
    stop("x should belong to one of the following classes: matrix, data.frame")
  } else if (all(!inherits(y, "matrix"), all(!inherits(y, "numeric")))) {
    stop("y should belong to one of the following classes: matrix, numeric")
  } else if (any(anyNA(x), any(is.nan(x)), any(is.infinite(x)))) {
    stop("x should not have missing, infinite or nan values")
  } else if (any(anyNA(y), any(is.nan(y)), any(is.infinite(y)))) {
    stop("y should not have missing, infinite or nan values")
  } else {
    if(inherits(y, "matrix")) {
      if (ncol(y)>1){
        stop("y should be a vector")
      }
      # Force to vector if input was a matrix
      y <- as.numeric(y)
    }
    len_y <- length(y)
    if (len_y != nrow(x)) {
      stop("y and x should have the same number of rows")
    }
  }
  
  # Checking the input for the number of models
  if(!is.null(n_models)){
    if (!inherits(n_models, "numeric")) {
      stop("n_models should be numeric")
    } else if (any(!n_models == floor(n_models), n_models <= 1)) {
      stop("n_models should be an integer, greater than one")
    }
  }
  
  # Checking the input for maximum number of variables per model
  if(!is.null(max_variables)){
    if (!inherits(max_variables, "numeric")) {
      stop("max_variables should be numeric")
    } else if (any(!(max_variables == floor(max_variables)), max_variables <= 1, max_variables > nrow(x))) {
      stop("max_variables should be an integer, greater than one and less than the sample size.")
    }
  }
  
  # Check keep value
  if(!inherits(keep, "numeric")) {
    stop("keep should be numeric.")
  } else if(!all(keep <= 1, keep > 0)) {
    stop("keep should be between 0 and 1.")
  }
  
  # Checking input for model criterion 
  if(!(model_criterion %in% c("RSS", "F-test")))
    stop("The shrinkage method must be one of: RSS or F-test.")
  
  # Checking input for stop criterion method
  if(!(stop_criterion %in% c("F-test", "pR2", "aR2", "R2", "Fixed")))
    stop("The stop criterion must be one of: F-test, pR2, aR2, R2 or Fixed.")
  
  # Checking input for stop parameter
  if(stop_criterion %in% c("F-test", "pR2", "aR2", "R2")){
    if(!inherits(stop_parameter, "numeric")) {
      stop("stop_parameter should be numeric.")
    } else if(!all(stop_parameter <= 1, stop_parameter > 0)) {
      stop("stop_parameter should be between 0 and 1.")
    }
  } else{
    if (!inherits(stop_parameter, "numeric")) {
      stop("stop_parameter should be numeric")
    } else if (any(!(stop_parameter == floor(stop_parameter)), stop_parameter <= 1, stop_parameter > nrow(x))) {
      stop("stop_parameter should be an integer, greater than one and less than the sample size.")
    }
  }

  # Check shrinkage parameter
  if(!(shrinkage %in% c(TRUE, FALSE)))
    stop("shrinkage should be TRUE or FALSE.")
  
  # Check alpha value
  if(!inherits(alpha, "numeric")) {
    stop("alpha should be numeric.")
  } else if(!all(alpha <= 1, alpha > 0)) {
    stop("alpha should be between 0 and 1.")
  }
  
  # Check shrinkage parameter
  if(!(include_intercept %in% c(TRUE, FALSE)))
    stop("include_intercept should be TRUE or FALSE.")
  
  # Check input for number of candidates for sparsity value
  if(!inherits(n_lambda, "numeric")) {
    stop("n_lambda should be numeric")
  } else if(any(!n_lambda == floor(n_lambda), n_lambda <= 0)) {
    stop("n_lambda should be a positive integer")
  }
  
  # Check tolerance
  if(!inherits(tolerance, "numeric")) {
    stop("tolerance should be numeric.")
  } else if(!all(tolerance < 1, tolerance > 0)) {
    stop("tolerance should be between 0 and 1.")
  }
  
  # Check maximum number of iterations
  if(!inherits(max_iter, "numeric")) {
    stop("max_iter should be numeric.")
  } else if(any(!max_iter == floor(max_iter), max_iter <= 0)) {
    stop("max_iter should be a positive integer.")
  }
  
  # Check input for number of folds
  if(!inherits(n_folds, "numeric")) {
    stop("n_folds should be numeric")
  } else if(any(!n_folds == floor(n_folds), n_folds <= 0)) {
    stop("n_folds should be a positive integer")
  }
  
  # Checking input for model weights method
  if(!(model_weights %in% c("Equal", "Stacking", "Proportional", "EN")))
    stop("The shrinkage method must be one of: Proportional, Stacking, Equal or EN")
}

# ------------------------------------------------
# Checking Input Data for cv.stepSplitReg Function
# ------------------------------------------------
Data_Check_CV <- function(x, y, n_models, max_variables, keep,
                          model_criterion,
                          stop_criterion, stop_parameter, 
                          shrinkage, alpha, include_intercept, 
                          n_lambda, tolerance, max_iter, n_folds, 
                          model_weights,
                          n_threads){
  
  # Checking the input for the design matrix (x) and the response vector (y)
  if (all(!inherits(x, "matrix"), !inherits(x, "data.frame"))) {
    stop("x should belong to one of the following classes: matrix, data.frame")
  } else if (all(!inherits(y, "matrix"), all(!inherits(y, "numeric")))) {
    stop("y should belong to one of the following classes: matrix, numeric")
  } else if (any(anyNA(x), any(is.nan(x)), any(is.infinite(x)))) {
    stop("x should not have missing, infinite or nan values")
  } else if (any(anyNA(y), any(is.nan(y)), any(is.infinite(y)))) {
    stop("y should not have missing, infinite or nan values")
  } else {
    if(inherits(y, "matrix")) {
      if (ncol(y)>1){
        stop("y should be a vector")
      }
      # Force to vector if input was a matrix
      y <- as.numeric(y)
    }
    len_y <- length(y)
    if (len_y != nrow(x)) {
      stop("y and x should have the same number of rows")
    }
  }
  
  # Checking the input for the number of models
  if(!is.null(n_models)){
    if (!inherits(n_models, "numeric")) {
      stop("n_models should be numeric")
    } else if (any(!n_models == floor(n_models), n_models <= 1)) {
      stop("n_models should be an integer, greater than one")
    }
  }
  
  # Checking the input for maximum number of variables per model
  if(!is.null(max_variables)){
    if (!inherits(max_variables, "numeric")) {
      stop("max_variables should be numeric")
    } else if (any(!(max_variables == floor(max_variables)), max_variables <= 1, max_variables > nrow(x))) {
      stop("max_variables should be an integer, greater than one and less than the sample size.")
    }
  }
  
  # Check keep value
  if(!inherits(keep, "numeric")) {
    stop("keep should be numeric.")
  } else if(!all(keep <= 1, keep > 0)) {
    stop("keep should be between 0 and 1.")
  }
  
  # Checking input for model criterion 
  if(!(model_criterion %in% c("RSS", "F-test")))
    stop("The shrinkage method must be one of: RSS or F-test.")
  
  # Checking input for stop criterion method
  if(!(stop_criterion %in% c("F-test", "pR2", "aR2", "R2", "Fixed")))
    stop("The stop criterion must be one of: F-test, pR2, aR2, R2 or Fixed.")
  
  # Checking input for stop parameter
  if(stop_criterion %in% c("F-test", "pR2", "aR2", "R2")){
    if(!inherits(stop_parameter, "numeric")) {
      stop("stop_parameter should be numeric.")
    } else if(!all(stop_parameter <= 1, stop_parameter > 0)) {
      stop("stop_parameter should be between 0 and 1.")
    }
  } else{
    if (!inherits(stop_parameter, "numeric")) {
      stop("stop_parameter should be numeric")
    } else if (any(!(stop_parameter == floor(stop_parameter)), stop_parameter <= 1, stop_parameter > nrow(x))) {
      stop("stop_parameter should be an integer, greater than one and less than the sample size.")
    }
  }
  
  # Check shrinkage parameter
  if(!(shrinkage %in% c(TRUE, FALSE)))
    stop("shrinkage should be TRUE or FALSE.")
  
  # Check alpha value
  if(!inherits(alpha, "numeric")) {
    stop("alpha should be numeric.")
  } else if(!all(alpha <= 1, alpha > 0)) {
    stop("alpha should be between 0 and 1.")
  }
  
  # Check shrinkage parameter
  if(!(include_intercept %in% c(TRUE, FALSE)))
    stop("include_intercept should be TRUE or FALSE.")
  
  # Check input for number of candidates for sparsity value
  if(!inherits(n_lambda, "numeric")) {
    stop("n_lambda should be numeric")
  } else if(any(!n_lambda == floor(n_lambda), n_lambda <= 0)) {
    stop("n_lambda should be a positive integer")
  }
  
  # Check tolerance
  if(!inherits(tolerance, "numeric")) {
    stop("tolerance should be numeric.")
  } else if(!all(tolerance < 1, tolerance > 0)) {
    stop("tolerance should be between 0 and 1.")
  }
  
  # Check maximum number of iterations
  if(!inherits(max_iter, "numeric")) {
    stop("max_iter should be numeric.")
  } else if(any(!max_iter == floor(max_iter), max_iter <= 0)) {
    stop("max_iter should be a positive integer.")
  }
  
  # Check input for number of folds
  if(!inherits(n_folds, "numeric")) {
    stop("n_folds should be numeric")
  } else if(any(!n_folds == floor(n_folds), n_folds <= 0)) {
    stop("n_folds should be a positive integer")
  }
  
  # Checking input for model weights method
  if(!(model_weights %in% c("Equal", "Stacking", "Proportional", "EN")))
    stop("The shrinkage method must be one of: Proportional, Stacking, Equal or EN")
  
  # Check input for number of threads
  if(!inherits(n_threads, "numeric")) {
    stop("n_threads should be numeric")
  } else if(any(!n_threads == floor(n_threads), n_threads <= 0)) {
    stop("n_threads should be a positive integer")
  }
}