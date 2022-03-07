# ----------------
# Stacking Models
# ----------------

# Function to return the individualized model predictions
models.predictions <- function(x, y, models, 
                               shrinkage, alpha, n_lambda_sparsity, tolerance, max_iter, n_folds){
  
  # Creating the folds
  folds <- create_folds(nrow(x), nfolds = n_folds)
  # Storing the number of models
  n.models <- length(unlist(models$cv_error))
  # Storing the number of folds
  n.folds <- length(folds)
  
  # Matrix of "leave-one-out predictions" 
  prediction.matrix <- matrix(nrow = nrow(x), ncol = n.models)
  
  if(shrinkage){
    
    # Looping over all models
    for(m in 1:n.models){
      
      # Looping over all folds
      for(k in 1:n.folds){
        
        # Training data
        x.train <- x[-folds[[k]], models$variables[[m]]]
        y.train <- y[-folds[[k]]]
        # Testing data
        x.test <- x[folds[[k]], models$variables[[m]]]
        
        # Training model and getting predictions
        model.fit <- SplitGLM::cv.SplitGLM(x.train, y.train, G = 1, alpha_s = alpha, 
                                           n_lambda_sparsity = n_lambda_sparsity, tolerance = tolerance, max_iter = max_iter, n_folds = n_folds)
        prediction.matrix[folds[[k]], m] <- predict(model.fit, newx = x.test)
      }
    }
  } else{
    
    # Looping over all models
    for(m in 1:n.models){
      
      # Looping over all folds
      for(k in 1:n.folds){
        
        # Training data
        x.train <- x[-folds[[k]], models$variables[[m]]]
        y.train <- y[-folds[[k]]]
        # Testing data
        x.test <- x[folds[[k]], models$variables[[m]]]
        
        # Training model and getting predictions
        model.fit <- SplitGLM::SplitGLM(x.train, y.train, G = 1, alpha_s = alpha, 
                                        lambda_sparsity = 0, lambda_diversity = 0, tolerance = tolerance, max_iter = max_iter)
        prediction.matrix[folds[[k]], m] <- predict(model.fit, newx = x.test)
      }
    }
  }
  
  # Returning the matrix with k-fold predictions
  return(prediction.matrix)
}

# Function to compute the optimal weights (prediction matrix is provided)
model.stacking.matrix <- function(y, prediction.matrix){
  
  # Optimal weights via the prediction matrix from NNLS
  optimal_weights <- nnls::nnls(prediction.matrix, y)$x  
  # Return optimal weights
  return(optimal_weights)
}







