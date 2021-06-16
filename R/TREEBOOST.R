#' Tuning and control parameters for the tree-based functional boosting algorithm
#' 
#' Tuning and control parameters for the TFBoost algorithm.
#'
#' Various tuning and control parameters for the \code{\link{TFBoost}} algorithm implemented in the
#' function 
#' 
#' @param make_prediction logical indicating whether to make predictions using \code{x_test} (defaults to \code{TRUE})
#' @param tree_control control parameters for the tree learner of TFBoost (defaults to \code{TREE.control()}, see \code{\link{TREE.control}})
#' @param loss loss function of TFBoost (character, 'l2' or 'lad', defaults to 'l2')
#' @param user_func list of the loss, gradient of the loss, and initialization function of TFBoost, required when \code{loss} parameter is missing (list, defaults to \code{NULL})
#' @param shrinkage shrinkage parameter in boosting (numeric, defaults to 0.05)
#' @param precision number of significant digits to keep when using validation error to calculate early stopping time (numeric, defaults to 4)
#' @param init_type type of initialization for TFBoost, at the mean or median of the training responses (character, 'mean' or 'median', defaults to 'mean')
#' @param nknot number of interior knots of the cubic B-spline basis (numeric, defaults to 3)
#' @param save_f logical indicating whether to save the function estimates at all iterations (defaults to \code{FALSE})
#' @param trace logical indicating whether to print the number of completed iterations for monitoring progress (defaults to \code{FALSE})
#' @param save_tree logical indicating whether to save the tree objects at all iterations, required when the user calls \code{TFBoost.predict} with the returned object from \code{TFBoost} (defaults to \code{FALSE})
#' @return A list of all input parameters
#'
#' @author Xiaomeng Ju, \email{xiaomeng.ju@stat.ubc.ca}
#' 
#' @export
#' 

TFBoost.control <- function(make_prediction = TRUE, tree_control = TREE.control(), loss = "l2", user_func = NULL, shrinkage  = 0.05, precision = 4, 
                               init_type = "mean", nknot = 3, save_f = FALSE, trace = FALSE, save_tree = FALSE){
  
  return(list(make_prediction =  make_prediction, tree_control = tree_control, loss = loss, user_func = user_func, shrinkage = shrinkage, precision = precision, 
              init_type = init_type, nknot = nknot, save_f = save_f, trace = trace, save_tree = save_tree))
}

#' Tree-based functional boosting 
#' 
#' This function implements a tree-based boosting algorithm for functional regression,
#'
#' This function implements  a tree-based boosting algorithm for functional regression (TFBoost).
#' This function uses the functions available in the \code{rpart} package to construct functional multi-index regression trees.
#' 
#' @param x_train functional predictor matrix in the training data (matrix/dataframe)
#' @param z_train scalar predictor matrix in the training data (matrix/dataframe, optional)
#' @param y_train scalar response vector in the training data (vector/dataframe)
#' @param x_val functional predictor matrix in the validation data (matrix/dataframe)
#' @param z_val scalar predictor matrix in the validation data (matrix/dataframe, optonal)
#' @param y_val scalar response vectorin the validation data (vector/dataframe)
#' @param x_test functional predictor matrix for test data (matrix/dataframe, optional, required when \code{make_prediction} in control is \code{TRUE})
#' @param z_test scalar predictor matrix for test data (matrix/dataframe, optional, required when \code{make_prediction} in control is \code{TRUE})
#' @param y_test scalar response vector for test data (vector/dataframe, optional, required when \code{make_prediction} in control is \code{TRUE} and \code{z_train} and \code{z_val} are provided)
#' @param grid common grid that the \code{x_train}, \code{x_val}, and \code{x_test} are evaluated at (vector)
#' @param t_range domain of the functional predictor, provided in the form of c(left end of the domain, right end of the domain) (vector)
#' @param niter number of boosting iterations (numeric)
#' @param control a named list of control parameters, as returned by \code{\link{TFBoost.control}}
#' 
#' @return A list with the following components:
#'
#' \item{B}{predicted values with model at the early stopping iteration using x_test (or x_test and z_test) as the predictors}
#' \item{loss_train}{a vector of training errors for all iterations}
#' \item{loss_val}{a vector of validation errors for all iterations}
#' \item{f_train_t}{predicted values with model at the early stopping iteration on the training predictors}
#' \item{f_val_t}{predicted values with model at the early stopping iteration on the validation predictors}
#' \item{f_t_test}{predicted values with model at the early stopping iteration on the test predictors (returned if make_prediction = TRUE in control)}
#' \item{early_stop}{early stopping iteration}
#' \item{err_train}{a vector of training mean-squared-errors for all iterations}
#' \item{err_val}{a vector of validation mean-squared-errors for all iterations}
#' \item{err_test}{a vector of test mean-squared-errors before and at the early stopping iteration (returned if make_prediction = TRUE in control)}
#' \item{grid}{\code{grid} form the input arguments}
#' \item{t_range}{\code{t_range} from the input arguments}
#' \item{init_vals}{a constant to initialize the function estimates for training, validation, and test sets}
#' \item{tree_obj}{a list of fitted functional multi-index trees (one per iteration and returned from \link{TREE})}
#' \item{alpha}{a vector of boosting step sizes (one per iteration)}
#' \item{control}{\code{control} from the input arguments}
#' \item{save_f_train}{a matrix of training function estimates at all iterations (returned if save_f = TRUE in control)}
#' \item{save_f_val}{a matrix of validation function estimates at all iterations (returned if save_f = TRUE in control)}
#' \item{save_f_test}{a matrix of test function estimates before and at the early stopping iteration (returned if save_f = TRUE and make_prediction = TRUE in control)}
#'
#' @author Xiaomeng Ju, \email{xiaomeng.ju@stat.ubc.ca}
#' @export

TFBoost <- function(x_train, z_train = NULL, y_train,  x_val,  z_val = NULL, y_val, x_test, z_test = NULL, y_test, grid, t_range, niter = 10, control = TFBoost.control()){
  
  my.envir <- list2env(control)
  control_names <- names(control)
  for(g in  control_names) {
    assign(g, get(g, envir=my.envir))
  }
  
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    oldseed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", oldseed, envir = .GlobalEnv))
  }
  
  n_train <- nrow(x_train)
  n_val <- nrow(x_val)
  
  if(save_f){
    save_f_train <- matrix(NA, n_train, niter )
    save_f_val <- matrix(NA, n_val, niter )
  }
  
  if(missing(t_range)){
     t_range <- c(min(grid), max(grid))
  }
  
#  if(!is.null(z_train) & is.matrix(z_train)){
#     z_train <- data.frame(z_train)
#  }
  
  dd <- 4; p <- ncol(x_train)
  grid0 <- seq(t_range[1],t_range[2], 1/(10*(p-1))) # in case of not evenly spaced
  knot <- quantile(grid0, (1:nknot)/(nknot+1) )
  delta <- sort(c(rep(range(grid0), dd), knot)) #exterior knots
  B<- spline.des(delta, grid, dd)$design
  B <- compute.orthonormal(B,grid, t_range)
  train_predictors <- t(apply(x_train, 1, function(xx){apply(B, 2, function(bb){riemman (xx*bb, grid, t_range)})}))
  val_predictors <- t(apply(x_val, 1, function(xx){apply(B, 2, function(bb){riemman (xx*bb, grid, t_range)})}))
  
  if(missing(x_test)) {
    make_prediction <- FALSE 
  } 
  
  tree.obj <- list()
  alpha <- rep(NA, niter)
  
  loss_train <- loss_val <- rep(NA, niter)
  err_train <- err_val <- rep(NA, niter)
  
  # initialize functions
  if(!missing(loss)){
    init_tmp <- init.boosting(loss, init_type) 
    init.func <- init_tmp$init.func
  }else{
    init_tmp <- user_func
  }

  func.grad <- init_tmp$func.grad
  func <- init_tmp$func
  
  f_train_t <- f_val_t <-  init_vals <- init.func(y_train)

  for(i in 1:niter){
    
    if(save_f){
      save_f_train[,i] <- f_train_t
      save_f_val[,i] <- f_val_t
    }
    
    if(trace){
      if(i%%100 ==0 )
      print(paste(i,"th iteration out of", niter, "iterations"))
    }

    u <- as.numeric(cal.neggrad(y_train, f_train_t, func.grad))
    tree.obj[[i]] <- TREE(x = train_predictors, y = u, z = z_train, newx = val_predictors, newz = z_val, random.seed = i, control = tree_control)
    h_train_t <-tree.obj[[i]]$pred_train
    h_val_t <- tree.obj[[i]]$pred_test

    alpha[i] <- cal.alpha(y_train = y_train, f_train_t = f_train_t, h_train_t = h_train_t, func = func, loss = loss)
    f_train_t <-  f_train_t + shrinkage*alpha[i]*h_train_t
    f_val_t <-  f_val_t + shrinkage*alpha[i]*h_val_t
    err_train[i] <- mse(f_train_t, y_train)
    err_val[i] <- mse(f_val_t, y_val)
    
    loss_train[i] <- mean(func(y_train -  f_train_t))
    loss_val[i] <- mean(func(y_val - f_val_t))

    if(i == 1){
      early_stop <- 1
      f_train_early <- f_train_t
      f_val_early <- f_val_t
    }else{
      if(round(loss_val[i], precision) < min(round(loss_val[1:(i-1)], precision))){
          early_stop <- i
          f_train_early <- f_train_t
          f_val_early <- f_val_t
      }   
    }
  }
  f_train_t <- f_train_early
  f_val_t <- f_val_early
  
  model <- list(B = B, loss_train=loss_train, loss_val = loss_val, f_train_t =  f_train_t, f_val_t = f_val_t, early_stop = early_stop, err_train = err_train, 
                 err_val = err_val, grid = grid,  t_range = t_range, init_vals = init_vals,  tree.obj = tree.obj, alpha = alpha, control = control)
  
  if(make_prediction){
    tmp_predict <- TFBoost.predict(model, newx = x_test, newy = y_test, newz = z_test)
    model <- c(model, list(f_test_t = tmp_predict$pred, err_test = tmp_predict$err_test))
    if(save_f){
      model <- c(model, list(save_f_test = tmp_predict$save_f_test))
    }
  }
      
  if(save_f){
    model <- c(model, list(save_f_train = save_f_train, save_f_val = save_f_val))
  }
  if(!save_tree){
    model$tree.obj <- NULL
  }
  return(model)
}
    
#' TFBoost.predict
#'
#' A function to make predictions and calculate test error given an object returned by TFBoost and test data
#'
#' A function to make predictions and calculate test error given an object returned by TFBoost and test data
#'
#'@param model an object returned by TFBoost
#'@param newx functional predictor matrix for test data (matrix/dataframe)
#'@param newy scalar predictor matrix for test data (matrix/dataframe, optional, requires if model is trained with z_train and z_val)
#'@param newz scalar response vector for test data (vector/dataframe)
#'@return A list with with the following components:
#'
#' \item{f_t_test}{predicted values with model at the early stopping iteration using x_test (or x_test and z_test) as the predictors}
#' \item{err_test}{a vector of test errors before and at the early stopping iteration (returned if newy is provided)}
#' \item{save_f_test}{a matrix of test function estimates at all iterations (returned if save_f = TRUE in control)}
#'
#' @author Xiaomeng Ju, \email{xmengju@stat.ubc.ca}
#' 
#' @export

TFBoost.predict <- function(model, newx, newy, newz = NULL){

  my.envir <- list2env(model$control)
  control_names <- names(model$control)
  for(g in  control_names) {
    assign(g, get(g, envir=my.envir))
  }
  
  my.envir <- list2env(model)
  control_names <- names(model)
  for(g in  control_names) {
    assign(g, get(g, envir=my.envir))
  }
 
  if(save_f){
    save_f_test <- matrix(NA, nrow(newx), early_stop)
  }

  test_predictors <- t(apply(newx, 1, function(xx){apply(B, 2, function(bb){riemman (xx*bb, grid, t_range)})}))
  f_test_t  <- init_vals
  err_test <- rep(NA,early_stop)
  
  for(i in 1: early_stop){
    if(save_f){
      save_f_test[,i] <- f_test_t
    }
    f_test_t <- f_test_t + shrinkage*alpha[i]* TREE.predict(tree.obj[[i]], newx =  test_predictors, newz = newz)$pred
    if(!missing(newy)){
      err_test[i] <- mse(f_test_t, newy)
    }
  }
  
  if( (missing(newy)) & (save_f == FALSE)){
    return(f_test_t)
  }else{
    res <- list(pred = f_test_t)
    if(!missing(newy)){
      res <- c(res, list(err_test = err_test))
    }
    if(save_f){
      res <- c(res,  list(save_f_test =  save_f_test))
    }
    return(res)
  }
}


# transform a basis matrix to an orthonormal basis matrix
compute.orthonormal <- function(B, grid, t_range){
  
  d <- ncol(B)
  Phi_i <- B
  Psi <- matrix(NA, nrow = nrow(B), ncol(B))
  Psi[,1] <-   B[,1]/ sqrt(riemman(B[,1]*B[,1], grid, t_range))
  
  for(i in 2:d){
    for(j in 1:(i-1)){
      Phi_i[,i] <-   Phi_i[,i]  -  riemman(Phi_i[,i]*Psi[,j], grid, t_range) * Psi[,j]
    }
    Psi[,i] <-   Phi_i[,i]/ sqrt(riemman(Phi_i[,i]*Phi_i[,i], grid, t_range))
  }
  return(Psi)
}

