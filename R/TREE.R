#' Tuning and control parameters for the functional multi-index tree
#' 
#' Tuning and control parameters for the functional multi-index tree. 
#'
#' Various tuning and control parameters for the functional multi-index tree algorithm implemented in the
#' function \code{\link{TREE}}
#' 
#' @param make_prediction logical indicating whether to make predictions using \code{newx} (defaults to \code{TRUE})
#' @param d the maximum depth of the functional multi-index tree (numeric, defaults to 1)
#' @param minbucket the minimum number of observations per node of the functional multi-index tree (numeric, defaults to 2)
#' @param tree_type type of the functional multi-index tree, 'A' or 'B' (character, defaults to 'A')
#' @param num_index number of indices for type a tree (numeric, defaults to 1)
#' @param nmulti number of random initial points  for performing optimization to fit a type A tree (numeric, defaults to 5)
#' @param nscreen number of random points from which to select the initial points to fit a type A tree (numeric, defaults to 30)
#' @param num_dir number of random directions sampled at every iteration for type b tree (numeric, defaults to 200)
#' 
#' @return A list of all input parameters
#'
#' @author Xiaomeng Ju, \email{xiaomeng.ju@stat.ubc.ca}
#' 
#' @export
#' 
TREE.control <- function(make_prediction = TRUE, d = 1, minbucket = 2, tree_type = "A",  nmulti= 3, nscreen = 50, num_dir = 20, num_index = 1){
  return(list(make_prediction = make_prediction, d = d, minbucket = minbucket, tree_type = tree_type, num_index = num_index, nmulti = nmulti, nscreen = nscreen,
              num_dir = num_dir))
}

#' Functional multi-index tree 
#' 
#' This function implements a algorithm for functional multi-index tree.
#'
#' This function implements a functional multi-index tree algorithm devloped based on functions available in the \code{rpart} package. 
#' 
#' @param x matrix of the functional predictor's basis projections in the training data (matrix/dataframe)
#' @param z matrix of the scalar predictors in the training data (matrix/dataframe)
#' @param y response vector in the training data (vector/dataframe)
#' @param newx  matrix of the functional predictor's basis projections in the test data (matrix/dataframe, optional)
#' @param newz  matrix of the scalar predictors in the test data (matrix/dataframe, optional)
#' @param newy response vector in the test data (matrix/dataframe, optional)
#' @param random.seed a seed that controls the randomness of type B tree (numerical)
#' @param control a named list of control parameters, as returned by \code{\link{TREE.control}}
#' 
#' @return A list with the following components:
#'
#' \item{beta_opt}{estimated index coefficients for type A tree, returned if \code{tree_type = "A"}}
#' \item{betas_selected}{the selected random directions for type B tree, returned if \code{tree_type = "B"}}
#' \item{pred_train}{predicted values with model on the training predictors}
#' \item{pred_test}{predicted values with model on the test predictors, returned if \code{make_predictions = TRUE}}
#' \item{control}{\code{control} from the input arguments}
#' \item{tree.model}{fitted tree estimator, an rpart object}
#' @author Xiaomeng Ju, \email{xiaomeng.ju@stat.ubc.ca}
#' 
#' @export

TREE <-function(x, y, z, newx, newy, newz, random.seed, control = TREE.control()) {
  
  my.envir <- list2env(control)
  control_names <- names(control)
  for(g in  control_names) {
    assign(g, get(g, envir=my.envir))
  }
  
  if(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    oldseed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", oldseed, envir = .GlobalEnv))
  }
  
  p <- ncol(x)

  if(tree_type == "A"){
    
    fval.min <- Inf
    set.seed(random.seed)
    optim.fn <-  function(param){tree.loss(param, x, y, z, d,  minbucket, num_index , type = "theta")}
    start_save <- screen.dir(nmulti, nscreen, p, num_index, optim.fn, max_iter = 10)$init_points

    for(i in 1:nmulti) {
      
      param <- start_save[((i-1)*num_index +1):(i*num_index),] 
      optim.return <- tryCatch(Nelder_Mead(fn = optim.fn, par= as.numeric(param), lower = c(rep(-pi/2,num_index), rep(0, (p-2)*num_index)), upper = c(rep(pi/2,num_index), rep(pi, (p-2)*num_index)), control = list(maxfun = 1000)), warning = function(a){return(1)})
      
      while(class(optim.return) == "numeric"){
        print("resample!")
        param <- start_save[((i-1)*num_index +1):(i*num_index),] <-  sample.dir(1, p, num_index, type = "theta")
        optim.return <- tryCatch(Nelder_Mead(fn = optim.fn, par = as.numeric(param), lower = c(rep(-pi/2,num_index), rep(0, (p-2)*num_index)), upper = c(rep(pi/2,num_index), rep(pi, (p-2)*num_index)), control = list(maxfun = 1000)), warning = function(a){return(1)}) 
      }
      
      par_tmp <- matrix(optim.return$par, ncol = p-1)
      
      if(optim.return$fval < fval.min) {
          if(num_index == 1){
            beta_opt <-  as.matrix(p_to_c(par_tmp))
          }else{
            beta_opt <- apply(par_tmp, 1, p_to_c) #ncol: number of index, nrow: p
          }
          fval.min <- optim.return$fval
        }
      }
    
    if(is.null(z)){
      index <- x %*% beta_opt
    }else{
      index <- cbind(x %*% beta_opt, z)
      colnames(index) <-  c(paste("X", 1:num_index, sep = ""), colnames(z))
    }
    dat <- data.frame(index, y = y)
    
    tree.model <- rpart(y~., data = dat, control = rpart.control(maxdepth = d, cp = 0, minbucket =  minbucket))
    model <- list(pred_train = predict(tree.model), beta_opt = beta_opt, control = control, tree.model = tree.model)
    
    if(make_prediction){
      if(!missing(newx)){
        tmp <-  TREE.predict(model, newx = newx, newz = newz)
        model$pred_test <- tmp$pred
        if(!missing(newy)){
          model$test_mse <- mean( (model$pred_test  - newy)^2)
        }
      }
    }
  }
  
  if(tree_type == "B"){
    
    set.seed(random.seed)
    betas <-  t(sample.dir(num_dir, p, num_index = 1, type = "beta"))
    
    if(is.null(z)){
      index <- x %*% betas
    }else{
      index <- cbind(x %*% betas, z)
      index <- matrix(index, ncol = num_dir + ncol(z))
      colnames(index) <-  c(paste("X", 1:num_dir, sep = ""), colnames(z))
    }
    dat <- data.frame(index, y = y)
    
    tree.model <- rpart(y~., data = dat, control = rpart.control(maxdepth = d,  cp = 0, minbucket =  minbucket))
    betas <-data.frame(betas)
    used.var <- setdiff(tree.model$frame$var, "<leaf>")
    used.var <- setdiff(used.var,colnames(z))
    
    if(length(used.var)==1){
      betas_selected <-data.frame(betas[, used.var])
      names(betas_selected) <- used.var 
    }else{
      betas_selected <-betas[, used.var]
    }
    
    model <- list(betas_selected =  betas_selected, pred_train = predict(tree.model), control = control, tree.model = tree.model)
    
    if(make_prediction){
      if(!missing(newx)){
        tmp <-  TREE.predict(model, newx = newx, newz = newz)
        model$pred_test <- tmp$pred
        if(!missing(newy)){
          model$test_mse <- mean((model$pred_test  - newy)^2)
        }     
      }
    }
  }
  return(model)
}

#' TREE.predict
#'
#' A function to make predictions given an object returned by TREE and test data
#'
#' A function to make predictions given an object returned by TREE and test data
#'
#'@param model an object returned by TREE
#'@param newx  matrix of the functional predictor's basis projections in the test data (matrix/dataframe)
#'@param newz  matrix of the scalar predictors in the test data (matrix/dataframe, required if the model includes scalar predictors) 
#'@return 
#'
#' \item{pred}{predicted values with model using x_test (or x_test and z_test) as the predictors}
#'
#' @author Xiaomeng Ju, \email{xmengju@stat.ubc.ca}
#' 
#' @export
#' 
TREE.predict <- function(model, newx, newz){

    if(model$control$tree_type == "A"){
      if(is.null(newz)){
        index <- as.matrix(newx %*% model$beta_opt)
      }else{
        if(!is.matrix(newz)){
          newz = as.matrix(newz, dimnames = list(NULL, names(newz)))
        }
        index <- cbind(as.matrix(newx %*% model$beta_opt), newz)
        colnames(index)[1:ncol(model$beta_opt)] <- c(paste("X", 1:ncol(model$beta_opt), sep = ""))
      }
    }

    if(model$control$tree_type %in% c("B") ){
      betas <- data.frame(matrix(0, ncol(newx), model$control$num_dir))
      betas[, names(model$betas_selected)] <- model$betas_selected
      if(is.null(newz)){
        index <- newx %*% as.matrix(betas)
      }else{
        if(!is.matrix(newz)){
          newz = as.matrix(newz, dimnames = list(NULL, names(newz)))
        }
        index <- cbind(newx %*% as.matrix(betas), newz)
        colnames(index)[1:model$control$num_dir] <- c(paste("X", 1:model$control$num_dir, sep = ""))
      }
    }
   
    dat <- data.frame(index)
    pred <- predict(model$tree.model, newdata = dat)
    return(list(pred = pred))
}


p_to_c <- function(theta){
  
  d <- length(theta)+1
  beta <- rep(NA, d)
  beta[1] <-cos(theta[1])
  tmp1 <- cumprod(sin(theta))
  
  if(d>2){
    for(i in 2:(d-1)){
      beta[i] <- cos(theta[i])*tmp1[i-1]
    }
  }
  beta[d] <- tmp1[d-1]
  
  return(beta)
}


c_to_p <- function(beta){
  
  d <- length(beta)
  theta <- rep(NA, d-1)
  tmp <- beta[d]
  
  for(i in (d-1):1){
    theta[i] <- atan(tmp/beta[i])
    if(theta[i] <0 & i> 1){
      theta[i] <-  theta[i] + pi
    }
    tmp <- tmp/sin(theta[i])
  }
  
  return(theta)
}



tree.loss<- function(param, x, y, z, d, minbucket, num_index, type = "theta") {
  
  param <- matrix(param, nrow = num_index) # every row is an index 
  
  if(type == "theta"){
    if(num_index == 1){
      beta <- p_to_c(param)
    }else{
      beta <- apply(param, 1,p_to_c) # apply to rows
    }
  }else{
    beta <- param
  }
  
  loss <- function(x, y, z,beta) {
    if(!is.null(z)){
      index <- cbind(x %*% beta, z)
    }else{
      index <- x %*% beta
    }
    dat <- data.frame(index, y = y)
    tree.model <- rpart(y~.,  data = dat, control = rpart.control(maxdepth = d,cp = 0, minbucket = minbucket))
    t.ret <- mean((y-predict(tree.model))^2)
    return(t.ret)
  }
  return(loss(x,y,z,beta))
}

sample.dir <- function(num_sample, p, num_index, type = "beta"){
  
  if(type == "theta"){
    tmppp <- matrix(runif(num_sample*(p-1)*num_index, 0,pi), nrow = num_sample*num_index)
    tmppp[,1] <- tmppp[,1] - pi/2
  }else{
    tmp <- matrix(rnorm(num_sample*p*num_index, 0,1), nrow = num_sample*num_index)
    tmpp <- apply(tmp, 1, function(x) { sqrt(sum(x^2))})
    tmppp <- apply(tmp, 2, function(x){x/tmpp})
    tmppp[,1] <- abs(tmppp[,1])
  }
  pars <- tmppp
  return(pars)
}

screen.dir <- function(nmulti, nscreen, p, num_index, optim.fn,  max_iter = 10){
  
  tmp <-  as.matrix(sample.dir(nscreen, p, num_index, type = "theta"))
  
  if(p == 2|| (num_index == 1 & nscreen == 1)){
    tmp <- t(tmp)
  }
  
  cal_nelder_mead <- function(i){
    par_val <- as.numeric(tmp[((i-1)*num_index +1):(i*num_index),])
    suppressWarnings(tmp.return <-Nelder_Mead(fn = optim.fn, par = par_val , lower = c(rep(-pi/2,num_index), rep(0, (p-2)*num_index)), upper = c(rep(pi/2,num_index), rep(pi, (p-2)*num_index)), control = list(maxfun = max_iter)))
    return(tmp.return) 
  }
  
  list_tmp <- lapply(1:nscreen, cal_nelder_mead)
  loss_sampled <- unlist(lapply(list_tmp, function(x){x$fval}))
  idx <- which(rank( loss_sampled, ties.method = "first") <= nmulti)
  idx_tmp  <- as.numeric(sapply(idx, function(x){((x-1)*num_index+1):(x*num_index)}))
  
  return(list(init_points = matrix(tmp[idx_tmp,], ncol = p-1)))
}


# set.seed(123)
# n <- 200; p <- 5
# x<- matrix(runif(n*p), ncol = p)
# y <- x[,1] + x[,2] + rnorm(n)
# train_idx <- sample(n,50)
# test_idx <- setdiff(1:n, train_idx)
# model <- TREE(x, y, newx = x[1:10,], nmulti = 3, random.seed = 42, control = TREE.control(d = 3, num_dir = 1000, method = "method3"))
# model$test_mse
