## for revision 
## fit our proposal; the tree-based method 
exp.tree.methods <- function(x_train, y_train, x_val, y_val, x_test, y_test, niter, grid, t_range, control.tree.list){
  
  model.tree.list <- list()
  err_trains <- err_tests <-   err_vals <-  matrix(NA, length(control.tree.list), 1000)  # save errors at each iterations 
  #res_trains <- res_tests <- res_vals <- list()
  err_test <- rep(NA, length(control.tree.list)) # one value at early stopping time 
  err_val <- rep(NA, length(control.tree.list)) # one value at early stopping time 
  
  early_stop <- rep(NA, length(control.tree.list)) # save the early stopping time 
  time_vec <- rep(NA, length(control.tree.list))   # save the time required for training 
  
  for(i in 1:length(control.tree.list)){
    
    print(paste(i, "out of", length(control.tree.list), "methods", control.tree.list[[i]]$tree_control$method, "shrinkage", control.tree.list[[i]]$shrinkage, "d",control.tree.list[[i]]$tree_control$d,
                "nknots", control.tree.list[[i]]$nknot, "num_index", control.tree.list[[i]]$tree_control$num_index))
    
    tt <- system.time(model.tree.list[[i]] <- TFBoost(x_train = x_train, 
                                                         y_train = y_train,  x_val = x_val, y_val = y_val,x_test =  x_test, y_test = y_test, grid = grid, t_range = t_range, niter = niter, control = control.tree.list[[i]]))
    
    time_vec[i] <- tt[3]
    err_trains[i,1:length(model.tree.list[[i]]$err_train)] <-model.tree.list[[i]]$err_train
    err_vals[i,1:length(model.tree.list[[i]]$err_val)] <-model.tree.list[[i]]$err_val
    
    err_tests[i,1:length(model.tree.list[[i]]$err_test)] <-model.tree.list[[i]]$err_test
    err_test[i] <- model.tree.list[[i]]$err_test[model.tree.list[[i]]$early_stop]
    err_val[i] <- model.tree.list[[i]]$err_val[model.tree.list[[i]]$early_stop]
    
    early_stop[i] <- model.tree.list[[i]]$early_stop
    #res_trains[[i]] <- model.tree.list[[i]]$f_train_t - y_train # not saving the residuals 
    #res_vals[[i]] <- model.tree.list[[i]]$f_val_t - y_val
    #res_tests[[i]] <- model.tree.list[[i]]$f_test_t - y_test
    model.tree.list[[i]] <- NULL
    
  }
  #return(list(res_trains = res_trains, res_vals = res_vals, res_tests = res_tests, time_vec  = time_vec, err_trains = err_trains, err_tests = err_tests, err_vals = err_vals, err_test = err_test, err_val = err_val, early_stop = early_stop))
  
  return(list(time_vec  = time_vec, err_trains = err_trains, err_tests = err_tests, err_vals = err_vals, err_test = err_test, err_val = err_val, early_stop = early_stop))
}

# perform experiment to run tree or the other competitors 
# case_id = 1 run trees
# case_id = 0 run competitors 

do.exp <- function(seed, g_func_no, SNR, x_type,  niter = 100, niter_FPPR = 10, niters_FAME = c(2,4,6,8,10),  control.tree.list, nknots_FPPR, nknots_FAME, nbasises_FGAM, thetas_FRF, case_id, methods) {
  
  
  dat_gen_control <- dat.generate.control(x_type = x_type, SNR = SNR, g_func_no = g_func_no, n_train = 400, n_val = 200, n_test = 1000)
  dat <- dat.generate(seed = seed, control = dat_gen_control)
  
  x_train <- dat$x$x_train;  x_val <- dat$x$x_val;  x_test <- dat$x$x_test
  y_train <- dat$y$y_train;  y_val <- dat$y$y_val;  y_test <- dat$y$y_test
  
  grid <- dat$tt
  
  dat2return <- NULL
  
  if(x_type == "ferraty"){
    t_range <- c(-1,1)
  }
  if(x_type == "mattern"){
    t_range <- c(0,1)
  }
  
  dat2return <- NULL
  
  if(case_id == 1){
    
    tree.tmp <- exp.tree.methods(x_train, y_train, x_val, y_val, x_test, y_test, niter,  grid, t_range, control.tree.list)
    dat2return <- c(dat2return, list(times_tree = tree.tmp$time_vec, err_trains_tree =tree.tmp$err_trains,
                                     err_tests_tree = tree.tmp$err_tests, err_vals_tree = tree.tmp$err_vals, 
                                     err_test_tree = tree.tmp$err_test, err_val_tree = tree.tmp$err_val, early_stop_tree = tree.tmp$early_stop))        
  }
  
  if(case_id == 0){
    
    if("FLM1" %in% methods){
      dat2return <- perform.FLM1(x_train, y_train, x_val, y_val, x_test, y_test,  t_range, grid, lambdas = seq(0, 2,0.05))
      print( dat2return$lambda.opt)
    }
    if("FLM2" %in% methods){
      dat2return <- c(dat2return, perform.FLM2(x_train, y_train, x_test, y_test, grid))
    }
    if("FAM" %in% methods){
      dat2return <- c(dat2return, perform.FAM(x_train, y_train, x_val, y_val, x_test, y_test,  t_range, grid))
    }
    if("FPPR" %in% methods){
      dat2return <- c(dat2return, perform.FPPR(x_train, y_train, x_val, y_val, x_test, y_test,  t_range, grid, niter_FPPR, nknots_FPPR))
    }
    if("FAME" %in% methods){
      dat2return <- c(dat2return, perform.FAME(x_train, y_train, x_val, y_val, x_test, y_test,  t_range, grid, niters_FAME, nknots_FAME))
    }
    if("FGAM" %in% methods){
      dat2return <- c(dat2return, perform.FGAM(x_train, y_train, x_val, y_val, x_test, y_test,  t_range, grid, nbasises_FGAM))
    }
    if("FRF" %in% methods ){
      dat2return <-  c(dat2return, perform.FRF(x_train, z_train = NULL, y_train,  x_val,  z_val = NULL, y_val, x_test, z_test = NULL, y_test, grid, t_range, m = 500, thetas = thetas_FRF))
    }
    if("RFGroove" %in% methods ){
      dat2return <-  c(dat2return, perform.RFGroove(x_train, z_train = NULL, y_train,  x_val,  z_val = NULL, y_val, x_test, z_test = NULL, y_test, grid, t_range, nbasisInit = 20,  m = 500))
    }
  }   
    dat2return$S <- dat$S
    
   return(dat2return)
}

## FLM1
perform.FLM1 <- function(x_train, y_train, x_val, y_val, x_test, y_test,  t_range, grid, lambdas = seq(0, 2,0.05)){
  
  
  ## get the panalty matrix
  basis_coef <- create.bspline.basis(rangeval = t_range, nbasis = 7)
  R <- getbasispenalty(basis_coef, Lfdobj=2)
  

  dd <- 4; p <- ncol(x_train); nknot  <- 3
  grid0 <- seq(t_range[1],t_range[2], 1/(10*(p-1))) # in case of not evenly spaced
  knot <- quantile(grid0, (1:nknot)/(nknot+1) )
  delta <- sort(c(rep(range(grid0), dd), knot)) #exterior knots
  B<- spline.des(delta, grid, dd)$design
  
  train_predictors <- t(apply(x_train, 1, function(xx){apply(B, 2, function(bb){riemman (xx*bb, grid, t_range)})}))
  val_predictors <- t(apply(x_val, 1, function(xx){apply(B, 2, function(bb){riemman (xx*bb, grid, t_range)})}))
  test_predictors <- t(apply(x_test, 1, function(xx){apply(B, 2, function(bb){riemman (xx*bb, grid, t_range)})}))
  
  # add an intercept 
  Z_tilde_train <- cbind(rep(1,nrow(train_predictors)), train_predictors)
  Z_tilde_val <- cbind(rep(1,nrow(val_predictors)), val_predictors)
  Z_tilde_test <- cbind(rep(1,nrow(test_predictors)), test_predictors)
  
  R_tilde <- rbind(rep(0,ncol(R)), R)
  R_tilde <- cbind(rep(0,nrow(R_tilde)), R_tilde)
  
  mse_vals <- rep(Inf, length(lambdas))
  b_mat <- matrix(NA, length(lambdas), ncol(Z_tilde_train))
  
  for(j in 1:length(lambdas)){
    lambda <- lambdas[j]
    if(is.singular.matrix(t(Z_tilde_train)%*%Z_tilde_train + lambda*R_tilde)){
      b_mat[j,] <- b <- pinv(t(Z_tilde_train)%*%Z_tilde_train + lambda*R_tilde)%*%t(Z_tilde_train)%*%y_train
      
    }else{
      b_mat[j,] <- b <- solve(t(Z_tilde_train)%*%Z_tilde_train + lambda*R_tilde)%*%t(Z_tilde_train)%*%y_train
    }
    tmp_val <- Z_tilde_val%*%b
    mse_vals[j] <- mean(  (tmp_val - y_val)^2)
  }
  
  pred_val <- Z_tilde_val%*%b_mat[which.min(mse_vals),]
  pred_test <- Z_tilde_test%*%b_mat[which.min(mse_vals),]
  pred_train <- Z_tilde_train%*%b_mat[which.min(mse_vals),]
  
  err_val <-  mean( (pred_val - y_val)^2)
  err_test <- mean( (pred_test - y_test)^2)
  err_train <- mean( (pred_train - y_train)^2)
  
  print("FLM1 done!")
  return(list(err_val_FLM1 = err_val, err_test_FLM1 = err_test, err_train_FLM1 = err_train, mse_vals_FLM1 = mse_vals, 
                   lambda.opt  = lambdas[which.min(mse_vals)]))
}


## FLM2
perform.FLM2 <- function(x_train, y_train, x_test, y_test, grid){
  
  n_train <- nrow(x_train)
  n_test <- nrow(x_test)
  
  Lt <- Lx <- newLx <- newLt <- list()
  for (i in 1:n_train) {
    Lt[[i]] <- grid
    Lx[[i]] <- x_train[i,]
  }
  for(i in 1:n_test){
    newLt[[i]]<- grid
    newLx[[i]] <- x_test[i,]
  }
  pca_est <- FPCA(Ly = Lx, Lt = Lt)
  train_predictors <- pca_est$xiEst
  test_predictors <-  predict(pca_est,newLy = newLx, newLt = newLt, xiMethod = "IN",K = ncol(pca_est$xiEst))$scores
  
  dat <- data.frame(x = train_predictors, y  = y_train)
  newdat <- data.frame(x = test_predictors, y  = y_test)
  model_FLM <- lm(y~., data = dat)
  err_test_FLM <- mse(predict(model_FLM, newdat), y_test)
  
  print("FLM2 done!")
  return(list(err_test_FLM2 = err_test_FLM, early_stop_FLM2 = ncol(train_predictors)))

}


## FAM
perform.FAM <- function(x_train, y_train, x_val, y_val, x_test, y_test,  t_range, grid){
    
  Lt <- Lx <- list()
    
    for (i in 1:nrow(x_train)) {
      Lt[[i]] <- grid
      Lx[[i]] <- x_train[i,]
    }
    
    LtTest <- LxTest <- list()
    
    for (i in 1: nrow(x_val)) {
      LtTest[[i]] <- grid
      LxTest[[i]] <- x_val[i,]
    }
    
    for (i in 1: nrow(x_test)) {
      LtTest[[i + nrow(x_val)]] <- grid
      LxTest[[i+ nrow(x_val)]] <- x_test[i,]
    }
    
    
    for (i in 1: nrow(x_train)) {
      LtTest[[i + nrow(x_val) + nrow(x_test)]] <- grid
      LxTest[[i+ nrow(x_val) + nrow(x_test)]] <- x_train[i,]
    }
    
    
    optns <- list(dataType='Dense', error=FALSE, verbose=TRUE)
    
    ## FPCA selects the number of components 
    alpha <- 1
    time_FAM <- system.time(model.FAM <- FAM_update(y_train, Lx, Lt,  newLx = LxTest, newLt = LtTest,alpha = alpha,   optns =optns))[3]
    
    pred_val <- model.FAM$mu +t(apply(model.FAM$fam[1:nrow(x_val), ], 1, sum))
    pred_test <- model.FAM$mu +t(apply(model.FAM$fam[(nrow(x_val)+1):(nrow(x_val)+nrow(x_test)), ], 1, sum))
    pred_train <- model.FAM$mu +t(apply(model.FAM$fam[(nrow(x_val)+nrow(x_test)+1):nrow(model.FAM$fam), ], 1, sum))

    
    list(time_FAM = time_FAM, err_train_FAM = mse(pred_train, y_train), err_val_FAM = mse(pred_val, y_val),err_test_FAM = mse(pred_test, y_test))
         
    print(alpha)
    print("FAM done!")
    return(list(time_FAM = time_FAM, err_train_FAM = mse(pred_train, y_train), err_val_FAM = mse(pred_val, y_val),err_test_FAM = mse(pred_test, y_test)))
}



## FPPR
perform.FPPR <- function(x_train, y_train, x_val, y_val, x_test, y_test,  t_range, grid, niter_FPPR, nknots_FPPR){
  
  time_FPPR <- system.time(model.FPPR <- FPPR(x_train, y_train, x_val, y_val, x_test, y_test, grid = grid, t_range = t_range,  niter = niter_FPPR, nknots = nknots_FPPR))
  err_trains_FPPR <- model.FPPR$err_train
  err_tests_FPPR <- model.FPPR$err_test
  err_vals_FPPR <-  model.FPPR$err_val
  err_val_FPPR  <- model.FPPR$err_val[model.FPPR$early_stop]
  err_test_FPPR <- model.FPPR$err_test[model.FPPR$early_stop]
  Js_FPPR <- model.FPPR$J
  early_stop_FPPR <- model.FPPR$early_stop
  print("FPPR done!")
  return(list(time_FPPR= time_FPPR[3],  Js_FPPR  =   Js_FPPR ,  
                                   err_test_FPPR =  err_test_FPPR, err_val_FPPR =  err_val_FPPR, err_trains_FPPR =  err_trains_FPPR, 
                                   err_tests_FPPR =  err_tests_FPPR, err_vals_FPPR =  err_vals_FPPR, early_stop_FPPR = early_stop_FPPR))

}


## FAME
perform.FAME <- function(x_train, y_train, x_val, y_val, x_test, y_test,  t_range, grid, niters_FAME, nknots_FAME){
  
  time_train <- rep(1:length(grid),length(y_train))
  time_val <- rep(1:length(grid),length(y_val))
  time_test <- rep(1:length(grid),length(y_test))

  curve_train <- unlist(lapply(1:length(y_train), FUN = function(x){rep(x, length(grid))}))
  curve_val <- unlist(lapply(1:length(y_val), FUN = function(x){rep(x, length(grid))}))
  curve_test <- unlist(lapply(1:length(y_test), FUN = function(x){rep(x, length(grid))}))

  data_train <- list(x = c(t(x_train)), y = y_train, time = time_train, curve = curve_train)
  data_val <- list(x = c(t(x_val)), y = y_val, time = time_val, curve = curve_val)
  data_test <- list(x =  c(t(x_test)), y = y_test, time = time_test, curve = curve_test)

  mses_FAME_val <- mses_FAME_test <- rep(Inf,length(niters_FAME))

  best_FMAE_mses <- Inf
  for(j in 1:length(niters_FAME)){
  
    err <- try(time_tmp_FAME <- system.time(model.FAME <- famefn(data_train,tol=.001, q = nknots_FAME+4,p = 4,r = niters_FAME[j], 
                                                               grid= grid, t_range = t_range, maxit=3,family=gaussian(), pert=0,prevfit=NULL, 
                                                               run=5, pc=F))[3])
    if(class(err) == "try-error"){
      break
    }else{
      mse_tmp <-  mse(fame.predict(model.FAME, newdata = data_val, type = "response"),y_val)
      mses_FAME_val[j] <- mse_tmp 
      mses_FAME_test[j] <-mse(fame.predict(model.FAME, newdata = data_test, type = "response"),y_test)
    
      if(mse_tmp < best_FMAE_mses){
        best_FAME_model <-model.FAME
        time_FAME <- time_tmp_FAME
        best_FMAE_mses <- mse_tmp
        early_stop_FAME <- j
      }
    }
  }

  pred_FAME <- fame.predict(best_FAME_model, newdata = data_test, type = "response")
  err_test_FAME <-mse(pred_FAME, y_test) 
  err_train_FAME <-   mse(fame.predict(best_FAME_model, newdata = data_train, type = "response"),y_train)
  print("FAME done!")
  return(list(time_FAME = time_FAME, early_stop_FAME  = early_stop_FAME,
                                 mses_FAME_val = mses_FAME_val, mses_FAME_test = mses_FAME_test,
                                 err_test_FAME = err_test_FAME, err_train_FAME = err_train_FAME))
}



## FGAM
perform.FGAM <- function(x_train, y_train, x_val, y_val, x_test, y_test,  t_range, grid, nbasises_FGAM){
  
  xx <- x_train
  yy <- y_train
  model.FGAM.list <- list()
  time_FGAM_list <- list()
  mses_FGAM <- rep(Inf, length(nbasises_FGAM))
  mses_FGAM_val <- mses_FGAM_test <- rep(Inf,length(nbasises_FGAM))
  
  
  for(j in 1:length(nbasises_FGAM)){   #actually only try 15 
    err <- try(
      time_FGAM_list[[j]] <- system.time(model.FGAM.list[[j]] <- fgam(yy ~ af(xx, splinepars=list(k=c(nbasises_FGAM[j],nbasises_FGAM[j]),m=list(c(2,2),c(2,2)))), gamma = 1.2, method="REML"))[3])
    if(class(err) == "try-error"){
      break
    }else{
      mses_FGAM[j] <- mse(predict(model.FGAM.list[[j]] ,newdata=list(xx =x_val),type='response'), y_val)
      mses_FGAM_val[j] <- mses_FGAM[j]
      mses_FGAM_test[j] <- mse(predict(model.FGAM.list[[j]] ,newdata=list(xx =x_test),type='response'), y_test)
    }
  }
  
  time_FGAM <- time_FGAM_list[[which.min(mses_FGAM)]]
  model.FGAM <- model.FGAM.list[[which.min(mses_FGAM)]]
  Js_FGAM <- nbasises_FGAM[which.min(mses_FGAM)]
  rm(model.FGAM.list)
  err_test_FGAM <-   mse(predict(model.FGAM ,newdata=list(xx =x_test),type='response'), y_test)
  err_train_FGAM <-   mse(predict(model.FGAM ,newdata=list(xx =x_train),type='response'), y_train)
  
  return(list(time_FGAM = time_FGAM, Js_FGAM  =  Js_FGAM , err_test_FGAM  = err_test_FGAM,
                                    err_train_FGAM =  err_train_FGAM, mses_FGAM_val = mses_FGAM_val,
                                    mses_FGAM_test = mses_FGAM_test))
}

FRF <- function(x_train, z_train = NULL, y_train,  x_pred,  z_pred = NULL, y_pred,  grid, t_range, m = 50, theta){
  
  train_pred <- rep(0, nrow(x_train))
  pred_pred <- rep(0, nrow(x_pred))  # the prediction set
  
  
  for(i in 1:m){
    
    if(i%%100 == 0){
      print(paste(i,"th iteration"))
    }
    
    set.seed(i)
    # draw bootstrap sample 
    idx <- sample(1:nrow(x_train), nrow(x_train), replace = TRUE)
    x_train_tree <- x_train[idx,]
    y_train_tree <- y_train[idx]
    
    
    if(!(is.null(z_train))){
      z_train_tree <- z_train[idx,]
    }
    
    r_l <- rexp(100, rate = theta)  # waiting times 
    
    tmp <- cut(grid,  c(t_range[1],t_range[1]+ cumsum(r_l)), include.lowest = TRUE)
    tmp <- as.factor(as.character(tmp))
    tmp_numeric <- as.numeric(tmp)
    
    train_predictors <- pred_predictors <- NULL
    
    for(j in tmp_numeric){
      if(sum(tmp_numeric == j) == 1){
        train_predictors <- cbind(train_predictors, x_train_tree[ ,   tmp_numeric == j])
        pred_predictors <- cbind(pred_predictors, x_pred[ ,   tmp_numeric == j])
      }else{
        train_predictors <- cbind(train_predictors, apply(x_train_tree[ , tmp_numeric == j], 1, mean))
        pred_predictors <- cbind(pred_predictors, apply(x_pred[,  tmp_numeric == j], 1, mean))
      }
    }
    
    if(!(is.null(z_train))){
      train_predictors <- cbind(train_predictors, z_train_tree)
      pred_predictors <- cbind(pred_predictors, z_pred_tree)
    }
    
    dat_train <- data.frame(train_predictors, y = y_train_tree)
    dat_pred <- data.frame(pred_predictors, y = y_pred)
    
    tree.model <- rpart(y~., data = dat_train)
    
    ## training error not useful 
    train_pred <-   train_pred + predict(tree.model, newdata = dat_train)
    pred_pred <-  pred_pred + predict(tree.model, newdata = dat_pred)
    
    #print(c(mse(y_train_tree, train_pred/i), mse(y_train_tree, predict(tree.model, newdata = dat_train)),mse(y_train_tree, mean(y_train_tree))))
    
    #print(c(mse(y_pred, pred_pred/i), mse(y_pred, predict(tree.model, newdata = dat_pred)),mse(y_pred, mean(y_train))))
  }
  
  train_pred <-  train_pred/m
  pred_pred <-  pred_pred/m
  
  return(list(err_train_FPR = mse(train_pred, y_train), err_pred_FPR = mse(pred_pred, y_pred)))
  
}

perform.FRF <- function(x_train, z_train = NULL, y_train,  x_val,  z_val = NULL, y_val, x_test, z_test = NULL, y_test, grid, t_range, m = 50, thetas){
  
  res <- list()
  val_errs <- NULL
  
  aa <- Sys.time()
  for(i in 1: length(thetas)){
    print(paste("testing theta", thetas[i]))
    res[[i]] <- FRF(x_train = x_train, z_train = z_train, y_train = y_train,  x_pred = x_val,  z_pred = z_val, y_pred = y_val,  grid = grid, t_range = t_range, m = m, theta = thetas[i])
    val_errs <- c(val_errs , res[[i]]$err_pred_FPR)
  }
  
  res<- FRF(x_train = x_train, z_train = z_train, y_train = y_train,  x_pred = x_test,  z_pred = z_test, y_pred = y_test,  grid = grid, t_range = t_range, m = m, theta = thetas[which.min(val_errs)])
  
  names(res)[2] <-  "err_test_FRF" 
  bb <- Sys.time()
  
  res$time_FRF <- bb- aa
  res$val_errs <-  val_errs
  return(res)
}

## used for RFGroove
fpca <- function (x, nbasisInit, propVar = 0.9, reconstruct = FALSE, 
                  varName = NULL, verbose = FALSE) {
  interval <- 1:ncol(x)
  nbasis <- ifelse(missing(nbasisInit), ncol(x)/4, nbasisInit)
  if (verbose) 
    cat(nbasis, " Spline basis coefficients\n")
  bsp <- create.bspline.basis(c(1, ncol(x)), nbasis = nbasis)
  basis <- eval.basis(interval, bsp)
  fdObj <- Data2fd(argvals = interval, y = t(x), basisobj = bsp)
  fpca <- pca.fd(fdObj, nharm = nbasis, centerfns = TRUE)
  nrPC <- pmax(which(cumsum(fpca$varprop) >= propVar)[1], 2)
  if (verbose) 
    cat(nrPC, "PCs selected\n")
  optimalDesign <- fpca$scores[, 1:nrPC]
  str <- ifelse(is.null(varName), "PC", paste(varName, "PC", 
                                              sep = "_"))
  colnames(optimalDesign) <- paste(str, 1:nrPC, sep = "")
  if (reconstruct) {
    basisMean <- eval.basis(interval, fpca$meanfd$basis)
    meanFunction <- as.numeric(basisMean %*% fpca$meanfd$coefs)
    smoothData <- t(basis %*% fpca$harmonics$coefs[, 1:nrPC] %*% 
                      t(optimalDesign))
    smoothData <- t(apply(smoothData, MARGIN = 1, FUN = function(z) z + 
                            meanFunction))
    lout <- list(design = optimalDesign, smoothData = smoothData)
  }
  else {
    lout <- list(design = optimalDesign, smoothData = NULL)
  }
  lout
}


# use wave for dimension reduction 
perform.RFGroove <- function(x_train, z_train = NULL, y_train,  x_val,  z_val = NULL, y_val, x_test, z_test = NULL, y_test, grid, t_range, nbasisInit = 4,  m = 500){
  
  
  aa <- Sys.time()
  xx <- rbind(x_train, x_val, x_test)
  verbose <- FALSE
  designMatrixList <- fpca(x = xx, nbasisInit = nbasisInit, verbose = verbose,  propVar = 0.99)$design
  
  
  train_predictors <- designMatrixList[1:nrow(x_train),]
  val_predictors <- designMatrixList[(nrow(x_train)+1):(nrow(x_train)+nrow(x_val)),]
  test_predictors <- designMatrixList[(nrow(x_train)+nrow(x_val)+1): nrow( designMatrixList),]
  
  if(!is.null(z_train)){
    train_predictors <- cbind(train_predictors, z_train)
    val_predictors <- cbind(val_predictors, z_train)
    test_predictors <- cbind(test_predictors, z_train)
  }
  
  forest.model <- randomForest(x = train_predictors,  y = y_train, xtest = test_predictors, ytest = y_test, ntree = m)
  bb <- Sys.time()
  
  
  return(list(time_RFGroove = bb - aa, err_test_RFGroove = mse(forest.model$test$predicted , y_test)))
}



