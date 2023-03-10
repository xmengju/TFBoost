# Programma di prova
# Versione settembre 2019

###############################################################
###############################################################

library(stats)
library(splines)

NZD <- function(a) {
  return(a)
  #ifelse(a<0,pmin(-.Machine$double.eps,a),pmax(.Machine$double.eps,a))
}

###################################################################
# Funzione obiettivo da minimizzare
# Si utilizza l'ottimizzatore OPTIM
###################################################################

fun.ob <- function(arg, Y, Z, B, bw)
  {
     arg <- as.matrix(arg)
     arg <- arg / sqrt( mean( (B %*% arg)^2 ) )
     X.arg <- Z %*% arg
     mean( (Y - ker.reg( as.vector(X.arg), as.vector(Y), as.vector(X.arg), bw)$hat.y )^2 )
  }


###################################################################
# Stima a nucleo (Nadaraya-Watson)
###################################################################

# Algoritmo con finestra globale fissata
ker.reg <- function(x, y, x.new, bw = NULL, loo = FALSE, kernel_type = "Gaussian")
 {
     x <- as.vector(x)
     x.new <- as.vector(x.new)
     y <- as.vector(y)
     if(is.null(bw)) bw <- diff(range(x)) * 0.2   # bandwidth = 20% range(x)

    # Estimation
     dist <- abs(outer(x, x, "-")) 
     u <- dist / bw 
     
     if(kernel_type == "Epanechnikov"){
       K.val <- 3/4 * (1 - u^2)       # kernel di Epanechnikov
       K.val[ K.val < 0 ] <- 0
       K.val[ K.val > 1 ] <- 0
     }
     if(kernel_type =="Gaussian"){
       K.val <- dnorm(u)
     }
     diag(K.val) <- 0               # elimino la i-ma coppia di osservazioni per la CV
     denom <- apply( K.val, 2, sum)   
     hat.y <- K.val %*% as.matrix(y) / NZD(denom)
     
     if(loo){
       diag(K.val)  = 0
       hat.y.loo <- t(K.val) %*% as.matrix(y) / NZD(denom)
     }

     # Prediction
     dist.new <- abs(outer(x, x.new, "-"))                         
     u.new <- dist.new / bw    
     
     if(kernel_type == "Epanechnikov"){
       K.val.new <- 3/4 * (1 - u.new^2)       # kernel di Epanechnikov
       K.val.new[ K.val.new < 0 ] <- 0
       K.val.new[ K.val.new > 1 ] <- 0
     }
     
     if(kernel_type =="Gaussian"){
       K.val.new <- dnorm(u.new) 
     }
     
     
     denom.new <- apply( K.val.new, 2, sum )  # need to find a way to calculate this ... 
     hat.y.new <- t(K.val.new) %*% as.matrix(y) / NZD(denom.new)
     idx  <- which(denom.new == 0)
     if(length(idx)>0){
       for(jj in idx){
           C <- -max((-u.new[,jj]^2/2))
           kk <- 1/(sqrt(2*pi))*exp( (-u.new[,jj]^2/2) + C)
           hat.y.new[jj] <- t(kk) %*% as.matrix(y)/sum(kk)
         }
       }
    
    
    
     if(loo){
       return( list( hat.y.loo = hat.y.loo, hat.y = hat.y, hat.y.new = hat.y.new ))
       
     }else{
         return( list( hat.y = hat.y, hat.y.new = hat.y.new ))
     }
     
 }


# Algoritmo con Cross-Validation

ker.reg.cv <- function(x, y, x.new = NULL, kernel_type = "Gaussian")
 {

   # Nadaraya-Watson estimator with GCV choice of bandiwidth
 
     x <- as.vector(x)
     if(is.null(x.new)) x.new <- x
     x.new <- as.vector(x.new)
     y <- as.vector(y)
     n <- length(y)
     dist <- abs(outer(x, x, "-"))                         # matrice delle distanze tra i punti disegno
     dist.1 <- dist[ row(dist) > col(dist) ]
     bw.seq <- 1.1 * quantile(dist.1, seq(.05, 0.5, length=10))     # bandwidth candidates 

   # Generalized Cross-Validation per selezione bandwidht

     CV <- 0
     bw.seq.corr <- bw.seq
     for(j in 1:length(bw.seq))
    {
        bw <- bw.seq[j]
        u <- dist / bw  
        if(kernel_type == "Epanechnikov"){
          K.val <- 3/4 * (1 - u^2)       # kernel di Epanechnikov
          K.val[ K.val < 0 ] <- 0
          K.val[ K.val > 1 ] <- 0
        }
        if(kernel_type == "Gaussian"){
          K.val <- dnorm(u)
        }
        
        diag(K.val) <- 0               # elimino i-ma osservazione per CV
        denom <- apply( K.val, 2, sum )
        while(sum(denom == 0) >= 1)
          {
             bw <- bw * 1.1
             u <- dist / bw  
             if(kernel_type == "Epanechnikov"){
               K.val <- 3/4 * (1 - u^2)       # kernel di Epanechnikov
               K.val[ K.val < 0 ] <- 0
               K.val[ K.val > 1 ] <- 0
             }
             if(kernel_type == "Gaussian"){
               K.val <- dnorm(u)
             }
             diag(K.val) <- 0               # elimino i-ma osservazione per CV
             denom <- apply( K.val, 2, sum )
             #logic <- sum(denom == 0)    
          }
        bw.seq.corr[j] <- bw
        hat.Mat <- K.val / denom                   
        hat.y <- hat.Mat %*% as.matrix(y)  
        CV[j] <- mean( (y - hat.y)^2 )      # Cross-Validation Error
     }
     
     bw.opt <- bw.seq.corr[ order(CV)[1] ]

     dist.new <- abs(outer(x, x.new, "-"))                         
     u.new <- dist.new / bw.opt   
     
     if(kernel_type == "Epanechnikov"){
       K.val.new <- 3/4 * (1 - u.new^2)       # kernel di Epanechnikov
       K.val.new[ K.val.new < 0 ] <- 0
       K.val.new[ K.val.new > 1 ] <- 0
     }
     if(kernel_type == "Gaussian"){
       K.val.new <- dnorm(u.new)
     }
     
     denom.new <- apply( K.val.new, 2, sum )
     while(sum(denom.new == 0) >= 1)
      {
          bw.opt <- bw.opt * 1.1
          u.new <- dist.new / bw.opt  
             
          if(kernel_type == "Epanechnikov"){
              K.val.new <- 3/4 * (1 - u.new^2)       # kernel di Epanechnikov
              K.val.new[ K.val.new < 0 ] <- 0
              K.val.new[ K.val.new > 1 ] <- 0
          }
          if(kernel_type == "Gaussian"){
              K.val.new <- dnorm(u.new)
          }
          diag(K.val.new) <- 0               # elimino i-ma osservazione per CV
          denom.new <- apply( K.val.new, 2, sum )
      }
     u.new <- dist.new / bw.opt  
     
     if(kernel_type == "Epanechnikov"){
       K.val.new <- 3/4 * (1 - u.new^2)       # kernel di Epanechnikov
       K.val.new[ K.val.new < 0 ] <- 0
       K.val.new[ K.val.new > 1 ] <- 0
     }  
     if(kernel_type == "Gaussian"){
       K.val.new <- dnorm(u.new)
     }
     
     denom.new <- apply( K.val.new, 2, sum )
     hat.y.new <- t(K.val.new) %*% as.matrix(y) / denom.new   
     
     # training error 
     u <- dist / bw.opt  
     
     if(kernel_type == "Epanechnikov"){
       K.val <- 3/4 * (1 - u^2)       # kernel di Epanechnikov
       K.val[ K.val < 0 ] <- 0
       K.val[ K.val > 1 ] <- 0
     }
     if(kernel_type == "Gaussian"){
       K.val <- dnorm(u)
     }
     
     denom <- apply( K.val, 2, sum )
     hat.y <- K.val %*% as.matrix(y) / denom
     diag(K.val) <- 0               # elimino i-ma osservazione per CV
     denom <- apply( K.val, 2, sum )
     hat.y.loo <- K.val %*% as.matrix(y) / denom
    
     list(hat.y.loo = hat.y.loo,  hat.y = hat.y, hat.y.new = hat.y.new, bw.opt = bw.opt )
 }



##############################################################
# SINGLE INDEX MODEL
##############################################################

fsim <- function(X.in, Y.in, X.out, Y.out, grid, t_range, nk = 5, d = 4, kernel_type = "Gaussian")
{
  
  p <- dim(X.in)[2]
  dd<- 4
  grid0 <- seq(t_range[1],t_range[2], 1/(10*(p-1))) # in case of not evenly spaced
  knot <- quantile(grid0, (1:nk)/(nk+1) )
  delta <- sort(c(rep(range(grid0), dd), knot)) #exterior knots
  B <- spline.des( delta, grid, dd )$design

  #B <- bs(grid, df = nk + 4,intercept=T)
  
 
  gamma <- t(t(c(rep(1, nk + d))))  
  gamma <- gamma / sqrt( mean( (B %*% gamma)^2 ) )
  #theta <- B %*% gamma 
  Z.in <- X.in %*% B / p 
  X.0 <- Z.in %*% gamma

  soglia <- 1e-2
  err2 <- 100
  d.err <- 100
  h <- 1
  while( (d.err >= soglia) && (h < 100))
  {
        #print(c(h,  err2[h], err2[h] - err2[h-1], err2[h] - err2[h-1]==0))
        h <- h +1
        Y.step <- Y.in
        bw<-  ker.reg.cv(as.vector(X.0), as.vector(Y.step), as.vector(X.0), kernel_type = kernel_type)$bw.opt
        ris <- optim( par=gamma, fn=fun.ob, Y=Y.step, Z=Z.in, B=B, bw=bw )
        gamma <- ris$par / sqrt( mean((B %*% ris$par)^2 ) ) 
        X.0 <- Z.in %*% gamma     
        ker.obj <- ker.reg(X.0, Y.step, X.0, bw, loo = TRUE, kernel_type  = kernel_type )
        Y.hat <- ker.obj$hat.y
        err2[h] <- mean( (Y.hat - Y.in)^2 ) / var(Y.in)
        if(err2[h] - err2[h-1] <= 0) d.err <- abs( err2[h] - err2[h-1] ) else d.err <- soglia + 1      
  }

  #theta <- B %*% gamma
  W.in <- ( X.in %*% B ) %*% gamma / p
  W.out <- ( X.out %*% B ) %*% gamma / p
  hat <- ker.reg.cv( as.vector(W.in), as.vector(Y.in), as.vector(W.out), kernel_type  = kernel_type )
  bw.opt <- hat$bw.opt
  hat.Y.in <- hat$hat.y 
  hat.Y.in.loo <- hat$hat.y.loo # loo error 
  hat.Y.out <- hat$hat.y.new

  RMSE.in <- mean( ( Y.in - hat.Y.in )^2 ) / var(Y.in)
  RMSE.out <- mean( ( Y.out - hat.Y.out )^2 ) / var(Y.out)
 
  return(list(B  = B, gamma = gamma, hat.Y.in = hat.Y.in, hat.Y.in.loo = hat.Y.in.loo, hat.Y.out = hat.Y.out, 
              RMSE.in = RMSE.in, RMSE.out = RMSE.out, 
               W.in = W.in, W.out = W.out, Y.in = Y.in,
              bw.opt = bw.opt))
  
}

predict.fsim <- function(model, X.out, Y.out, kernel_type ){
  
  p <- dim(X.out)[2]
  #theta <- model$B %*% model$gamma
  W.out <- ( X.out %*% model$B ) %*% model$gamma / p
  ker.obj <- ker.reg(as.vector(model$W.in), as.vector(model$Y.in), as.vector(W.out), bw = model$bw.opt, kernel_type =kernel_type )
  hat.Y.out <- ker.obj$hat.y.new
  if(!missing(Y.out))  {
    err_test<- mse(hat.Y.out, Y.out)
    return(list(pred = hat.Y.out, err_test = err_test))
  }else{
    return(pred = hat.Y.out)
  }
}

## how to store the SIM objects .... 
FPPR <- function(x_train, y_train, x_val, y_val, x_test, y_test, grid, t_range, niter, nknots = 5:8, make_prediction = TRUE, precision = 6, kernel_type  = "Gaussian")
{
  
  SIM_objs <- list()
  n_train <- nrow(x_train)
  f_train_t <-  f_val_t <- 0
  SIM_list <- list()
  
  err_train <- err_val <- J <- rep(NA, niter)
  u_train <- y_train
  for(i in 1:niter){
    print(paste(i, "th iteration"))
    AIC <- rep(NA, length(nknots))
    for( j in 1:length(nknots)){
      
      SIM_objs[[j]] <-  fsim(X.in = x_train, Y.in = u_train, X.out = x_val, Y.out = y_val, grid = grid, t_range = t_range,  nk = nknots, d = 4, kernel_type  = kernel_type)
      cv_error <- mean((SIM_objs[[j]]$hat.Y.in.loo - y_train)^2)
      AIC[j] <- n_train*log(cv_error) + 2*(nknots[j])  
    }
    
    J[i] <-which(AIC == min(AIC))
    SIM_list[[i]] <- SIM_objs[[which(AIC == min(AIC))]]
    
    f_train_t <-  f_train_t + as.numeric(SIM_list[[i]]$hat.Y.in)
    f_val_t <-  f_val_t + predict.fsim(SIM_list[[i]], x_val, y_val, kernel_type = kernel_type)$pred
    u_train  <- y_train - f_train_t
    # save error and loss 
    err_train[i] <- mse(f_train_t, y_train)
    err_val[i] <- mse(f_val_t, y_val)
    
    if(i == 1){
      early_stop <- 1
      f_train_early <- f_train_t
      f_val_early <- f_val_t
    }else{
      if(round(err_val[i], precision) < min(round(err_val[1:(i-1)], precision))){
        early_stop <- i
        f_train_early <- f_train_t
        f_val_early <- f_val_t
      }   
    }
  }
  
  model <- list(J = J, kernel_type = kernel_type, early_stop = early_stop, err_val = err_val, err_train = err_train, SIM_list = SIM_list )
  if(make_prediction){
    if(!missing(y_test)){
      tmp_predict <- FPPR.predict(model, newx = x_test, newy = y_test)
      model <- c(model, list(f_val_t =   f_val_early,  f_train_t =   f_train_early, f_test_t = tmp_predict$pred, err_test = tmp_predict$err_test))
    }else{
      tmp_predict <- FPPR.predict(model, newx = x_test)
      model <- c(model, list(f_val_t =   f_val_early,  f_train_t =   f_train_early,f_test_t = tmp_predict$pred))
    }
  }
  
  return(model)
}

FPPR.predict <- function(model, newx, newy){
  
  list2env(setNames(model,paste0(names(model))), envir = environment()) 
  
  f_test_t  <- 0
  err_test <- rep(NA,early_stop)
  
  
  for(i in 1: early_stop){
    
    if(!missing(newy)){
      f_test_t <- f_test_t + predict.fsim(SIM_list[[i]], X.out = newx, Y.out = newy, kernel_type  = kernel_type )$pred
      err_test[i] <- mse(f_test_t, newy)
    } else{
      f_test_t <- f_test_t + predict.fsim(SIM_list[[i]], X.out = newx, Y.out = newy, kernel_type  = kernel_type )$pred
      
    }
  }
  
  if(!missing(newy)){
    return(list(pred = f_test_t, err_test = err_test))
  }else{
    return(pred = f_test_t)
  }
}



