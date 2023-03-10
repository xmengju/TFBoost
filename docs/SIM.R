NZD <- function(a) {
  ifelse(a<0,pmin(-.Machine$double.eps,a),pmax(.Machine$double.eps,a))
}

Ker.gaussian <- function(x){
   dnorm(x)
}


# options for optim.method ("Nelder-Mead", "BFGS", "CG")
SIM.control <- function(optim.method = "Nelder-Mead", optim.maxattempts = 10, optim.reltol=  10*sqrt(.Machine$double.eps), 
                        optim.abstol = 10*.Machine$double.eps, optim.maxit = 1000, save.data = TRUE,   
                       Ker.type  = "gaussian", cv_error = FALSE){
  if(Ker.type == "gaussian"){
    Ker = Ker.gaussian
    kernel = gaussK
  }
  
  return(list(optim.method = optim.method,  optim.maxattempts = optim.maxattempts, optim.reltol = optim.reltol, 
              optim.abstol = optim.abstol, optim.maxit = optim.maxit, save.data = save.data, Ker.type = Ker.type, kernel = kernel, Ker = Ker, cv_error = cv_error))
  
}

constant.loss <- function(param,x,y, kernel, Ker) {
  
  beta <-   c(1, param) 
  loss <- function(x, y, beta) {
    index <- as.numeric(x %*% beta)
    #h <- thumbBw(index, y, deg = 0, kernel)
    #h <- glkerns(x = index, y = y)$bandwidth
    h <- suppressWarnings(regCVBwSelC(x = index, y = y, deg = 0, kernel = kernel))
    #print(c(glkerns(x = index, y = y)$bandwidth,  thumbBw(x = index, y = y, deg = 0, kernel), regCVBwSelC(x = index, y = y, deg = 0, kernel = kernel)))
    index1 <- outer(index, index, "-")
    diag(index1) <- Inf
    index2 <-as.matrix(sweep(index1,1,h,FUN="/"))
    G<-t(Ker(-index2))
    rsum<-rowSums(G,na.rm = TRUE)
    S =G/NZD(rsum)
    y_pred <- S%*%y
    t.ret <- mean((y-y_pred)^2)
    return(t.ret)
  }

  return(loss(x,y,beta))
}

linear.loss <- function(param,x,y, kernel, Ker) {
  
  beta <-   c(1, param) 

  loss <- function(x, y, beta) {
    
    index <- as.numeric(x %*% beta)
    
    #h <- thumbBw(index, y, deg = 1, kernel)
    #kre1 <- np::npreg(y~ index, newdata= data.frame(index), bws = h, regtype = "ll", ckertype="gaussian", ckerorder=2)
    h <- suppressWarnings(regCVBwSelC(x = index, y = y, deg = 1, kernel = kernel))
    index1 <- outer(index, index, "-")
    #print(c(glkerns(x = index, y = y)$bandwidth,  find.h(index1, y, Ker, type = "linear"),  thumbBw(x = index, y = y, deg = 1, kernel), regCVBwSelC(x = index, y = y, deg = 1, kernel = kernel)))
    index2 <- index1/h
    S0 = colMeans(Ker(-index2))
    S1 = colMeans(index1*Ker(-index2))
    S2 = colMeans(index1^2*Ker(-index2))
    S00 <- matrix(1, nrow = nrow(x))%*%S0
    S11 <- matrix(1, nrow = nrow(x))%*%S1
    S22 <- matrix(1, nrow = nrow(x))%*%S2
    WW <- ((S22 - S11 *index1)/(S22*S00 - S11^2)) *Ker(-index2)/nrow(x)
    pred <- t(WW)%*%y
    idx <- which(diag(WW)==1)
    #print(c(h, idx))
    if(length(idx)>0){
      #print("diagnal!")
      LOOV <- mean( ((pred[-idx] - y[-idx])/(1-diag(t(WW))[-idx])) ^2)
      #print(c(length(idx), LOOV))
    }else{
      LOOV <- mean(((pred - y)/(1-diag(t(WW)))) ^2)
    }
    #LOOV <- mean( ((pred - y)/(1-diag(t(WW)))) ^2)
    #print(c(h, LOOV))
    return(LOOV)
  }
  
    return(loss(x,y,beta))

}

# if lambda not equal 0, then it fits single index kernel ridge regression
SIM <- function(x, y, newx, nmulti, random.seed = 42, type = "linear",  control = SIM.control()) {
    
     list2env(setNames(control,paste0(names(control))), envir = environment()) 
    
     if(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        oldseed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        on.exit(assign(".Random.seed", oldseed, envir = .GlobalEnv))
      }
  
      set.seed(random.seed)
      MaxPenalty <- .Machine$double.xmax
      p <- ncol(x)
      n <- nrow(x)

      fval.min <- .Machine$double.xmax
      fval.value <- numeric(nmulti)
      
      if(type == "constant"){
        optim.fn <-  function(param){constant.loss(param,x,y, kernel, Ker)}
      }
      
      if(type == "linear"){
        optim.fn <-  function(param){linear.loss(param,x,y,kernel,Ker)}
      }
      
      converg <- 0
      for(i in 1:nmulti) {
        
        optim.control <- list(abstol=optim.abstol, reltol=optim.reltol, maxit=optim.maxit)
        if(i == 1) {
          ols.fit <- lm(y~x-1)
          beta_tilde <- coef(ols.fit)
          beta_tilde[which(is.na(beta_tilde))] <- 0
          fit <- fitted(ols.fit)
        } else {
          beta_tilde <- coef(ols.fit)
          beta_tilde[which(is.na(beta_tilde))] <- 0
          beta_tilde <- runif(p, min=0.5, max=1.5)*beta_tilde
        }
          
        if(sum(abs(beta_tilde)) == 0){
          beta_tilde <- runif(length(beta_tilde))
        }

        
        optim.parm <- beta_tilde[2:p]/sign(beta_tilde[1])
        optim.return <- optim(optim.parm,fn=optim.fn,gr=NULL,method=optim.method,control=optim.control)
        param <-  optim.return$par
        attempts <- 0
        
        if(optim.return$convergence == 0){
          converg <- 1
        }
        # if does not converge
        while((optim.return$convergence != 0) && (attempts <= optim.maxattempts)) {
          
          print(c("attempts",attempts))
          attempts <- attempts + 1
          
          if(optim.return$convergence == 1){
            if(optim.control$maxit < (2^32/10))
                optim.control$maxit <- 10*optim.control$maxit
           # else
            #      stop(paste("optim failed to converge after optim.maxattempts = ", optim.maxattempts, " iterations."))
          }
              
          if(optim.return$convergence == 10){
            beta_tilde <- coef(ols.fit)
            beta_tilde[which(is.na(beta_tilde))] <- 0
            beta_tilde <- runif(p, min=0.5, max=1.5)*beta_tilde
            optim.parm <- beta_tilde[2:p]/sign(beta_tilde[1])
          }
          
          #print(c(beta_tilde[1], optim.parm))
          optim.return <- optim(optim.parm,fn=optim.fn,gr=NULL,method=optim.method,control=optim.control)
          
          if(optim.return$convergence == 0){
            converg <- 1
            break
          }
        }
            
        fval.value[i] <- optim.return$value
        if(optim.return$value < fval.min) {
            param <-  optim.return$par
            fval.min <- optim.return$value
            best <- i
          }
      }
      
      if(converg == 0){
        warning("SIM not converging")
      }
      beta = c(1,param)
      index <- as.numeric(x %*% beta)
      
      if(type == "constant"){
        #h <- thumbBw(index, y, deg = 0, kernel)
        h <- suppressWarnings(regCVBwSelC(x = index, y = y, deg = 0, kernel = kernel))
        #h <- glkerns(x = index, y = y)$bandwidth
      }
      if(type == "linear"){
        h <- suppressWarnings(regCVBwSelC(x = index, y = y, deg = 1, kernel = kernel))
        #h <- thumbBw(index, y, deg = 1, kernel)
      }
      
      res <- list(beta = beta, h = h, fval = fval.min, x = x, y = y, Ker = Ker, type = type)
     
      sim_train_pred <- SIM.predict(SIM_obj = res, newx = x, cv_error = cv_error)
      res$pred_train <-    sim_train_pred$pred
      
      if(cv_error == TRUE){
        res$LOO <-  sim_train_pred$LOO
      }
        
      if(!missing(newx)){
          res$pred_test <- SIM.predict(SIM_obj = res, newx = newx, cv_error = FALSE)$pred
      }
          
      if(!save.data){
          res$x <- NULL
          res$y <- NULL
      }
    return(res)
}

SIM.predict <- function(SIM_obj, newx, cv_error = FALSE){
  
    list2env(setNames(SIM_obj,paste0(names(SIM_obj))), envir = environment()) 
    index <- as.numeric(as.matrix(x) %*% beta)
    new_index <- as.numeric(as.matrix(newx) %*% beta)

    index1 <- outer(index, new_index, "-")
    index2 <-index1/h
    
    if(type == "constant"){
      G<- t(Ker(-index2))
      rsum<-rowSums(G,na.rm = TRUE)
      S =G/NZD(rsum)
      pred <- S%*%y
      
      if(cv_error == TRUE){
        G<-t(Ker(-index2))
        diag(G) <- 0 
        rsum<-rowSums(G,na.rm = TRUE)
        S =G/NZD(rsum)
        pred_tmp <- S%*%y
        LOO <- mean((pred_tmp - y) ^2)
      }
    }
    
    if(type == "linear"){
      S0 = colMeans(Ker(-index2))
      S1 = colMeans(index1*Ker(-index2))
      S2 = colMeans(index1^2*Ker(-index2))
      S00 <- matrix(1, nrow = nrow(x))%*%S0
      S11 <- matrix(1, nrow = nrow(x))%*%S1
      S22 <- matrix(1, nrow = nrow(x))%*%S2
      WW <- ((S22 - S11 *index1)/(S22*S00 - S11^2)) *Ker(-index2)/nrow(x)
      pred <- t(WW)%*%y
      if(cv_error == TRUE){
        idx<- which(diag(WW)==1)
        if(length(idx)>0){
          LOO <- mean( ((pred[-idx] - y[-idx])/(1-diag(t(WW))[-idx])) ^2)
        }else{
         LOO <- mean( ((pred - y)/(1-diag(WW))) ^2)
        }
      }
    }
    
    if(cv_error == TRUE){
      return(list(pred = pred, LOO = LOO))
    }else{
      return(list(pred = pred))
      
    }
}



find.h <-function(index1, y, Ker, type = "constant"){

  if(type == "constant"){
    ff <- function(h){
      index2 <-as.matrix(sweep(index1,1,h,FUN="/"))
      G<-t(Ker(-index2))
      diag(G) <- 0 
      rsum<-rowSums(G,na.rm = TRUE)
      S <- G/NZD(rsum)
      y_pred <- S%*%y
      t.ret <- mean((y-y_pred)^2)
      return(t.ret)
    }
  }else{
    ff <- function(h){
      index2 <-as.matrix(sweep(index1,1,h,FUN="/"))
      S0 = colMeans(Ker(-index2))
      S1 = colMeans(index1*Ker(-index2))
      S2 = colMeans(index1^2*Ker(-index2))
      S00 <- matrix(1, nrow = nrow(index2))%*%S0
      S11 <- matrix(1, nrow = nrow(index2))%*%S1
      S22 <- matrix(1, nrow = nrow(index2))%*%S2
      WW <- ((S22 - S11 *index1)/(S22*S00 - S11^2)) *Ker(-index2)/nrow(index2)
      pred <- t(WW)%*%y
      LOO <- mean( ((pred - y)/(1-diag(WW))) ^2)
      return(LOO)
    }
  }
  
  plot(seq(0.01,5,0.01), unlist(lapply(seq(0.01,5,0.01),FUN =  ff)))
  optimize_res <- optimize(f = ff, interval = c(0,3))
  return(optimize_res$minimum)
}





# 
# sim.krr.gasser.loss <- function(param,x,y, lambda, Ker){
#   
#   beta <-  c(1, param) #beta1 = 1
#   
#   loss <- function(x,y, beta, h) {
#     
#     index <- as.numeric(x %*% beta)
#     index1 <- abs(outer(index, index, "-"))
#     h <- glkerns(x = index, y = y)$bandwidth
#     index2 <-as.matrix(sweep(index1,1,h,FUN="/"))
#     G<-Ker(index2)
#     tmp_alpha_r <- try(qr.solve(G  + lambda *diag(nrow(G))), silent = T)
#     alpha_r <- tmp_alpha_r%*%y
#     y_pred <- G%*%alpha_r
#     t.ret <- sum((y-y_pred)^2) + lambda * t(alpha_r)%*%G%*%alpha_r
#     return(t.ret)
#   }
#   
#   return(loss(x,y,beta, h))
#   
# }
# 
# sim.krr.cv.loss <- function(param,x,y, lambda, Ker){
#   
#   beta <-  c(1, param[1:(length(param)-1)]) #beta1 = 1
#   h <- param[length(param)]
#   
#   loss <- function(x,y, beta, h) {
#     
#     index <- as.numeric(x %*% beta)
#     index1 <- abs(outer(index, index, "-"))
#     index2 <-as.matrix(sweep(index1,1,h,FUN="/"))
#     G<-Ker(index2)
#     tmp_alpha_r <- try(qr.solve(G  + lambda *diag(nrow(G))), silent = T)
#     alpha_r <- tmp_alpha_r%*%y
#     y_pred <- G%*%alpha_r
#     t.ret <- sum((y-y_pred)^2) + lambda * t(alpha_r)%*%G%*%alpha_r
#     return(t.ret)
#   }
#   
#   if(h > 0) {
#     return(loss(x,y,beta, h))
#   }else {
#     return(MaxPenalty)
#   }
# }
# 
# sim.nw.fixed.loss <- function(param,x,y, h, Ker) {
#   
#   beta <- param
#   loss <- function(x,y,beta) {
#     
#     index <- as.numeric(x %*% beta)
#     index1 <- abs(outer(index, index, "-"))
#     diag(index1) <- Inf
#     index2 <-as.matrix(sweep(index1,1,h,FUN="/"))
#     G<-Ker(index2)
#     rsum<-rowSums(G,na.rm = TRUE)
#     S =G/NZD(rsum)
#     y_pred <- S%*%y
#     t.ret <- mean((y-y_pred)^2)
#     return(t.ret)
#   }
#   
#   return(loss(x,y,beta))
# }
# 
# sim.krr.fixed.loss <- function(param,x,y,h, lambda, Ker){
#   
#   beta <- param
#   loss <- function(x,y,beta) {
#     index <- as.numeric(x %*% beta)
#     index1 <- abs(outer(index, index, "-"))
#     index2 <-as.matrix(sweep(index1,1,h,FUN="/"))
#     G<-Ker(index2)
#     tmp_alpha_r <- qr.solve(G  + lambda *diag(nrow(G)))
#     alpha_r <- tmp_alpha_r%*%y
#     y_pred <- G%*%alpha_r
#     t.ret <- sum((y-y_pred)^2) + lambda * t(alpha_r)%*%G%*%alpha_r
#     return(t.ret)
#   }
#   
#   return(loss(x,y,beta))
#   
# }
# sim.nw.gasser.loss <- function(param,x,y, Ker) {
#   
#   beta <-   c(1, param) # beta1 = 1
#   #beta <-   param 
#   loss <- function(x, y, beta) {
#     
#     index <- as.numeric(x %*% beta)
#     index1 <- abs(outer(index, index, "-"))
#     diag(index1) <- Inf
#     #h <- glkerns(x = index, y = y)$bandwidth
#     h <- find.h(index1, y)
#     index2 <-as.matrix(sweep(index1,1,h,FUN="/"))
#     G<-Ker(index2)
#     rsum<-rowSums(G,na.rm = TRUE)
#     S =G/NZD(rsum)
#     y_pred <- S%*%y
#     #t.ret <- mean((y-y_pred)^2) + MaxPenalty  * abs(sqrt(sum(beta^2)) -1)
#     t.ret <- mean((y-y_pred)^2)
#     
#     return(t.ret)
#   }
#   
#     return(loss(x,y,beta))
# 
# }


#Silverman <- function(val, n){
#  return( max(1.059224*apply(as.matrix(val),2,sd)*n^(-1/5), .Machine$double.eps))
#}

