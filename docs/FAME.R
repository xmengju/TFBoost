"fameafit" <- function(family,gamfit,a,x,S,alpha,sigmax,W,data,p,r,Theta,pc){
  
    if (family$family=="gaussian"){
      w <- 1/mean(gamfit$res^2)
    }else{
      w <- 1
    }
  return (optim(par=a,fn=famecalclike,method="BFGS",x=data$x-S%*%alpha,sigmax=sigmax,
  W=W,S=S,data=data,gamfit=gamfit,p=p,family=family,w=w,r=r,Theta=Theta,pc=pc,control =
  list(maxit=1)))
}

"famecalclike" <- function(a,x,gamfit,sigmax,data,S,W,p,family,w,r,Theta,pc){
  
    a <- matrix(a,length(a)/p,p)
    if (pc){
      Theta <- matrix(t(svd(var(t(W%*%t(a))))$u[,1:r])%*%W,p,r)
    }
    if("z"%in% names(data)){
      newdata <- data.frame(a=a%*%Theta, z = data$z)
    }else{
      newdata <- data.frame(a=a%*%Theta)
    }
    sum(family$dev(data$y,predict.gam(gamfit,newdata,type="response"),w))+ sum((x -((S%*%W)*a[data$curve,])%*%rep(1,p))^2)/sigmax
}


"famecalclike2" <- function(theta,a,data,p,r,family,model){
  
    rot <- matrix(theta,p,r)
    
    if("z"%in%names(data)){
      data1 <- data.frame(y=data$y,a=a%*%rot, z = data$z)
    }else{
      data1 <- data.frame(y=data$y,a=a%*%rot)
    }
  
    sum(gam(model,data=data1,family=family)$dev)
}

"fame.predict" <- function(obj, newdata, type = "response"){
  
  N <- length(table(newdata$curve))
  a <- matrix(NA, N, obj$p)  
  W <- obj$W
  S <- obj$BaseS[newdata$time,]
  alpha <- obj$alpha
  
  for(i in 1:N){
      Si <- S[newdata$curve==i,]
      X <- solve(t(Si%*%W)%*%Si%*%W)%*%t(Si%*%W)
      a[i,] <- X%*%(newdata$x[newdata$curve==i]-Si%*%alpha)
  }
  
  if("z" %in% names(newdata)){
    newdat <- data.frame(a=a%*%obj$Theta, z = newdata$z)
  }else{
    newdat <- data.frame(a=a%*%obj$Theta)
  }

  tmp <- predict.gam(obj$gamfit,newdat,type=type)
  if(type == "terms"){
    pred <- attr(tmp, "constant") + t(apply(t(tmp), 2, cumsum)) 
  }else{
    pred <- tmp
  }
  return(pred)
}

"famegamfit" <- function(Theta,a,data,p,r,family,model,maxit=5,W,pc){
  
  like1 <- -famecalclike2(Theta,a=a,data=data,p=p,r=r,family=family,model=model)/2
  if(pc){
    tempTheta <- matrix(t(svd(var(t(W%*%t(a))))$u[,1:r])%*%W,p,r)
    like2<- -famecalclike2(tempTheta,a=a,data=data,p=p,r=r,family=family,model=model)/2
    if(like1<like2){
      Theta <- tempTheta
    }
  }else{
    Theta <- matrix(optim(par=Theta,fn=famecalclike2,a=a,data=data,p=p,r=r,family=family,model=model,control=list(maxit=maxit))$par,p,r)
  }
  
  if("z"%in% names(data)){
    data1 <- data.frame(y=data$y,a=a%*%Theta, z = data$z)
  }else{
    data1 <- data.frame(y=data$y,a=a%*%Theta)
  }
  gamfit <- gam(model,data=data1,family=family)
  list(Theta=Theta,gamfit=gamfit)
}

"fameinit2" <- function(run=10,data,S,p,a,N,q,r,family,model,pc){
  
  current.best <- -10^10
  
  for (j in 1:run){
    
    a <- matrix(rnorm(N*p),N,p)  
    
    for(i in 1:10){
      
      xfit <- famexfit(q,S,N,data,a,p)
      alpha <- xfit$alpha
      W <- xfit$W
      sigmax <- xfit$sigmax
      
      for(i in 1:N){  ## every curve 
        Si <- S[data$curve==i,]
        X <- solve(t(Si%*%W)%*%Si%*%W)%*%t(Si%*%W)  ## projected curve 
        a[i,] <- X%*%(data$x[data$curve==i]-Si%*%alpha) ## a is the representation in lower dimension space 
      }
    }
    
    for(i in 1:5){
      #print(paste("initialization", i, "out of", 5))
      
      Theta <- matrix(rnorm(p*r),p,r)
      yfit <- famegamfit(Theta,a,data,p,r,family,model,maxit=10,W,pc)
      Theta <- yfit$Theta
      gamfit <- yfit$gamfit  
      like <- -(famecalclike(a=a,x=data$x-S%*%alpha,sigmax=sigmax,
                             W=W,S=S,data=data,gamfit=gamfit,p=p,family=family,w=1,r,Theta,pc)+log(sigmax)*length(data$x))/2
      if(like > current.best){
        newa <- a
        current.best <- like
        newTheta <- Theta
        newW <- W
      }
    }
  }
  list(a=newa,Theta=newTheta,W=newW)
}


"famexfit" <- function(q,S,N,data,a,p){
  
  suppressWarnings(largeS <-
    matrix(rep(c(rep(1,q),rep(0,(N)*q)),N),N,N*q,byrow=T)[data$curve,]*matrix(rep(S,N),length(data$curve),N*q)) ## not matrix multiplication! 
  tempA <- outer(diag(q),a)
  tempA <- matrix(aperm(tempA,c(1,3,2,4)),N*q,q*p)
  A <- largeS%*%tempA
  rm(largeS)
  rm(tempA)
  lmfit <- lm(data$x~S+A-1)
  alpha <- lmfit$coef[1:q] # coefficients for the basis 
  W <- matrix(lmfit$coef[-(1:q)],q,p) # seems to be weight for random combinations of columns 
  sigmax <- mean(lmfit$res^2)
  list(W=W,sigmax=sigmax,alpha=alpha)
} 


# data: a list with four or five components described as follows 
# - x (vector): contains each of the evaluated points of the predictor curves i.e. (x_1, x_2, …, x_N) 
# where x_i is a vector containing the observations from curve i (note each curve can have observations at different time points).
# - time (vector) It is the same length as x. Each component corresponds to the time point that the corresponding component of x is measured at.
# This is with respect to grid, which is a component you also supply to famefn (by default grid is 1:100)
# - curve (vector) This is also the same length as x. It will be of the form (1,1,1…,1,2,2,…2…) etc. assuming the first points in x come from curve 1, 
# then the next points from curve 2 etc.
# - y (vector) It has one element for each of the N predictors.
# - z (vector) This argument is optional and represents a length N scalar predictor  

# other arguments
# p: a lower dimension that the basis coefficients can be projected into
# q: cubic B-spline basis dimension (>4)
# r: the number of non-linear projections 
# t_range: range to define the B-spline basis
# maxit: maximum number of iterations
# pert: a small perturbation to get the initialization to work if the time points that the x’s are observed at is not dense enough to get the initial least squares estimates.
# By default pert=0

"famefn" <- function(data,tol=.001,q,p,r,grid=1:100, t_range = c(0,1), maxit=3,family=gaussian(), pert=0,prevfit=NULL,run=4,pc=F){
  
  like.old <- 1  
  like.new <- 2
  N <- length(table(data$curve))
  
  if("z" %in% names(data)){
    model <- as.formula(paste("y",paste(paste( "s(a.", (1:r),", bs = 'cs')",sep="",collapse="+"), "+s(z, bs = 'cs')"),sep="~"))
  }else{
    model <- as.formula(paste("y",paste( "s(a.", (1:r),", bs = 'cs')",sep="",collapse="+"),sep="~"))
  }
  
  if(r==1){
    model <- y~s(a,bs = "cs")
  }
  
  nknot <- q - 4  # number of interior knots of B spline
  dd <- 4; # cubic 
  grid0 <- seq(t_range[1],t_range[2], 1/(10*(length(grid)-1))) 
  knot <- quantile(grid0, (1:nknot/(nknot+1)))
  delta <- sort(c(rep(range(grid0), dd), knot)) #add exterior knots
  B<- spline.des(delta, grid, dd)$design
  BaseS <- svd(B)$u
  S <- BaseS[data$time,] #Basis 
  
  if(is.null(prevfit)){
    
    prevfit <- fameinit2(run,data,S,p,a,N,q,r,family,model, pc)
    print("init done!")
  }
  
  a <- prevfit$a
  Theta <- prevfit$Theta
  W <- prevfit$W
  iter <- 1
  
  while ((abs(like.new-like.old)/abs(like.old)>tol) & (iter <= maxit)){
    
    yfit <- famegamfit(Theta,a,data,p,r,family,model,20,W,pc)
    Theta <- yfit$Theta
    gamfit <- yfit$gamfit
    xfit <- famexfit(q,S,N,data,a,p)
    W <- xfit$W
    alpha <- xfit$alpha
    sigmax <- xfit$sigmax
    afit <- fameafit(family,gamfit,a,x,S,alpha,sigmax,W,data,p,r,Theta,pc)
    a <- afit$par
    print(paste("Iteration ",iter))
    like.old <- like.new
    like.new <- -(afit$value+log(sigmax)*length(data$x))/2
    print(paste("Log likelihood = ", like.new, "sigmax = ", sigmax))
    iter <- iter+1
  }
  list(a=a,W=W,gamfit=gamfit,sigmax=sigmax,alpha=alpha,Theta=Theta,family=family,p=p,r=r,BaseS=BaseS,q=q,model=model)
}



