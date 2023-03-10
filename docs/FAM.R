FAM_update <- function (Y, Lx, Lt, nEval = 51, newLx = NULL, newLt = NULL, 
                 bwMethod = 0, alpha = 0.7, supp = c(-2, 2), optns = NULL) 
{
  if (is.null(optns) == TRUE) {
    optns <- list()
  }
  n <- length(Y)
  
  if(!missing(optns)){
    tmpFPCA <- FPCA(Ly = Lx, Lt = Lt, optns = optns)
  }else{
    tmpFPCA <- FPCA(Ly = Lx, Lt = Lt)
  }
  XiStd <- t(t(tmpFPCA$xiEst)/sqrt(tmpFPCA$lambda))
  d <- length(tmpFPCA$lambda)
  estLambda <- tmpFPCA$lambda
  estEigen <- tmpFPCA$phi
  workGrid <- tmpFPCA$workGrid
  N <- xiStdGrid <- c()
  
  predobj <- predict(object = tmpFPCA, newLy = newLx, 
                     newLt = newLt, K = d,  xiMethod = "IN")
  xiStdGrid <- predobj$scores %*% diag(1/sqrt(tmpFPCA$lambda))
  N <- nrow(xiStdGrid)
  
  h <- c()
  for (j in 1:d) {
    if (bwMethod > 0) {
      h[j] <- suppressWarnings(CVLwls1D(y = (Y - mean(Y)), 
                                        t = XiStd[, j], kernel = "gauss", npoly = 1, nder = 0, 
                                        dataType = "DenseWithMV", kFolds = bwMethod))
    }
    else {
      h[j] <- suppressWarnings(GCVLwls1D1(yy = (Y - mean(Y)), 
                                          tt = XiStd[, j], kernel = "gauss", npoly = 1, 
                              nder = 0, dataType = "DenseWithMV")$bOpt)
    }
  }
  h <- alpha * h
  fam <- matrix(nrow = N, ncol = d)
  for (j in 1:d) {
    xiTmp <- sort(xiStdGrid[, j])
    fitTmp <- Lwls1D(bw = h[j], kernel_type = "gauss", xin = sort(XiStd[, 
                                                                        j]), yin = (Y[order(XiStd[, j])] - mean(Y)), xout = xiTmp, 
                     npoly = 1, nder = 0)
    fam[, j] <- fitTmp[match(xiStdGrid[, j], xiTmp)]
    fam[, j] <- fam[, j] - mean(fam[, j])
  }
  yMean <- mean(Y)
  xiGrid <- xiStdGrid %*% diag(sqrt(estLambda))
  bw <- h * sqrt(estLambda)
  phi <- estEigen
  fit <- list(mu = yMean, fam = fam, xi = xiGrid, bw = bw, 
              lambda = estLambda, phi = phi, workGrid = workGrid)
  return(fit)
}




FAM_FD <- function(x_train_fd, y_train, x_test_fd, y_test, niter){
  
  N <- length(y_train)
  
  pca_est <- pca.fd(x_train_fd, nharm  = niter)
  d <- max(which(cumsum(pca_est$varprop)<=0.95))
  #train_scores <-  inprod(x_train_fd,pca_est$harmonics[1:d])
  train_scores <-  pca_est$scores[,1:d]
  #train_scores <-   train_scores - matrix(1, nrow = N) %*%inprod(pca_est$meanfd, pca_est$harmonics[1:d])
  XiStd  = t(t(train_scores)/sqrt(pca_est$values[1:d]))

  h <- c()
  for (j in 1:d) {
      h[j] <- suppressWarnings(GCVLwls1D1(y=(y_train-mean(y_train)),tt=XiStd[,j],kernel='gauss',npoly=1,nder=0,dataType='Dense')$bOpt)
  } 
  
  fam <- matrix(nrow=N,ncol=d)
  xiStdGrid <- XiStd
  for (j in 1:d) {
    xiTmp <- sort(xiStdGrid[,j])
    fitTmp <- Lwls1D(bw=h[j],kernel_type='gauss',xin=sort(XiStd[,j]),yin=(y_train[order(XiStd[,j])]-mean(y_train)),xout=xiTmp,npoly=1,nder=0)
    fam[,j] <- fitTmp[match(xiStdGrid[,j],xiTmp)]
    fam[,j] <- fam[,j] - mean(fam[,j])
  }
  yMean <- mean(y_train)
  
  
  ## make predictions of test data 
  test_scores <-  inprod(x_test_fd,pca_est$harmonics[1:d]) -matrix(1, nrow = ncol(x_test_fd$coefs)) %*%inprod(pca_est$meanfd, pca_est$harmonics[1:d])
  test_XiStd  = t(t(test_scores)/sqrt(pca_est$values[1:d]))
  
  fam_test <- matrix(nrow=length(y_test),ncol=d)
  xiStdGrid <- test_XiStd 
  for (j in 1:d) {
    xiTmp <- sort(xiStdGrid[,j])
    fitTmp <- Lwls1D(bw=h[j],kernel_type='gauss',xin=sort(XiStd[,j]),yin=(y_train[order(XiStd[,j])]-mean(y_train)),xout=xiTmp,npoly=1,nder=0)
    fam_test[,j] <- fitTmp[match(xiStdGrid[,j],xiTmp)]
    fam_test[,j] <- fam_test[,j] - mean(fam_test[,j])
  }
  
  f_trains <- matrix(NA, N, d)
  f_tests <- matrix(NA, length(y_test), d)
  
  for(i in 1:d){
    if(i == 1){
      f_trains[,i] <- yMean + fam[,1]
      f_tests[,i] <- yMean + fam_test[,1]
    }else{
      f_trains[,i] <- yMean + apply(fam[,1:i], 1, sum) 
      f_tests[,i] <- yMean + apply(fam_test[,1:i], 1, sum) 
    }
  }
 
  err_train <-  apply(f_trains - y_train%*% matrix(rep(1,d), ncol = d), 2, FUN = function(x) mean(x^2))
  err_test <-  apply(f_tests - y_test%*% matrix(rep(1,d), ncol = d), 2, FUN = function(x) mean(x^2))
  
  
    
  model <- list(early_stop = d, err_train = err_train, err_test = err_test, f_test_t = f_tests[,d], f_train_t = f_trains[,d] )
  return(model)
  
}

FAM_dense <- function(Y,Lx,Lt,nEval=51,newLx=NULL,newLt=NULL,bwMethod=0,alpha=0.7,supp=c(-2,2),optns=NULL){
  
  
  if (is.null(optns)==TRUE) {
    optns <- list() 
  }
  
  n <- length(Y)
  
  tmpFPCA <- FPCA(Ly=Lx, Lt=Lt, optns=optns)
  
  XiStd <- t(t(tmpFPCA$xiEst)/sqrt(tmpFPCA$lambda))
  d <- length(tmpFPCA$lambda)
  
  estLambda <- tmpFPCA$lambda
  estEigen <- tmpFPCA$phi
  workGrid <- tmpFPCA$workGrid
  
  N <- xiStdGrid <- c()
  if (is.null(newLx)==TRUE | is.null(newLt)==TRUE) {
    
    if (nEval==0) {
      N <- nrow(XiStd)
      xiStdGrid <- XiStd
    } else {
      N <- nEval
      xiStdGrid <- matrix(rep(seq(supp[1],supp[2],length.out=N),d),nrow=N,ncol=d)
    }
    
  } else {
    predobj <-  predict.FPCA(object = tmpFPCA,newLy=newLx,newLt=newLt,K=d)
    xiStdGrid <- predobj$scores%*%diag(1/sqrt(tmpFPCA$lambda))
    N <- nrow(xiStdGrid)
  }
  
  h <- c()
  for (j in 1:d) {
    if (bwMethod>0) {
      h[j] <- suppressWarnings(CVLwls1D(y=(Y-mean(Y)),t=XiStd[,j],kernel='epan',npoly=1,nder=0,dataType='Dense',kFolds=bwMethod))
    } else {
      h[j] <- suppressWarnings(GCVLwls1D1(yy=(Y-mean(Y)),tt=XiStd[,j],kernel='epan',npoly=1,nder=0,dataType='Dense')$bOpt)
    }
  } 
  
  h <- alpha*h
  
  fam <- matrix(nrow=N,ncol=d)
  for (j in 1:d) {
    xiTmp <- sort(xiStdGrid[,j])
    fitTmp <- Lwls1D(bw=h[j],kernel_type='epan',xin=sort(XiStd[,j]),yin=(Y[order(XiStd[,j])]-mean(Y)),xout=xiTmp,npoly=1,nder=0)
    fam[,j] <- fitTmp[match(xiStdGrid[,j],xiTmp)]
    
    fam[,j] <- fam[,j] - mean(fam[,j])
  }
  yMean <- mean(Y)
  
  xiGrid <- xiStdGrid%*%diag(sqrt(estLambda))
  bw <- h*sqrt(estLambda)
  phi <- estEigen
  
  fit <- list(mu=yMean, fam=fam, xi=xiGrid, bw=bw, lambda=estLambda, phi=phi, workGrid=workGrid)
  
  return(fit)
  
}
