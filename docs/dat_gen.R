dat.generate.control <- function(n_train = 200,  n_val = 200, n_test = 200, g_func_no = 1, SNR = 10,  x_type = "ferraty", t_dim = 100){
  
  sd_g_list <- readRDS("mad_g_list.rds")
  pca_f <- readRDS("eigen_f.rds")
  eigen_f <- pca_f[[1]]; mu_f <- pca_f[[2]]
  
  mattern_tmp <- readRDS("mattern.rds")
  mu_m <-  mattern_tmp$mu
  phis <- mattern_tmp$phis
  eigen_m <- list(phis[[1]], phis[[2]])
  
  if(x_type == "ferraty"){
    eigenf <- eigen_f
    muf <- mu_f
    rr <- c(-1,1)
  }else{
    eigenf <- eigen_m
    muf <- mu_m
    rr <- c(0,1)
  }
  
  switch(g_func_no,
         "1" = {
           g.func <- function(x,t){
             (x -muf(t)) *(eigenf[[1]](t)+ eigenf[[2]](t)) 
           }
           int_func <- function(x, t_eval){
             tmp <- riemman(g.func(x, t_eval), t_eval, rr = rr)
             sign(tmp)*(abs(tmp))^{1/3}
           }
         },
         "2" = {
           g.func <- function(x,t){
             x*log(abs(x))/2
           }
           int_func <- function(x, t_eval){ 5*exp(-abs(riemman(g.func(x, t_eval), t_eval, rr = rr)))}
         },
         "3" = {
           g.func <- function(x,t){
             x^2*sin(2*pi*t)
           }
           int_func <- function(x, t_eval){
             tmp <- riemman(g.func(x, t_eval), t_eval, rr = rr)
             return(5/(1+exp(-2*tmp)))
          }
         },
         "4" = {
           g.func.1 <-  function(x,t){
             cos(2*pi*t^2)*x
           }
           g.func.2 <-  function(x,t){
             sin(x)
           }
           
           if(x_type == "ferraty"){
             threshold <- 0
             rr_1 <- c(-1,0)
             rr_2 <- c(0,1)
           }else{
             threshold <- 0.5             
             rr_1 <- c(0,0.5)
             rr_2 <- c(0.5,1)
           }
           int_func <- function(x, t_eval){5*(sqrt(abs(riemman(g.func.1(x[t_eval<threshold], t_eval[t_eval<threshold]), t_eval[t_eval<threshold], rr = rr_1))) +  
                                                + sqrt(abs(riemman(g.func.2(x[t_eval>=threshold], t_eval[t_eval>=threshold]), t_eval[t_eval>=threshold], rr = rr_2))))}
          },
         "5" = {
           g.func <- function(x,t){
             (sin((3*pi)*t/2) + sin(0.5*pi*t))*x
           }
           int_func <- function(x, t_eval){riemman(g.func(x, t_eval), t_eval, rr =rr)}
         }
         
  )
  return(list(n_train = n_train, n_val = n_val, n_test = n_test, g_func_no = g_func_no,  int_func = int_func,  SNR = SNR, t_dim = t_dim, x_type = x_type,   sd_g_list =   sd_g_list))
}


### generate the data  control$sd_g_list is used to control the SNR .... need to double check this .... 

dat.generate <- function(seed, control = dat.generate.control()){
  
  n_train <- control$n_train;  n_val <- control$n_val;  n_test <- control$n_test; 
  t_dim <- control$t_dim  # time dimension 
  n <- n_train + n_val + n_test
  
  g_func_no <- control$g_func_no
  int_func <- control$int_func
  SNR <- control$SNR
  
  sd_g <- control$sd_g_list[[control$x_type]][[control$g_func_no]]
  
  set.seed(seed)
  
  if(control$x_type == "ferraty"){
    t_eval <- seq(-1,1, length.out = 1000)
    a <- runif(n)
    b <- runif(n)
    c <- runif(n, min = -1, max = 1)
    d <- runif(n, min = -2*pi, max = 2*pi)
    x.fun <- function(i){
      return(function(t){a[i] + b[i]*t^2 + c[i]*exp(t) + sin(d[i]*t)})
    }
    tmp <- lapply(1:n, function(i){x.fun(i)(t_eval)})
    x_func <- matrix(unlist(tmp), ncol = 1000, byrow = T)
    
  }else{
    t_eval <- seq(0,1, length.out = 1000)
    #xi_matrix <- matrix(runif(n*4, min = -sqrt(3),max = sqrt(3)), ncol = 4)
    xi_matrix <- matrix(rnorm(n*4), ncol = 4)
    x_func <- generate_mattern(xi_matrix)(t_eval)
  }
  
  g <- unlist(lapply(1:n,  FUN = function(i) {int_func(x_func[i,], t_eval)}))
  S <- sd_g/sqrt(SNR)
  y <- g + S*rnorm(n)
  print(paste("constant is ", round(S, 3)))
  
  train_idx <-1:n_train
  val_idx <-  (n_train+1): (n_train + n_val)
  test_idx <- (n_train + n_val+1):n
  
  tt <- seq(1, 1000, by = 10)
  x_train <- x_func[train_idx,tt]
  x_val <- x_func[val_idx,tt]
  x_test <- x_func[test_idx,tt]
  
  y_train <- y[train_idx]
  y_val <- y[val_idx]
  y_test <- y[test_idx]
  
  dat <- list()
  dat$x <- list(x_train = x_train, x_val = x_val, x_test = x_test)
  dat$y <- list(y_train = y_train, y_val = y_val, y_test = y_test)
  dat$t_eval <- t_eval
  dat$tt <- t_eval[tt]
  dat$g <- g
  dat$S <- S
  return(dat)
}



## compute eigen values and sd of g 
## for mattern we use the true values
cal_g_sd_eigen <- function(g_func_no = 1, x_type = "ferraty", cal_eigen = TRUE){
  
  set.seed(123)
  
  if(!cal_eigen){
    control <- dat.generate.control(g_func_no = g_func_no, x_type = x_type)
  }
  n <- 3000
  
  if(x_type == "ferraty"){
    t_eval <- seq(-1,1, length.out = 1000)
    a <- runif(n)
    b <- runif(n)
    c <- runif(n, min = -1, max = 1)
    d <- runif(n, min = -2*pi, max = 2*pi)
    x.fun <- function(i){
      return(function(t){a[i] + b[i]*t^2 + c[i]*exp(t) + sin(d[i]*t)})
    }
    tmp <- lapply(1:n, function(i){x.fun(i)(t_eval)})
    x_func <- matrix(unlist(tmp), ncol = 1000, byrow = T)
    rr <- c(-1,1)
  }else{
    t_eval <- seq(0,1, length.out = 1000)
    xi_matrix <- matrix(rnorm(n*4), ncol = 4)
    x_func <- generate_mattern(xi_matrix)(t_eval)
    rr <- c(0,1)
  }
  
  if(cal_eigen == TRUE){
    Lt <- Lx <-  list()
    for (i in 1:n) {
      Lt[[i]] <- t_eval
      Lx[[i]] <- x_func[i,]
    }
    pca_est <- FPCA(Ly = Lx, Lt = Lt, optns= list(dataType = "Dense", error = FALSE, maxK = 2, methodXi = "IN"))
    mu_func <- approxfun(pca_est$workGrid, pca_est$mu)
    eigen_func <- list()
    eigen_func[[1]] <- approxfun(pca_est$workGrid,pca_est$phi[,1])
    eigen_func[[2]] <- approxfun(pca_est$workGrid,pca_est$phi[,2])

    return(list(mu_func = mu_func, eigen_func = eigen_func))
    
  }else{
    
    g <- unlist(lapply(1:n,  FUN = function(i) {control$int_func(x_func[i,], t_eval)}))
    
    return(list(mad_g = mad(g), sd_g = sd(g), g = g))
  }
}



mattern.cov <- function(d, nu=5/2, rho=1, sigma=1) {
  tmp <- sqrt(2*nu)*d/rho
  a <- besselK(x=tmp, nu=nu) #, expon.scaled = FALSE)
  b <- tmp^nu
  return( sigma^2 * 2^(1-nu) * b * a / gamma(nu) )
}


generate_mattern <- function(xi_matrix){
  
  K <- 4 # generate a 4-dim "Mattern"
  # eigenvalues of the covariance function
  #lambdas <- c(0.83, 0.08, 0.029, 0.015)
  lambdas <- c(0.8, 0.3, 0.2, 0.1)
  # eigenfunctions
  rrho <- 3
  sigma <- 1
  nu <- 1/3
  # "compute" eigenfunctions numerically
  x <- seq(0, 1, length=500) 
  mm <- matrix(NA, length(x), length(x))
  for(i in 1:length(x)) mm[i,] <- mattern.cov(abs(x[i]-x), nu=nu, 
                                              rho=rrho, sigma=sigma)
  diag(mm) <- sigma^2
  mm.svd <- svd(mm)$u
  delta <- diff(x)[1]
  phis <- vector('list', K)
  for(j in 1:K) phis[[j]] <- approxfun(x, mm.svd[,j]/sqrt(delta))
  
  mu <- function(tt){
    2*sin(tt*pi)*exp(1-tt)
  }
  
  
  fun <- function(tt){
    xx <- matrix(1,nrow = nrow(xi_matrix), ncol = 1)%*%mu(tt)
    for(j in 1:K){
      xx <- xx + sqrt(lambdas[j]) *as.matrix(xi_matrix[,j])%*%(phis[[j]](tt))
    }
    return(xx)
  }
  
  return(fun)
}



