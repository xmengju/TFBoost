#' @import rpart stats
mse <- function(a,b){
  return(mean((a-b)^2))
}

func.l2 <- function(x, cc = NULL) {
  return((x)^2/2)
}

func.l2.grad<- function(x, cc = NULL) {
  return(x)
}

func.lad <- function(x, cc = NULL) {
  return(abs(x))
}

func.lad.grad<- function(x, cc = NULL) {
  return(sign(x))
}


func.tukey <- function(r, cc= 4.685) {
  w <- as.numeric(abs(r) <= cc)
  v <- w*(1 - (1 - (r/cc)^2)^3)  +(1-w)*1
  return(v)
}


func.tukey.grad <- function(r, cc = 4.685) {
  w <- as.numeric(abs(r) <= cc )
  gg <- w*6*r*(1 - (r/cc)^2)^2/(cc^2)  +(1-w)*0
  return(gg)
}

func.huber <- function(r, cc= 0.98) {
  res <- r^2
  res[abs(r) >cc] <- 2*cc*abs(r)[abs(r) >cc] - cc^2
  return (res)
}

func.huber.grad <- function(r, cc = 0.98) {
  res <- r
  res[abs(r) > cc] = sign(r)[abs(r) > cc]*cc
  return(res)
}

init.boosting <- function(loss, init_type = "mean"){
  switch (loss,
          l2 = {
            func <- func.l2
            func.grad <- func.l2.grad
          },
          lad = {
            func <- func.lad
            func.grad <- func.lad.grad
          },
          huber = {
            func <- func.huber
            func.grad <- func.huber.grad
          },
          tukey = {
            func <- func.tukey
            func.grad <- func.tukey.grad
          }
  )
  
  switch(init_type,
          mean = {
            init.func <- function(x){
              mean(x)
            }
          },
          median = {
            init.func <- function(x){
              median(x)
            }
          }
  )
  
  return (list(func = func, func.grad = func.grad, init.func = init.func))
}


cal.neggrad <- function(x1,x2,func.grad){
  
  return(func.grad(x1 - x2))
}

cal.alpha <- function(f_train_t, h_train_t, y_train, func, loss) {

      ff = function(a,r,h){
        return(mean(func(r - a*h)))
      }
       
    if(loss == "l2"){
      return(1)
    }else{
      obj_val <- Inf
      min_val <- 0
      for(upper in c(1,5,10)){
        tmp <- optimize(ff, lower = -1, upper = upper, r = y_train - f_train_t, h = h_train_t)
        if(tmp$objective < obj_val){
          obj_val <- tmp$objective
          min_val <- tmp$minimum
        }
      }
    }
    return(min_val)
}

riemman <- function(uu, tt, rr = c(-1,1)){
  tmp <- sum((uu[-1] + uu[-length(uu)])/2 *(diff(tt)))
  if(rr[1] < tt[1]){
    tmp <- tmp + (tt[1] - rr[1])*uu[1]
  }
  if(rr[2]>tt[length(tt)]){
    tmp <- tmp + (rr[2] - tt[length(tt)])*uu[length(uu)]
  }
  return(tmp)
}
