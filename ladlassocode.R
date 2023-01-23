##package
library(L1pack)
library(ggplot2)
library(caret)
library(quantreg)
library(perryExamples)
library(glmnet) 
library(caret)




rho_a <- function (u, a = 2.68) ifelse(abs(u) < a, abs(u), a)
skipped_median2 <- function (x, a = 2.68, prt = FALSE) {
  # individual intercept : algorithm 2
  m <- median(x)
  sig <- mad(x, center = m)
  iter <- 1; dif <- Inf
  repeat {
    if (prt == TRUE) cat(iter, ",", m, "\n")
    z <- (x-m)/sig
    r <- ifelse(abs(z) > a, z, 0)
    m.new <- median(x[r==0])
    if (abs(m-m.new) < 1e-4 || iter == 200) break
    m <- m.new; iter <- iter + 1
  }
  return(list(m = m, iter = iter))
}


lcad.r <- function(u, a){ifelse(abs(u)>a,u,0)}
lcad.r <- Vectorize(lcad.r)
##code for without lasso
iialcad <- function(y,x.data,a){
  x <- x.data[,-1]
  beta <- lad(y ~ x,data = data.frame(y = y,x = x))$coefficients
  small.r = y - x.data%*%beta
  sigma = mad(small.r, center = skipped_median2(small.r)$m)
  b = 0
  small.r.hat <- c()
  repeat{
    z = (y - x.data%*%beta)/sigma
    small.r.hat <- lcad.r(u = z, a = a)
    po <- which(small.r.hat == 0)
    y <- y[po]
    x.data <- x.data[po,]
    x <- x.data[,-1]
    new.beta <- lad(y ~ x)$coefficients
    judge.loop <-  sqrt(sum((beta - new.beta)^2))
    beta <- new.beta
    b = b+1
    if(judge.loop < 0.001|b == 1000)break
  }
  return(list(iialcad.beta = beta,a = a,b = b))
}



##code with lasso
iialcad_lasso <- function(y, x.data, a){
  x <- x.data[,-1]
  beta <- lad(y ~ x,data = data.frame(y = y,x = x))$coefficients
  small.r = y - x.data%*%beta
  sigma = mad(small.r, center = skipped_median2(small.r)$m)
  b = 0
  small.r.hat <- c()
  repeat{
    z = (y - x.data%*%beta)/sigma
    small.r.hat <- lcad.r(u = z, a = a)
    po <- which(small.r.hat == 0)
    y <- y[po]
    x.data <- x.data[po,]
    x <- x.data[,-1]
    result <- cv.hqreg(x,y,type.measure = "mae")
    best.L.idx <- which(result$lambda == result$lambda.min)
    new.beta <- result$fit$beta[,best.L.idx]
    judge.loop <-  sqrt(sum((beta - new.beta)^2))
    beta <- new.beta
    b = b + 1
    if(judge.loop < 0.001|b == 1000)break
  }
  
  return(list(iialcad.beta = beta,a = a,b = b))
}


LCAD <- function(y, x.data, a, method = "lasso"){
  if(method == "lasso"){
    real_result <- iialcad_lasso(y, x.data, a)
  }else{
    real_result <- iialcad(y, x.data, a)
  }
  return(real_result)
}


