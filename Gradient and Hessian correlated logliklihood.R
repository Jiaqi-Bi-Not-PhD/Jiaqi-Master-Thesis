### Correlated log-likelihood newton-raphson step: inner loop gradient ###
### Note that Sigma_inv should be computed within the outer loop ###
gradient <- function(X = X, event = status1, theta = theta, dist = base.dist) {
  expz <- exp(theta[nbase+dim(X):length(theta)])
  grad_beta <- sum(event %*% X, na.rm = TRUE) + sum(- H %*% X %*% expz, na.rm = TRUE)
  grad_z <- sum(event) + sum(-H %*% expz, na.rm = TRUE) - Sigma_inv %*% log(expz)
  grad_alpha <- sum(event/bparms1[1]) + sum(-H %*% expz / bparms[1], na.rm = TRUE)
  grad_lambda <- sum(event * (1/bparms[2] + log(time))) + sum(-H %*% expz * log(time), na.rm = TRUE)
  
  return(list(grad_beta = grad_beta,
              grad_z = grad_z,
              grad_alpha = grad_alpha, 
              grad_lambda = grad_lambda))
}

hessian <- function(X = X, time = time, event = status1, theta = theta, dist = base.dist, Sigma_inv) {
  expz <- exp(theta[(nbase + 1):length(theta)])
  
  hessian_beta <- -t(X) %*% X * H * expz
  hessian_z <- diag(-H * expz) - Sigma_inv
  hessian_alpha <- -sum(event/bparms1[1]^2)
  hessian_lambda <- -sum(event/bparms[2]^2) - sum(H * expz * log(time)^2)
  
  return(list(hessian_beta = hessian_beta,
              hessian_z = hessian_z,
              hessian_alpha = hessian_alpha, 
              hessian_lambda = hessian_lambda))
}

corr_frailty_llhd <- function(X, Y, theta, cuts=NULL, nbase, data, design, base.dist, frailty.dist, agemin, vec=FALSE) {
  if(!design %in% c("pop", "pop+"))  stop("Frailty model is only available for POP or POP+ design.")
  
  nb <- nbase # number of baselines parameters
  
  if(base.dist == "lognormal") bparms1 <- c(theta[1], exp(theta[2]))
  else bparms1 <- exp(theta[1:nbase])
  nx <- dim(X)[2]
  beta <- theta[(nb+1):(nb+nx)]
  xbeta <- c(X%*%theta[(nb+1):(nb+nx)])
  sigma <- exp(theta[nb+nx+1])
  
  
  time <- Y[,1] - agemin
  cuts0 <- cuts - agemin
  status1 <- Y[,2]
  eventsum <- sum(status1)
  print(eventsum)
  
  ip <- data$proband == 1
  
  bhaz <- hazards(base.dist, time, bparms1, cuts=cuts0)
  bcumhaz <- cumhaz(base.dist, time, bparms1, cuts=cuts0)
  
  H <- bcumhaz*exp(xbeta)
  logh <- log(bhaz) + xbeta
  term1 <- status1*logh 
  df1 <- data$df1[ip]
  Hsum <- sum(H, na.rm = TRUE)
  
}

?grad
