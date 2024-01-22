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

library(parallel)

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


