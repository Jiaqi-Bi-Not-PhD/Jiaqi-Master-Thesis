##### Log-Normal Frailty Weibull #####
lognormal_single_corr <- function(X, Y, theta, cuts=NULL, nbase, data, design, base.dist, frailty.dist, agemin, vec=FALSE)
{
  
  if(!design %in% c("pop", "pop+"))  stop("Frailty model is only available for POP or POP+ design.")
  
  nb <- nbase # number of baselines parameters
  
  if(base.dist == "lognormal") bparms1 <- c(theta[1], exp(theta[2]))
  else bparms1 <- exp(theta[1:nbase])
  
  nx <- dim(X)[2]
  xbeta <- c(X%*%theta[(nb+1):(nb+nx)])
  sigma <- exp(theta[length(theta)])
  
  
  time <- Y[,1] - agemin
  cuts0 <- cuts - agemin
  status1 <- Y[,2]
  
  ip <- data$proband == 1
  
  bhaz <- hazards(base.dist, time, bparms1, cuts=cuts0)
  bcumhaz <- cumhaz(base.dist, time, bparms1, cuts=cuts0)
  
  H <- bcumhaz*exp(xbeta)
  logh <- log(bhaz) + xbeta
  loglik <- status1*logh 
  df1 <- data$df1[ip]
  Hfam <- aggregate(H, list(data$famID), sum)[,2]
  
  ##### Sort Hfam to run mgauss.hermite for each family
  
  family_list <- split(data, data$famID)
  kinship_matrices <- list()
  for (family_id in names(family_list)) {
    family_data <- family_list[[family_id]]
    kinship_matrix <- with(family_data, kinship(id = indID, dadid = fatherID, momid = motherID, sex = rep(2, nrow(family_data))))
    kinship_matrices[[family_id]] <- kinship_matrix
  }
  
  family_ids <- names(kinship_matrices)
  sorted_families <- kinship_matrices[order(family_ids)]
  
  sorted_Hfam <- Hfam[order(data$famID)]
  sorted_df1 <- df1[order(data$famID)]
  
  #### numerical integral #### 
  gfun <- function(x, L, d, H) {(L %*% x) ^ d %*% exp(-(L %*% x) * H)}
  results <- list()
  
  for(i in seq_along(sorted_families)) {
    family_id <- family_ids[i]
    Sigma <- sorted_families[[family_id]]
    Hfam_value <- sorted_Hfam[i]
    df1_value <- sorted_df1[i]
    L <- chol(Sigma)
    print(L)
    mu <- rep(0, nrow(Sigma))
    pts <- mgauss.hermite(5, mu = mu, sigma = Sigma, prune = 0)
    logdL <- log(apply(pts$points, 1, gfun, d = df1, H = Hfam, L = L) %*% pts$weights)
    
    results[[family_id]] <- logdL
  }
  
  
  # Ascertainment correction (design = pop, pop+)
  
  #cagep <- data$currentage[ip]-agemin
  #xbeta.p <- xbeta[ip]
  #bcumhaz.p <- cumhaz(base.dist, cagep, bparms1, cuts=cuts0)
  
  #H.p <- bcumhaz.p*exp(xbeta.p)
  #logasc <- log(1-laplace(frailty.dist, H.p, sigma))
  #logasc[logasc == -Inf] <- 0
  #slogasc <- sum(logasc[logasc!=-Inf], na.rm=T) 
  #sloglik <- sum(loglik[loglik!=-Inf], na.rm=T) + sum(logdL[logdL!=-Inf], na.rm = T) - slogasc
  
  
  return(-results)
}
X <- as.matrix(data.frame(brca1_prs_cca$mgeneI), 
               nrow=nrow(brca1_prs_cca), 
               ncol = 1)
Y <- as.matrix(data.frame(brca1_prs_cca$timeBC, brca1_prs_cca$BC), 
               nrow = nrow(brca1_prs_cca), 
               ncol = 2)
lognormal_single_corr(X = X, Y = Y, theta = c(1/41.41327,1,0, 0), 
                 nbase = 2, data = brca1_prs_cca,
                 design = "pop", base.dist = "Weibull", frailty.dist = "lognormal",
                 agemin = 18)
