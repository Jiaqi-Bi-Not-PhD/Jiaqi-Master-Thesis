# 1. Data generation
# 2. Likelihood
# 3. Penetrance Functions
# 4. Standard errors for the parameters
# 5. function for the coverage
# 6. Standard errors for the penetrance
# 7. time-dependent C-index (Penetrance for each datasets have to be estiamted first)
#   7.1 Kaplan-Meier Estimators for different time points
#   7.2 All possible pairs
#   7.3 Time-dependent Uno's C-index
# 8. Brier score

# NOTE: Penetrance for each datasets have to be estiamted first

####################################################################
# Parameters
# theta[1:4]:baseline; theta[5:6]: gene effect for cause 1 and cause 2
# theta[7:9]: frailties (k0, k1, k2); theta[10]:beta; theta[11:12]: eta, eta0 for CO

#################################################################### 
# 1. Data generation
# under Permanent Exposure (PE) TVC
cumhaz_pe <- function(a, b, theta) {
  lambda <- theta[1]
  rho <- theta[2]
  beta <- theta[3]
  
  H0 <- function(t) (lambda * t)^rho
  return(exp(beta) * (H0(b) - H0(a)))
}

# under Cox and Oakes (CO) TVC
cumhaz_co <- function(a, b, theta) {
  lambda <- theta[1]
  rho <- theta[2]
  beta <- theta[3]
  eta <- theta[4]
  eta0 <- theta[5]
  
  h0 <- function(t) lambda*rho*(lambda*t)^(rho-1)
  integrand <- function(t) h0(t) * exp((beta * exp(-(t - a) * eta) + eta0) * (t > a))
  int <- try(integrate(integrand, lower = a, upper = b), silent = TRUE)
  if(inherits(int, "try-error")) {
    warning(as.vector(int))
    integrated <- NA_real_
  } else {
    integrated <- int$value
  }
  return(integrated)
}

# under B-spline (BS) TVC
cumhaz_bs <- function(a, b, theta, ik, degree, bndy, intc) {
  lambda <- theta[1]
  rho <- theta[2]
  beta <- theta[3:length(theta)]
  
  h0 <- function(t) lambda*rho*(lambda*t)^(rho-1)
  integrand <- function(t) h0(t) * 
    exp(c(t(bSpline(t-a, knots = ik, degree = degree, Boundary.knots = c(0, bndy),
                    intercept = intc) %*% beta)) * (t > a))
  int <- try(integrate(integrand, lower = a, upper = b), silent = TRUE)
  if(inherits(int, "try-error")) {
    warning(as.vector(int))
    integrated <- NA_real_
  } else {
    integrated <- int$value
  }
  return(integrated)
}

comp_surv_dist <- function(t, base_theta, genebeta1, genebeta2, z1, z2, tvc, 
                           tvcbeta, tt, eta, res, ik = FALSE, degree = FALSE, 
                           bndy = FALSE, intc = FALSE) {
  ind1 <- which(t <= tvc)
  ind2 <- which(t > tvc)
  
  cumhaz1 <- rep(0, length(t))
  
  if (length(ind1) > 0) {
    cumhaz1[ind1] <- (base_theta[1] * t[ind1])^base_theta[2] * exp(genebeta1) * z1
  }
  if (length(ind2) >0) {
    if (tt == "PE") {
      int1 <- cumhaz_pe(a = tvc, b = t[ind2], theta = c(base_theta[1:2], tvcbeta))
    } else if (tt == "CO") {
      int1 <- sapply(t[ind2], cumhaz_co, a = tvc, theta = c(base_theta[1:2],tvcbeta,
                                                            eta[1],eta[2]))
    } else if (tt == "BS") {
      int1 <- sapply(t[ind2], cumhaz_bs, a = tvc, theta = c(base_theta[1:2], tvcbeta),
                     ik = ik, degree = degree, bndy = bndy, intc = intc)
    } 
    cumhaz1[ind2] <- ((base_theta[1] * tvc)^base_theta[2] + int1) * exp(genebeta1) * z1
  }

  cumhaz2 <- (base_theta[3] * t)^base_theta[4] * exp(genebeta2) * z2
  
  return(exp(-cumhaz1 - cumhaz2) - res)
}

comp_inv_surv <- function(base_theta, z1, z2, tvcbeta, tt, eta, val, ik = FALSE, 
                          degree = FALSE, bndy = FALSE, intc = FALSE) {
  uniroot(comp_surv_dist, lower = 0, upper = 250, base_theta = base_theta, 
          z1 = z1, z2 = z2, tvcbeta = tvcbeta, tt = tt, eta = eta, genebeta1 = val[1], 
          genebeta2 = val[2], tvc = val[3], res = val[4], ik = ik, degree = degree, 
          bndy = bndy, intc = intc, extendInt = "yes")$root
}

comp_survp_dist <- function(t, currentage, base_theta, genebeta1, genebeta2, z1, z2, 
                            tvc, tvcbeta, tt, eta, res, ik = FALSE, degree = FALSE,
                            bndy = FALSE, intc = FALSE) {
  A <- comp_surv_dist(t = t, base_theta = base_theta, genebeta1 = genebeta1,
                      genebeta2 = genebeta2, z1 = z1, z2 = z2, tvc = tvc, 
                      tvcbeta = tvcbeta, tt = tt, eta = eta, res = 0, ik = ik, 
                      degree = degree, bndy = bndy, intc = intc)
  B <- comp_surv_dist(t = currentage, base_theta = base_theta, genebeta1 = genebeta1, 
                      genebeta2 = genebeta2, z1 = z1, z2 = z2, tvc = tvc, 
                      tvcbeta = tvcbeta, tt = tt, eta = eta, res = 0, ik = ik, 
                      degree = degree, bndy = bndy, intc = intc)
  return((A - B)/(1 - B) - res) 
}

comp_inv_survp <- function(val, base_theta, z1, z2, tvcbeta, tt, eta, ik = FALSE, 
                           degree = FALSE, bndy = FALSE, intc = FALSE) {
  out <- try(uniroot(comp_survp_dist, lower = 0, upper = 100000, base_theta = base_theta, 
                     z1 = z1, z2 = z2, tvcbeta = tvcbeta, tt = tt, eta = eta, 
                     genebeta1 = val[1], genebeta2 = val[2], currentage = val[3],
                     tvc = val[4], res = val[5], ik = ik, degree = degree,
                     bndy = bndy, intc = intc)$root)
  if(is.null(attr(out, "class"))) {return(out)
  } else {print(c(base_theta, val))}
}

pgene <- function(g, pg, a.freq = 0.0014) {
  qAA <- a.freq^2
  qAa <- 2 * a.freq * (1 - a.freq)
  qaa <- (1 - a.freq)^2
  
  if(length(g) == 1) g <- rep(g, length(pg))
  if(length(pg) == 1) pg <- rep(pg, length(g))
  re <- 0
  re[g == 1] <- cbind(qAA, 1, 0.5, 0, 0.5, 0.25, 0,0,0,0)[pg[g == 1] + 1]
  re[g == 2] <- cbind(qAa, 0, 0.5, 1, 0.5, 0.5, 0.5, 1, 0.5, 0)[pg[g == 2] + 1]
  re[g == 3] <- cbind(qaa, 0,0,0,0, 0.25,0.5, 0, 0.5,1)[pg[g == 3] + 1]
  return(re)
}

parents_g <- function(gene, g, allelefreq){
  ## input gene as a vector of (1,2,3) : AA==1, Aa==2, aa==3
  if(length(allelefreq) == 1) {allelefreq <- c(allelefreq,0)}
  qgene <- allelefreq
  AAq <- qgene^2
  Aaq <- 2 * qgene * (1 - qgene)
  aaq <- (1 - qgene)^2
  
  tmp.g <- matrix(c(1,1, 1,2, 1,3, 2,1, 2,2, 2,3, 3,1, 3,2, 3,3), ncol=2, byrow=T)
  
  tmp.prob<-c(AAq[g]^2, AAq[g]*Aaq[g], AAq[g]*aaq[g], 
              Aaq[g]*AAq[g], Aaq[g]^2, Aaq[g]*aaq[g], 
              aaq[g]*AAq[g], aaq[g]*Aaq[g], aaq[g]^2)
  
  tmp <- sample(1:9, 1, prob = pgene(gene, 1:9) * tmp.prob)
  
  return(tmp.g[tmp,])
}

kids_g <- function(n_kids, p_gene) { 
  ptmp <- sum((p_gene - c(1,0)) * c(3,1))
  return(sample(1:3, n_kids, replace = TRUE, prob = pgene(1:3, ptmp)))
}

haz <- function(t, base_theta) {
  haz <- base_theta[1] * base_theta[2] * (base_theta[1] * t)^(base_theta[2] - 1)
  return(haz)
}

cumhaz <- function(t, base_theta) {
  chaz <- (base_theta[1] * t)^base_theta[2]
  return(chaz)
}


# Penetrance function from 0 to t when t <= tvc
comp_pen1 <- function(t, status, theta, genebeta1, genebeta2) {
  pen <- integrate(comp_integ_pen1, lower = 0, upper = t, status = status,
                   theta = theta, genebeta1 = genebeta1, 
                   genebeta2 = genebeta2)$value
  return(pen)
}

comp_integ_pen1 <- function(u, status, theta, genebeta1, genebeta2) {
  bhaz1 <- haz(u, base_theta = theta[1:2]) * theta[5]
  bhaz2 <- haz(u, base_theta = theta[3:4]) * theta[6]
  
  H1 <- cumhaz(u, base_theta = theta[1:2]) * exp(genebeta1) * theta[5]
  H2 <- cumhaz(u, base_theta = theta[3:4]) * exp(genebeta2) * theta[6]
  integ_pen <- (status == 1) * bhaz1 * exp(genebeta1) * exp(-H1-H2) + 
    (status == 2) * bhaz2 * exp(genebeta2) * exp(-H1-H2)
  return(integ_pen)
}

# Penetrance function from 0 to t when t > tvc
comp_pen2 <- function(t, tvc, status, theta, genebeta1, genebeta2, tvcbeta, tt, 
                      eta = FALSE, ik = FALSE, degree = FALSE, 
                      bndy = FALSE, intc = FALSE) {
  pen <- integrate(comp_integ_pen2, lower = tvc, upper = t, tvc = tvc, status = status,
                   theta = theta, genebeta1 = genebeta1, genebeta2 = genebeta2,
                   tvcbeta = tvcbeta, tt = tt, eta = eta, ik = ik, degree = degree,
                   bndy = bndy, intc = intc)$value
  return(pen)
}

comp_integ_pen2 <- function(u, tvc, status, theta, genebeta1, genebeta2, tvcbeta, tt, 
                            eta = FALSE, ik = FALSE, degree = FALSE, 
                            bndy = FALSE, intc = FALSE) {
  bhaz1 <- haz(u, base_theta = theta[1:2]) * theta[5]
  bhaz2 <- haz(u, base_theta = theta[3:4]) * theta[6]
  
  if (tt == "PE") {
    tvcb1 <- tvcbeta
    H1 <- (cumhaz(u, base_theta = theta[1:2]) - cumhaz(tvc, base_theta = theta[1:2])) *
      exp(genebeta1 + tvcbeta) * theta[5]
    H2 <- cumhaz(u, base_theta = theta[3:4]) * exp(genebeta2) * theta[6]
  } else if (tt == "CO") {
    tvcb1 <- tvcbeta * exp(- (u - tvc) * eta[1]) + eta[2]
    H1 <- exp(genebeta1) * theta[5] *
      (cumhaz(tvc, base_theta = theta[1:2]) + 
         sapply(u, cumhaz_co, a = tvc, theta = c(theta[1:2], tvcbeta, eta)))
    H2 <- cumhaz(u, base_theta = theta[3:4]) * exp(genebeta2) * theta[6]
  } else if (tt == "BS") {
    bs <- bSpline(u - tvc, knots = ik, degree = degree,
                  Boundary.knots = c(0, bndy), intercept = intc)
    tvcb1 <- c(t(bs %*% tvcbeta))
    H1 <- exp(genebeta1) * theta[5] *
      (cumhaz(tvc, base_theta = theta[1:2]) + 
         sapply(u, cumhaz_bs, a = tvc, theta = c(theta[1:2], tvcbeta), 
                ik = ik, degree = degree, bndy = bndy, intc = intc))
    H2 <- cumhaz(u, base_theta = theta[3:4]) * exp(genebeta2) * theta[6]
  }
   
  integ_pen <- (status == 1) * bhaz1 * exp(genebeta1 + tvcb1) * exp(-H1-H2) + 
    (status == 2) * bhaz2 * exp(genebeta2) * exp(-H1-H2)
  return(integ_pen)
}

comp_fgene <- function(affage, tvc, variation, parms, vbeta, alpha1, alpha2, 
                       pg=0, m_carrier=0, dominant.m=TRUE, aq, tt, ik = FALSE, 
                       degree = FALSE, bndy = FALSE, intc = FALSE){
  # returns 2x3 matrix , first row for the major gene, 
  #	second row for the second gene
  pAA <- pAa <- paa <- 0
  AAq <- Aaq <- aaq <- 0
  AAq[1] <- pgene(1, pg=pg[1], a.freq=aq[1])
  Aaq[1] <- pgene(2, pg=pg[1], a.freq=aq[1])
  aaq[1] <- pgene(3, pg=pg[1], a.freq=aq[1])
  
  av = c(alpha1,alpha2)
  if(tt=="PE") {
    bz = vbeta[3]
    eta = FALSE  
  } else if(tt=="CO") {
    bz = vbeta[3]
    eta = vbeta[4:5]
  } else if(tt=="BS") {
    bz = vbeta[3:length(vbeta)]
    eta = FALSE
  } 
  
  Ft <- 0
  for(i in c(0,1)){
    xbeta1 <- i*vbeta[1] #+alpha1
    xbeta2 <- i*vbeta[2] #+alpha2
    Ft[i+1] <- 1 - comp_surv_dist(affage, base_theta=parms, genebeta1=xbeta1, genebeta2=xbeta2, 
                                  z1=alpha1, z2=alpha2, tvc=tvc, 
                                  tvcbeta=bz, tt=tt, eta=eta, res=0, 
                                  ik = ik, degree = degree, bndy = bndy, intc = intc) 
    
#    if(affage <= tvc){
#      Ft[i+1] <- comp_pen1(affage, status = 1, theta = c(parms,av), 
#                           genebeta1 = xbeta1, genebeta2 = xbeta2) + 
#        comp_pen1(affage, status = 2, theta = c(parms,av), 
#                  genebeta1 = xbeta1, genebeta2 = xbeta2)
#    }else{
#      Ft[i+1] <- comp_pen1(tvc, status = 1, theta = c(parms,av), 
#                           genebeta1 = xbeta1, genebeta2 = xbeta2) + 
#        comp_pen2(affage, tvc, status = 1, theta = c(parms,av), genebeta1 = xbeta1, 
#                  genebeta2 = xbeta2, tvcbeta = bz, tt=tt, eta=eta, ik = ik, 
#                  degree = degree, bndy = bndy, intc = intc) + 
#        comp_pen1(affage, status = 2, theta = c(parms,av), 
#                  genebeta1 = xbeta1, genebeta2 = xbeta2)
#    }
#  }
  
  pAA <- Ft[2]*AAq[1]
  if(dominant.m) pAa <- Ft[2]*Aaq[1]
  else pAa <- Ft[1]*Aaq[1]
  paa <- Ft[1]*aaq[1]
  
  return(cbind(pAA, pAa, paa)/c(pAA + pAa + paa))
}

# Generate the data for the simulation
comp_familyStructure <- function(i, cumind, m_carrier, depend, vbeta, base_theta, 
                                 dominant.m = TRUE, dominant.s = TRUE, allelefreq, mrate,
                                 agemin, page, tvcage, tvctype, frailty, ik = FALSE, 
                                 degree = FALSE, bndy = FALSE, intc = FALSE) {
  data_i <- numeric()
  indID <- c(cumind + 1, cumind + 2)
  motherID <- c(0,0)
  fatherID <- c(0,0)
  gender <- c(1,0) # 0 - female, 1 - male
  proband <- c(0,0) # 0 - non-proband, 1 - proband
  generation <- c(1,1) # generation number, 0 in second generation means married from outside
  relation <- c(4,4) # 1 = proband, 2 = sibling, 3 = child, 4 = parent, 5 = sibchild, 
  # 6 = husband, 7 = sibspouse  
  
  # number of siblings in second generation be 2 to 5
  # truncated negative binomial (>=2 and <=5) with prob=sibprob=0.4551, and rr=1
  sec_num <- sample(c(2,3,4,5), 1, replace = TRUE, prob = c(0.4991, 0.2720, 0.1482, 0.0807))
  
  # generate gender by prob = 0.5
  sec_gender <- sample(c(1,0), sec_num, replace = TRUE, prob = c(0.5, 0.5))
  sec_num <- length(sec_gender)
  sec_proband <- sample(sec_num, 1, replace = TRUE)
  num_mem <- 2 * sec_num + cumind + 2
  
  thi_num <- rep(0, sec_num)
  for (j in 1:sec_num) {
    if (j == sec_proband) {
      proband <- c(proband, c(1,0))
      relation <- c(relation, c(1,6))
    } else {
      proband <- c(proband, c(0,0))
      relation <- c(relation, c(2,7))
    }
    
    indID <- c(indID, 2 * j + cumind + 1, 2 * j + cumind + 2)
    fatherID <- c(fatherID, c(cumind + 1, 0))
    motherID <- c(motherID, c(cumind + 2, 0))
    
    if (sec_gender[j] == 0) {
      gender <- c(gender, c(0,1))
    } else {
      gender <- c(gender, c(1,0))
    }
    
    generation <- c(generation, c(2,0))
    
    thi_num[j] <- sample(c(0, 1, 2), 1, replace = TRUE, prob = c(0.1089, 0.3420, 0.5491))
    #thi_num <- sample(c(2,3,4,5), 1, replace = TRUE, prob = c(0.4991, 0.2720, 0.1482, 0.0807))
    #thi_num <- sample(c(2,3,4), 1, replace = TRUE, prob = c(0.4991, 0.2720, 0.2289))
    
    if (thi_num[j] != 0) {
      for (k in 1:thi_num[j]) {
        proband <- c(proband, 0)
        indID <- c(indID, num_mem + k)
        
        if (j == sec_proband) {
          relation <- c(relation, 3)
        } else {
          relation <- c(relation, 5)
        }
        
        if (gender[indID == (2 * j + cumind + 1)] == 0) {
          fatherID <- c(fatherID, 2 * j + cumind + 2)
          motherID <- c(motherID, 2 * j + cumind + 1)
        } else {
          fatherID <- c(fatherID, 2 * j + cumind + 1)
          motherID <- c(motherID, 2 * j + cumind + 2)
        }
        
        # generate gedner
        gender <- c(gender, sample(c(1,0), 1, replace = TRUE, prob = c(0.5, 0.5)))
        generation <- c(generation, 3)
      }
    }
    
    num_mem <- num_mem + thi_num[j]
  }
  
  num_ind_i <- length(indID)
  famID <- rep(i, num_ind_i)
  ageonset <- rep(0, num_ind_i)
  censor_age <- rep(0, num_ind_i)
  status <- rep(0, num_ind_i)
  tvc <- rep(0, num_ind_i)
  affected <- rep(0, num_ind_i)
  disgene.m <- rep(0, num_ind_i)
  disgene.s <- rep(0, num_ind_i)
  parentsG.m <- rep(0, num_ind_i) # parents genotype: it can be 1 to 9, 0 if founders
  pos <- c(1:num_ind_i)
  
  # -----------------------------------------------------------
  ### Familial correlation is generated by gamma frailty
  if (frailty == "gamma") {
    z0 <- 0
  } else if (frailty == "cgamma") {
    z0 <- rgamma(1, shape = depend[1], scale=1/depend[1])
  }
  #z0 <- rgamma(1, shape = depend[1], scale=1/depend[1])
  z1 <- rgamma(1, shape = depend[2], scale = 1/(depend[1] + depend[2])) + 
    z0 * depend[1]/(depend[1] + depend[2])
  z2 <- rgamma(1, shape = depend[3], scale = 1/(depend[1] + depend[3])) + 
    z0 * depend[1]/(depend[1] + depend[3])
  
  ### Generating time-varying covariate
  tvc <- rnorm(length(pos), mean = tvcage-agemin, sd = 2.5) #mean=15
  
  ### Generating the current age
  # Generating the current age of the proband first
  gen_pos <- pos[proband == 1]
  a_f <- rnorm(1, mean = page, sd = 5)
  censor_age[gen_pos] <- a_f
  
  gen_pos <- pos[generation == 1]
  censor_age[gen_pos] <- rnorm(length(gen_pos), mean = (a_f + 20), sd = 2.5)
  min_age1<- min(censor_age[gen_pos])
  if(min_age1 < agemin) stop("agemin is too large.")
  
  gen_pos <- pos[generation == 0]
  censor_age[gen_pos] <- rnorm(length(gen_pos), mean = a_f, sd = 2.5)
  
  gen_pos <- pos[generation == 2]
  for (j in 1:length(gen_pos)) {
    if (gen_pos[j] != pos[proband == 1]) {
      sec_age <- rtruncnorm(1, a = agemin, b = min_age1-14, mean = page, sd = 2.5) # use page instead of a_f
      censor_age[gen_pos[j]] <- sec_age
    }
    
    if (thi_num[j] != 0) {
      thi_pos <- pos[fatherID == indID[gen_pos[j]] | motherID == indID[gen_pos[j]]]
      min_age2 <- min(censor_age[indID == fatherID[thi_pos[1]] | 
                                   indID == motherID[thi_pos[1]]])
      for (k in 1:length(thi_pos)) {
        censor_age[thi_pos[k]] <- rtruncnorm(1, a = agemin, b = min_age2, 
                                             mean = (min_age2 - 20), sd = 2.5)
        #censor_age[thi_pos[k]] <- rnorm(1, mean = min_age2 - 20 - k, sd = 1.5) 
      }
    }
  }
  
  # -----------------------------------------------------------
  # Generating the genotypes
  if (length(allelefreq) == 1) {allelefreq <- c(allelefreq, 0)}
  qgene <- allelefreq
  AAq <- qgene^2
  Aaq <- 2 * qgene * (1 - qgene)
  aaq <- (1 - qgene)^2
  G <- cbind(AAq, Aaq, aaq)
  
  # Generating the proband's genotype first given its age at onset and gender
  prob_age <- censor_age[proband == 1]
  prob_gender <- gender[proband == 1]
  prob_tvc <- tvc[proband == 1]
  prob_ID <- indID[proband == 1]
  
  prob_gene <- comp_fgene(affage = prob_age - agemin, tvc = prob_tvc, variation = "frailty",
                          parms = base_theta, vbeta = vbeta, alpha1 = z1, alpha2 = z2, pg = 0,
                          m_carrier = m_carrier, dominant.m = dominant.m, aq = allelefreq, 
                          tt = tvctype, ik = ik, degree = degree, bndy = bndy, intc = intc)
  
  if (any(is.na(prob_gene))) {
    return(data.frame())
  } else{
    
    if (m_carrier == 1) {
      if (dominant.m) {
        gg <- 1:2
      } else {
        gg <- 1
      }
      prob_G <- sample(gg, 1, replace = TRUE, prob = prob_gene[1,gg])
    } else {
      prob_G <- sample(c(1,2,3), 1, replace = TRUE, prob = prob_gene[1,])
    }
    
    # generating proband's parents' genotype based on proband's genotype 
    G1 <- parents_g(prob_G, g = 1, allelefreq = allelefreq)
    
    # generating kids' genotype given the proband's genotype
    n_sibs <- sum(proband == 0 & generation == 2)
    sib_G <- kids_g(n_sibs, G1)
    
    disgene.m[generation == 1] <- G1 
    disgene.m[proband == 1] <- prob_G
    disgene.m[proband == 0 & generation == 2] <- sib_G
    
    # genotypes of first and second gene for generation=0 which is founder by 
    # random selection among AA, Aa, aa 
    disgene.m[generation == 0] <- sample(c(1,2,3), sum(generation == 0), replace=TRUE, 
                                         prob = G[1,])
    
    # Generating genotypes of third generation
    if (sum(thi_num) != 0) {
      for (j in indID[generation == 3]) {
        m.g <- disgene.m[indID == motherID[indID == j]]
        f.g <- disgene.m[indID == fatherID[indID == j]]
        disgene.m[indID == j] <- kids_g(1, c(m.g, f.g))
      }
    }
    
    if(dominant.m) {majorgene <- ifelse(disgene.m == 3, 0, 1)
    } else {majorgene <- ifelse(disgene.m == 1, 1, 0)}
    
    # -----------------------------------------------------------
    # simulate the age at on set for the family member
    vbeta1 <- vbeta[1]   ## log relative risk of the major gene for cause1
    vbeta2 <- vbeta[2]   ## log relative risk of the major gene for cause2
    
    ## log relative risk of the beta for tvc for cause1
    if(tvctype == "BS") {
      tvcbeta <- vbeta[3:length(vbeta)]
    } else if (tvctype == "CB") {
      tvcbeta <- vbeta[3:6]
    } else {
      tvcbeta <- vbeta[3]
    } 
    
    # gen = generation 1,2,3
    gen <- ifelse(generation == 2 | generation == 0, 2, ifelse(generation == 1, 1, 3))
    affected <- (proband == 1)
    
    genebeta1 <- majorgene * vbeta1
    genebeta2 <- majorgene * vbeta2
    
    if (tvctype == "CO") {
      eta <- vbeta[4:5]
    } else {
      eta <- FALSE
    } 
    
    # -----------------------------------------------------------
    ## generate ageonset
    uni <- runif(num_ind_i, 0, 1)
    
    ageonset <- apply(cbind(genebeta1, genebeta2, tvc, uni), 1, 
                      comp_inv_surv, base_theta = base_theta, z1 = z1, z2 = z2, 
                      tvcbeta = tvcbeta, tt = tvctype, eta = eta, ik = ik, 
                      degree = degree, bndy = bndy, intc = intc) + agemin
    
    t <- ageonset - agemin
    itx <- (t > tvc)
    tvcb1 <- rep(0, num_ind_i)
    if (sum(itx) > 0) {
      tx <- (t - tvc)[itx]
      if (tvctype == "PE") {
        tvcb1[itx] <- tvcbeta
      } else if (tvctype == "CO") {
        tvcb1[itx] <- tvcbeta * exp(- tx * eta[1]) + eta[2] 
      } else if (tvctype == "BS") {
        bs <- bSpline(tx, knots = ik, degree = degree, 
                      Boundary.knots = c(0, bndy), intercept = intc)
        tvcb1[itx] <- c(t(bs %*% tvcbeta))
      } 
    }
    
    h1 <- z1 * haz(t, base_theta[1:2]) * exp(genebeta1 + tvcb1)
    h2 <- z2 * haz(t, base_theta[3:4]) * exp(genebeta2)
    true_status <- ifelse(runif(num_ind_i) < h1/(h1 + h2), 1, 2) 
    
    # currentage <- ifelse(censor_age > 100, 100, censor_age)
    currentage <- censor_age
    time <- pmin(currentage, ageonset)
    true_time <- ageonset # no censoring
    status <- ifelse(currentage >= ageonset, true_status, 0)
    
    tvc <- tvc + agemin
    tvc_status <- ifelse(time > tvc, 1, 0)
    true_tvc_status <- ifelse(true_time > tvc, 1, 0)
    
    mgene <- majorgene
    fsize <- length(famID)
    naff <- sum(status != 0)
    true_naff <- sum(true_status != 0)
    data_i <- cbind(famID, indID, gender, motherID, fatherID, proband, generation,
                    majorgene = disgene.m, secondgene = disgene.s, ageonset, currentage,
                    time, status, mgene, relation, fsize, naff, tvc, true_status, 
                    tvc_status, true_tvc_status, true_naff)
    return(data_i)
  }
}

comp_familyDesign <- function(n, m_carrier, dominant.m = TRUE, dominant.s = TRUE, 
                              depend, vbeta, base_theta, allelefreq, mrate, agemin, page, tvcage,
                              tvctype, frailty, ik = FALSE, degree = FALSE, 
                              bndy = FALSE, intc = FALSE) {
  data <- numeric()
  cumind <- 0
  i <- 1
  j <- 0
  
  while (i <= n) {
    j <- j + 1
    data_i <- comp_familyStructure(i, cumind = cumind, m_carrier = m_carrier,
                                   depend = depend, vbeta = vbeta, base_theta = base_theta,
                                   dominant.m = TRUE, dominant.s = TRUE,
                                   allelefreq = c(0.0021, 0.2), mrate = mrate, 
                                   agemin = agemin, page=page, tvcage=tvcage, tvctype = tvctype, frailty = frailty,
                                   ik = ik, degree = degree, bndy = bndy, intc = intc)
    if (is.null(attr(data_i, "class"))) {
      until <- (data_i[data_i[,"proband"]==1, "status"] == 1 | 
                   data_i[data_i[,"proband"]==1, "status"] == 2)
      
      if (!is.null(dim(data_i))) {
        if (nrow(data_i) > 0) {
          if (until) {
            data <- rbind(data, data_i)
            cumind <- cumind + nrow(data_i)
            i <- i + 1
          }
        }
      }
    }
  }
  return(data)
}

comp_simfam <- function(n, theta, tvctype, page, tvcage, frailty, ik = FALSE, degree = FALSE, bndy = FALSE,
                        intc = FALSE) {
  # vbeta: frailties, beta,..
  if (tvctype == "PE" | tvctype == "SN" | tvctype == "SIN") {vbeta <- c(theta[5:6], theta[10])}
  else if (tvctype == "CO") {vbeta <- c(theta[5:6], theta[10], exp(theta[11]), theta[12])}
  else if (tvctype == "BS") {vbeta <- c(theta[5:6], theta[10:length(theta)])}
  else if (tvctype == "CB") {vbeta <- c(theta[5:6], theta[10:13])}
  #else if (tvctype == "SN") {vbeta <- c(theta[5:6], theta[10])}
  
  simdata = comp_familyDesign(n = n, m_carrier = 1, depend = exp(theta[7:9]), 
                              vbeta = vbeta, base_theta = exp(theta[1:4]), 
                              allelefreq = 0.0021, mrate = 0, agemin = 16, page=page, tvcage=tvcage,
                              tvctype = tvctype, frailty = frailty, ik = ik, degree = degree,
                              bndy = bndy, intc = intc)
  
  simdata <- data.frame(simdata)
  simdata <- simdata[simdata$time > 16,] # include age > 16
  names(simdata)[1:5] <- c("famID", "indID", "gender", "moid", "faid")
  
  fsize <- aggregate(simdata$famID, by = list(simdata$famID), length)[,2]
  df1 <- aggregate(simdata$status == 1, by = list(simdata$famID), sum)[,2]
  df2 <- aggregate(simdata$status == 2, by = list(simdata$famID), sum)[,2]
  simdata$fsize <- rep(fsize, fsize) # family size
  simdata$df1 <- rep(df1, fsize) # number of subjects with event 1
  simdata$df2 <- rep(df2, fsize) # number of subjects with event 2
  return(simdata)
}

######################################
# Sampling n families
sampled_fam_fun <- function(data, num_fam, replicates) {
  sampled_fam <- matrix(NA, nrow = num_fam, ncol = replicates)
  for (i in 1:replicates) {
    sampled_fam[,i] <- sample(unique(data$famID), size = num_fam, replace = FALSE)
  }
  return(sampled_fam)
}

#################################################################### 
# 2. Likelihood
fp_fun <- function(x1, x2, k0, k1, k2, H1,H2,d1, d2){
  gamma(d1+1)*gamma(d2+1)/(gamma(k0)*gamma(k1)*gamma(k2))*(k0+k1)^(-d1)*(k0+k2)^(-d2)*
    gamma(k0+x1+x2)*(1+H1/(k0+k1)+H2/(k0+k2))^(-k0-x1-x2)*
    gamma(k1+d1-x1)*(1+H1/(k0+k1))^(-k1-d1+x1)/gamma(x1+1)/gamma(d1-x1+1)*
    gamma(k2+d2-x2)*(1+H2/(k0+k2))^(-k2-d2+x2)/gamma(x2+1)/gamma(d2-x2+1)
}
fp_fun1 <- function(vec, x2, k0, k1, k2) {
  sapply(0:vec[1], fp_fun, x2=x2, k0=k0, k1=k1, k2=k2, H1=vec[3], H2=vec[4], 
         d1=vec[1], d2=vec[2])
}

fp_fun2 <- function(vec, k0, k1, k2) {
  sum(sapply(0:vec[2], fp_fun1, k0=k0, k1=k1, k2=k2, vec=vec))
}


cumhaz_i <- function(ab, theta, tvctype, ik = FALSE, degree = FALSE, bndy = FALSE,
                     intc = FALSE) {
  a <- ab[1]
  b <- ab[2]
  if (tvctype == "PE") {
    cumhaz <- cumhaz_pe(a, b, theta)
  } else if (tvctype == "CO") {
    cumhaz <- cumhaz_co(a, b, theta)
  } else if (tvctype == "BS") {
    cumhaz <- cumhaz_bs(a, b, theta, ik, degree, bndy, intc)
  } else if (tvctype == "CB") {
    cumhaz <- cumhaz_cb(a, b, theta)
  } else if (tvctype == "SN") {
    cumhaz <- cumhaz_sn(a, b, theta)
  } else if (tvctype == "SIN") {
    cumhaz <- cumhaz_sin(a, b, theta)
  }
  return(cumhaz)
}

cumhaz_tvc <- function(a, b, theta, tvctype, ik = FALSE, degree = FALSE, 
                       bndy = FALSE, intc = FALSE) {
  unlist(lapply(1:length(a), function(i) 
    cumhaz_i(ab = c(a[i], b[i]), theta = theta, tvctype = tvctype, ik = ik, 
             degree = degree, bndy = bndy, intc = intc)))
}

comp_llik_sim <- function(theta, data, agemin, tvctype, frailty, ik = FALSE, 
                          degree = FALSE, bndy = FALSE, intc = FALSE) {
  data = data[data$currentage >= agemin,]
  
  # base parameters
  lambda1 = exp(theta[1])
  rho1 = exp(theta[2])
  lambda2 = exp(theta[3])
  rho2 = exp(theta[4])
  
  # vbeta for gene 
  beta_gen1 = theta[5]
  beta_gen2 = theta[6]
  
  # correlated gamma frailties
  if (frailty == "gamma") {
    if (tvctype == "PE") {
      beta_tvc <- theta[9]
    } else if (tvctype == "CO") {
      beta_tvc <- theta[9]
      eta <- exp(theta[10])
      eta0 <- theta[11]
    } else if (tvctype == "BS") {
      beta_tvc <- theta[9:length(theta)]
    } 
    
    k0 <- 0
    k1 <- exp(theta[7])
    k2 <- exp(theta[8])
    w1 <- k0 + k1
    w2 <- k0 + k2
  } else if (frailty == "cgamma") {
    if (tvctype == "PE") {
      beta_tvc <- theta[10]
    } else if (tvctype == "CO") {
      beta_tvc <- theta[10]
      eta <- exp(theta[11])
      eta0 <- theta[12]
    } else if (tvctype == "BS") {
      beta_tvc <- theta[10:length(theta)]
    }
    
    k0 <- exp(theta[7])
    k1 <- exp(theta[8])
    k2 <- exp(theta[9])
    w1 <- k0 + k1
    w2 <- k0 + k2
  }
  
  time0 <- data$time - agemin
  status <- data$status
  mgene <- data$mgene
  tvc <- data$tvc - agemin
  itvc <- data$tvc_status
  ip <- which(data$proband == 1)
  famID <- data$famID
  
  # number of individuals in the unique family
  df1 <- data$df1[!duplicated(data$famID)] # 1
  df2 <- data$df2[!duplicated(data$famID)] # 2
  
  bhaz1 <- lambda1 * rho1 * (lambda1 * time0)^(rho1 - 1)
  bhaz2 <- lambda2 * rho2 * (lambda2 * time0)^(rho2 - 1)
  
  logh1 <- log(bhaz1) + mgene * beta_gen1
  logh2 <- log(bhaz2) + mgene * beta_gen2
  
  if (tvctype == "PE") {
    tvc_vec <- itvc * beta_tvc
    cumhaz1 <- cumhaz_tvc(tvc[itvc == 1], time0[itvc == 1], 
                          c(lambda1, rho1, beta_tvc), "PE")
  } else if (tvctype == "CO") {
    tvc_vec <- itvc * (beta_tvc * exp(- eta * (time0 - tvc)) + eta0)
    cumhaz1 <- cumhaz_tvc(tvc[itvc == 1], time0[itvc == 1], 
                          c(lambda1, rho1, beta_tvc, eta, eta0), "CO")
  } else if (tvctype == "BS") {
    # by including intercept, bSpline gives 1 when tx equals to 0
    if (intc == TRUE) {
      bs <- matrix(0, nrow = nrow(data), ncol = length(ik) + degree + 1)
    } else {
      bs <- matrix(0, nrow = nrow(data), ncol = length(ik) + degree)
    }
    tx <- time0[itvc == 1] - tvc[itvc == 1]
    bs[itvc == 1,] <- bSpline(tx, knots = ik, degree = degree, 
                              Boundary.knots = c(0, bndy), intercept = intc)
    tvc_vec <- itvc * c(t(bs %*% beta_tvc))
    cumhaz1 <- cumhaz_tvc(tvc[itvc == 1], time0[itvc == 1], c(lambda1, rho1, beta_tvc),
                          "BS", ik = ik, degree = degree, bndy = bndy, intc = intc)
  } 
  
  llik1 <- - sum((logh1 + tvc_vec)[status == 1], na.rm = TRUE)
  llik2 <- - sum(logh2[status == 2], na.rm = TRUE)
  
  H1 <- (lambda1 * time0)^rho1 * exp(beta_gen1 * mgene)
  H1[itvc == 1] <- ((lambda1 * tvc[itvc==1])^rho1 + cumhaz1) *
    exp(beta_gen1 * mgene[itvc == 1])
  H2 <- (lambda2 * time0)^rho2 * exp(beta_gen2 * mgene)
  Hfam1 <- aggregate(H1, by = list(famID), FUN = sum)[,2]
  Hfam2 <- aggregate(H2, by = list(famID), FUN = sum)[,2]
  
  if (frailty == "gamma") {
    llik3 <- -sum(lfactorial(k1+df1-1)-(df1-1)*log(k1)-lfactorial(k1) + (-k1-df1)*log(1+Hfam1/k1), na.rm=T) -
      sum(lfactorial(k2+df2-1)-(df2-1)*log(k2)-lfactorial(k2) + (-k2-df2)*log(1+Hfam2/k2), na.rm=T)
  } else if (frailty == "cgamma") {
    llik3 <- log(apply(cbind(df1, df2, Hfam1, Hfam2), 1, fp_fun2, k0=k0, k1=k1, k2=k2))
    llik3 <- - sum(llik3, na.rm = TRUE)
  }
  
  llik <- llik1 + llik2 + llik3
  
  # Ascertainment correction by design="pop+"
  cage_p <- data$currentage[ip] - agemin
  status_p <- data$status[ip]
  tvc_p <- tvc[ip]
  itvc_p <- ifelse(tvc_p < cage_p, 1, 0)
  
  if (tvctype == "PE") {
    cumhaz1_p <- cumhaz_tvc(tvc_p[itvc_p == 1], cage_p[itvc_p == 1], 
                            c(lambda1, rho1, beta_tvc), "PE")
  } else if (tvctype == "CO") {
    cumhaz1_p <- cumhaz_tvc(tvc_p[itvc_p == 1], cage_p[itvc_p == 1], 
                            c(lambda1, rho1, beta_tvc, eta, eta0), "CO")
  } else if (tvctype == "BS") {
    cumhaz1_p <- cumhaz_tvc(tvc_p[itvc_p == 1], cage_p[itvc_p == 1], 
                            c(lambda1, rho1, beta_tvc), "BS", ik = ik, degree = degree,
                            bndy = bndy, intc = intc)
  }
  
  H1_p <- (lambda1 * cage_p)^rho1 * exp(beta_gen1)
  H1_p[itvc_p == 1] <- ((lambda1 * tvc_p[itvc_p == 1])^rho1 + cumhaz1_p) * exp(beta_gen1)
  H2_p <- (lambda2 * cage_p)^rho2 * exp(beta_gen2)
  llik_p <- - sum(log(1 - (1 + H1_p/w1 + H2_p/w2)^(-k0) * 
                        (1 + H1_p/w1)^(-k1)*(1 + H2_p/w2)^(-k2)), na.rm = TRUE)
  
  # Final value of the log likelihood
  final_llik <- llik - llik_p 
  return(final_llik)
}

# vector form of the likelihood
comp_llik_sim_vec <- function(theta, data, agemin, tvctype, frailty, ik = FALSE, degree = FALSE,
                              bndy = FALSE, intc = FALSE) {
  data = data[data$currentage >= agemin,]
  
  # base parameters
  lambda1 = exp(theta[1])
  rho1 = exp(theta[2])
  lambda2 = exp(theta[3])
  rho2 = exp(theta[4])
  
  # vbeta for gene 
  beta_gen1 = theta[5]
  beta_gen2 = theta[6]
  
  if (frailty == "gamma") {
    if (tvctype == "PE") {
      beta_tvc <- theta[9]
    } else if (tvctype == "CO") {
      beta_tvc <- theta[9]
      eta <- exp(theta[10])
      eta0 <- theta[11]
    } else if (tvctype == "BS") {
      beta_tvc <- theta[9:length(theta)]
    } 
    
    k0 <- 0
    k1 <- exp(theta[7])
    k2 <- exp(theta[8])
    w1 <- k0 + k1
    w2 <- k0 + k2
  } else if (frailty == "cgamma") {
    if (tvctype == "PE") {
      beta_tvc <- theta[10]
    } else if (tvctype == "CO") {
      beta_tvc <- theta[10]
      eta <- exp(theta[11])
      eta0 <- theta[12]
    } else if (tvctype == "BS") {
      beta_tvc <- theta[10:length(theta)]
    }
    
    k0 <- exp(theta[7])
    k1 <- exp(theta[8])
    k2 <- exp(theta[9])
    w1 <- k0 + k1
    w2 <- k0 + k2
  }
  
  time0 <- data$time - agemin
  status <- data$status
  mgene <- data$mgene
  tvc <- data$tvc - agemin
  itvc <- data$tvc_status
  ip <- which(data$proband == 1)
  famID <- data$famID
  
  # number of individuals in the unique family
  df1 <- data$df1[!duplicated(data$famID)] # 1
  df2 <- data$df2[!duplicated(data$famID)] # 2
  
  bhaz1 <- lambda1 * rho1 * (lambda1 * time0)^(rho1 - 1)
  bhaz2 <- lambda2 * rho2 * (lambda2 * time0)^(rho2 - 1)
  
  logh1 <- log(bhaz1) + mgene * beta_gen1
  logh2 <- log(bhaz2) + mgene * beta_gen2
  
  if (tvctype == "PE") {
    tvc_vec <- itvc * beta_tvc
    cumhaz1 <- cumhaz_tvc(tvc[itvc == 1], time0[itvc == 1], 
                          c(lambda1, rho1, beta_tvc), "PE")
  } else if (tvctype == "CO") {
    tvc_vec <- itvc * (beta_tvc * exp(- eta * (time0 - tvc)) + eta0)
    cumhaz1 <- cumhaz_tvc(tvc[itvc == 1], time0[itvc == 1], 
                          c(lambda1, rho1, beta_tvc, eta, eta0), "CO")
  } else if (tvctype == "BS") {
    # by including intercept, bSpline gives 1 when tx equals to 0
    if (intc == TRUE) {
      bs <- matrix(0, nrow = nrow(data), ncol = length(ik) + degree + 1)
    } else {
      bs <- matrix(0, nrow = nrow(data), ncol = length(ik) + degree)
    }
    tx <- time0[itvc == 1] - tvc[itvc == 1]
    bs[itvc == 1,] <- bSpline(tx, knots = ik, degree = degree, 
                              Boundary.knots = c(0, bndy), intercept = intc)
    tvc_vec <- itvc * c(t(bs %*% beta_tvc))
    cumhaz1 <- cumhaz_tvc(tvc[itvc == 1], time0[itvc == 1], c(lambda1, rho1, beta_tvc),
                          "BS", ik = ik, degree = degree, bndy = bndy, intc = intc)
  } 
  
  llik1 <- (logh1 + tvc_vec) * (status == 1)
  llik1 <- - aggregate(llik1, by = list(famID), FUN = sum)[,2]
  llik2 <- logh2  * (status == 2)
  llik2 <- - aggregate(llik2, by = list(famID), FUN = sum)[,2]
  
  H1 <- (lambda1 * time0)^rho1 * exp(beta_gen1 * mgene)
  H1[itvc == 1] <- ((lambda1 * tvc[itvc==1])^rho1 + cumhaz1) *
    exp(beta_gen1 * mgene[itvc == 1])
  H2 <- (lambda2 * time0)^rho2 * exp(beta_gen2 * mgene)
  Hfam1 <- aggregate(H1, by = list(famID), FUN = sum)[,2]
  Hfam2 <- aggregate(H2, by = list(famID), FUN = sum)[,2]
  
  if (frailty == "gamma") {
    llik3 <- -(lfactorial(k1+df1-1)-(df1-1)*log(k1)-lfactorial(k1) + (-k1-df1)*log(1+Hfam1/k1)) -
      (lfactorial(k2+df2-1)-(df2-1)*log(k2)-lfactorial(k2) + (-k2-df2)*log(1+Hfam2/k2))
  } else if (frailty == "cgamma") {
    llik3 <- -log(apply(cbind(df1, df2, Hfam1, Hfam2), 1, fp_fun2, k0=k0, k1=k1, k2=k2))
  }
  
  llik <- llik1 + llik2 + llik3
  
  # Ascertainment correction by design="pop+"
  cage_p <- data$currentage[ip] - agemin
  status_p <- data$status[ip]
  tvc_p <- tvc[ip]
  itvc_p <- ifelse(tvc_p < cage_p, 1, 0)
  
  if (tvctype == "PE") {
    cumhaz1_p <- cumhaz_tvc(tvc_p[itvc_p == 1], cage_p[itvc_p == 1], 
                            c(lambda1, rho1, beta_tvc), "PE")
  } else if (tvctype == "CO") {
    cumhaz1_p <- cumhaz_tvc(tvc_p[itvc_p == 1], cage_p[itvc_p == 1], 
                            c(lambda1, rho1, beta_tvc, eta, eta0), "CO")
  } else if (tvctype == "BS") {
    cumhaz1_p <- cumhaz_tvc(tvc_p[itvc_p == 1], cage_p[itvc_p == 1], 
                            c(lambda1, rho1, beta_tvc), "BS", ik = ik, degree = degree,
                            bndy = bndy, intc = intc)
  }
  
  H1_p <- (lambda1 * cage_p)^rho1 * exp(beta_gen1)
  H1_p[itvc_p == 1] <- ((lambda1 * tvc_p[itvc_p == 1])^rho1 + cumhaz1_p) * exp(beta_gen1)
  H2_p <- (lambda2 * cage_p)^rho2 * exp(beta_gen2)
  llik_p <- - (log(1 - (1 + H1_p/w1 + H2_p/w2)^(-k0) * 
                     (1 + H1_p/w1)^(-k1)*(1 + H2_p/w2)^(-k2)))
  
  # Final value of the log likelihood
  final_llik <- llik - llik_p 
  return(final_llik)
}


#################################################################### 
# 3. Standard errors for the parameters
f_theta_se <- function(fit, agemin, tvctype, frailty, data, ik = FALSE, degree = FALSE, 
                       bndy = FALSE, intc = FALSE) {
  theta <- fit$par
  I <- fit$hessian
  
  if (tvctype == "PE" | tvctype == "CO" | tvctype == "SN" | tvctype == "SIN") {
    #I <- hessian(comp_llik_sim, theta, data = data, agemin = agemin, tvctype = tvctype)
    var_cov <- solve(I)
    
    # Robust
    u <- jacobian(comp_llik_sim_vec, theta, data = data, agemin = agemin, tvctype = tvctype, frailty = frailty)
    #u <- na.omit(u)
    J <- t(u) %*% u
    robust_var <- var_cov %*% J %*% var_cov
    robust_se <- sqrt(diag(robust_var))
    se <- sqrt(diag(var_cov))
    return(list(var_cov, robust_var, robust_se, se))
  } else {
    #I <- hessian(comp_llik_sim, theta, data = data, agemin = agemin, tvctype = tvctype, 
    #             ik = ik, degree = degree, bndy = bndy, intc = intc)
    var_cov <- solve(I)
    # Robust
    u <- jacobian(comp_llik_sim_vec, theta, data = data, agemin = agemin, tvctype = tvctype, frailty = frailty,
                  ik = ik, degree = degree, bndy = bndy, intc = intc)
    #u <- na.omit(u)
    J <- t(u) %*% u
    robust_var <- var_cov %*% J %*% var_cov
    robust_se <- sqrt(diag(robust_var))
    se <- sqrt(diag(var_cov))
    
    return(list(var_cov, robust_var, robust_se, se))
  }
}

#################################################################### 
# 4. function for the coverage
sim_cover <- function(par_se, true_par, type) {
  par <- par_se[1:(length(par_se)/2)]
  se <- par_se[(length(par_se)/2+1):(length(par_se))]
  coverage <- vector()
  
  if (type == TRUE) {
    n <- length(par_se)/2
  } else {
    n <- 9
  }
  
  for (i in 1:n) {
    ci <- par[i] + c(-1, 1) * 1.96 * se[i]
    coverage[i] <- (ci[1] <= true_par[i]) * (true_par[i] <= ci[2])
  }
  return(coverage)
}



