familyDesign2 <-
  function(depend = 1, affectnum = 0,
           base.dist = "Weibull", base.parms = c(0.016, 3),
           var_names = NULL, vbeta= c(-1.126, 2.55, 1.6), 
           agemin=20, data, corr_type = "frailty", IBD = NULL)
  {
    data.new<-numeric()
    n <- length(unique(data$famID))
    cumind<-0
    i<- 1
    j<- 0
    
    if(!is.null(IBD) & !is.positive.definite(IBD)) stop("IBD matrix must be a positive definite matrix")

    while (i <= n) {
      pedigreeFam <- data[data$famID == i,]
      j <- j + 1
      dat <- try(familyStructure2(i,cumind=cumind, depend = depend, 
                                base.dist=base.dist, base.parms = base.parms, 
                                data = pedigreeFam, corr_type = corr_type, 
                                IBD = IBD[data$famID == i, data$famID == i],
                                var_names = var_names, vbeta=vbeta, agemin=agemin))
      
      if(affectnum == 1) until <- ifelse(sum(dat$status[dat$proband==1]) == 1, TRUE, FALSE)
      else if(affectnum == 3) until <- ifelse(sum(dat$status[dat$generation %in% c(1,2)]) == 3, TRUE, FALSE)
      else until <- TRUE

      if(!is.null(dim(dat))){
        if(nrow(dat)>0 ){
          if(until){
            data.new<-rbind(data.new, dat)
            cumind<-cumind+nrow(dat)
            i<-i+1
          }    
        }
      }
      
   } # close while
    
    data.new
  }
