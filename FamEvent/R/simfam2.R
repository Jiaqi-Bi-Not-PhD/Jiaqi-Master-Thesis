simfam2 <- function(inputdata = NULL, IBD = NULL, design="pop", variation="none", 
           depend=NULL, base.dist="Weibull", base.parms=c(0.016,3), 
           var_names = c("gender", "mgene"), vbeta=c(1,1), 
           agemin=20, hr=NULL)
{
 # kinship and IBD only    
    
  if(!any(is.element(variation, c("kinship", "IBD")))) stop("Unrecognized variation; variation should be either kinship or IBD")
  if(length(variation)==2 & length(depend)!=2) stop("depend should be a vector of two elements.")
  if(length(variation)==1 & length(depend)!=1) stop("depend should be a single value.")
  

  if(agemin > 70) warning("agemin might be too high.")
  
  if(!is.element(base.dist, c("Weibull", "loglogistic", "Gompertz", "lognormal", "logBurr"))) 
      stop("Unrecognized base.dist; base.dist should be one of \"Weibull\", \"loglogistic\", \"Gompertz\", \"lognormal\", or \"logBurr\". ")
    
    if(is.null(depend)) stop("depend should be specified.")
    else if(any(depend<0)) stop("depend should be > 0")
    
    if(!is.element(design, c("pop","pop+","cli","cli+","noasc"))) stop("Unrecognized design; should be one of pop, pop+, cli, cli+, twostage, or noasc")  
    
    if(length(vbeta) != length(var_names)) stop(paste("vbeta should be a vector of length", length(var_names)))
    if(base.dist=="logBurr" & length(base.parms)!=3) stop("base.parms should be a vector of length 3")
    else if(base.dist!="logBurr" & length(base.parms)!=2) stop("base.parms should be a vector of length 2")
    
    if(length(variation)==1){
      if(!is.element(variation, c("kinship", "IBD"))) stop("variation should be one of none, frailty, secondgene, kinship, IBD, or both kinship and IBD")
    }else{
      if(!all(is.element(variation, c("kinship", "IBD")))) stop("variation should be either kinship or IBD, or both kinship and IBD")
    }

    if(design=="pop") {affectnum=1; m.carrier=0}
    if(design=="pop+"){affectnum=1; m.carrier=1}
    if(design=="cli") {affectnum=3; m.carrier=0}
    if(design=="cli+"){affectnum=3; m.carrier=1}
    if(design=="twostage"){
      affectnum = 2
      m.carrier = 0 #proband is not necessary to be a carrier.
      if(hr==0) stop("Please specify the sampling rate of high risk families (0<hr<1)")
    }
    else if(design == "noasc"){affectnum=0; m.carrier=0}
    

    dat <- data.frame(familyDesign2(depend=depend, affectnum=affectnum, 
                            base.dist=base.dist, base.parms=base.parms,
                            var_names=var_names, vbeta=vbeta, agemin=agemin, 
                            data = inputdata, corr_type = variation, IBD = IBD))
    
    class(dat) <- c("simfam","data.frame")
    attr(dat, "design") <- design
    attr(dat, "variation") <- variation
    attr(dat, "frailty.dist") <- "lognormal"
    attr(dat, "base.dist") <- base.dist
    attr(dat, "base.parms") <- base.parms
    attr(dat, "vbeta") <- vbeta
    attr(dat, "agemin") <- agemin
    return(dat)
  }
