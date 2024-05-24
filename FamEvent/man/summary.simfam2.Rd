\name{summary.simfam2}
\alias{summary.simfam2}

\title{
Summary method for \code{simfam2}
}
\description{
Provides a summary of simulated data.
}
\usage{
\method{summary}{simfam2}(object, digits = max(3, getOption("digits") - 3), \ldots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object class of \code{'simfam2'} generated from \code{simfam2} function
}
 \item{digits}{Number of significant digits to use when printing.}
 \item{...}{Further arguments passed to or from other methods.}


}
\value{
Displays a summary of simulated data and returns the following values:
\item{num.fam}{Number of families simulated.}
\item{avg.num.affected}{Average number of affected individuals per family.}
\item{avg.num.carriers}{Average number of mutation carriers per family.}
\item{avg.family.size}{Average family size.}
\item{ave.ageonset}{Average age of onset for affected individuals.}
}

\author{
Yun-Hee Choi	
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{simfam2}}
}
\examples{
set.seed(4321)

data <- simfam(N.fam = 10, design = "noasc", variation = "none",
         base.dist = "Weibull", base.parms = c(0.016, 3), vbeta = c(1, 1))

IBDmatrix <- diag(1, dim(data)[1])
data <- data[ , c(1:7, 11, 14)]

fam2 <- simfam2(inputdata = data, IBD = IBDmatrix, design = "pop", 
        variation = c("kinship","IBD"), depend = c(1, 1), 
        base.dist = "Weibull", base.parms = c(0.016, 3),
        var_names = c("gender", "mgene"), vbeta = c(1,1),
        agemin=20) 

summary(fam2)
\dontrun{
Study design:          pop 
Baseline distribution: Weibull 
Frailty distribution:  lognormal with kinship and IBD matrices 
Number of families:                     50 
Average number of affected per family:  1.24 
Average number of carriers per family:  1.3 
Average family size:                    17.02 
Average age of onset for affected:      40.08 
Sampling weights used:                  1
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ summary }
