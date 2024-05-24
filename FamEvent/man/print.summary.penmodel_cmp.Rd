\name{print.summary.penmodel_cmp}
\alias{print.summary.penmodel_cmp}

\title{
  Print method for \code{summary.penmodel_cmp} of a fitted competing risks penetrance model.
}
\description{
  Prints a short summary of parameter and penetrance estimates of a \code{'summary.penmodel_cmp'} object.
}
\usage{
\method{print}{summary.penmodel_cmp}(x, digits = max(3, getOption("digits") - 3),
  signif.stars=TRUE, \ldots)
}
\arguments{
  \item{x}{
    An object class of \code{'summary.penmodel_cmp'}, a result of a call to \code{\link{summary.penmodel_cmp}}.}
  \item{digits}{Number of significant digits to use when printing.}
  \item{signif.stars}{Logical; if \code{TRUE}, provides stars to hightlight significant p-values. Default is \code{TRUE}.}
  \item{...}{Further arguments passed to or from other methods.}
}

\value{
Prints a summary of parameter estimates, their standard errors, \eqn{t}-statistics and corresponding two-sided \eqn{p}-values and additionally indicates significance stars if \code{signif.stars} is \code{TRUE}.

Also prints penetrance estimates for each event by age 70 specific to gender and mutation-status subgroups along with their standard errors and 95\% confidence intervals.

Returns an object of class \code{'summary.penmodel_cmp'}.

}
\author{
  Yun-Hee Choi	
}

\seealso{
    \code{\link{penmodel_cmp}, \link{print.penmodel_cmp}, \link{summary.penmodel_cmp}}
}

\keyword{ print }
