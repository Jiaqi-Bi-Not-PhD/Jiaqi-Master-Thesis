\name{print.penmodel_cmp}
\alias{print.penmodel_cmp}

\title{
Print method for \code{penmodel_cmp}.
}
\description{
 Prints a summary of parameter estimates of a fitted competing risk penetrance model.
}
\usage{
\method{print}{penmodel_cmp}(x, digits = max(3, getOption("digits") - 3), \ldots)
}
\arguments{
  \item{x}{
An object class of \code{'penmodel_cmp'}, a fitted model by \code{\link{penmodel_cmp}} function.
}
 \item{digits}{Number of significant digits to use when printing.}
 \item{...}{Further arguments passed to or from other methods.}
}
\value{
Prints a short summary of the model and model fit.

Returns an object of class \code{'penmodel_cmp'}.
}

\author{
Yun-Hee Choi	
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{penmodel_cmp}, \link{summary.penmodel_cmp}, \link{print.summary.penmodel_cmp}, \link{plot.penmodel_cmp}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ print }
