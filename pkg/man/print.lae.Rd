\name{print.lae}
\alias{print.lae}
\title{
Printing results of Lasso Averaging Estimation
}
\description{
Prints stimates and standard errors of Lasso (averaging) estimation and OLS estimation and reports variable importance and weights.
}
\usage{
\method{print}{lae}(x, ...)
}
\arguments{
  \item{x}{
An object of class \code{"lae"}.
}
  \item{\dots}{
Further arguments to be passed.
}
}
\author{
Michael Schomaker
}
\seealso{
\code{\link{lae}} for Lasso Averaging Estimation.
}
\examples{
library(lasso2)
data(Prostate)
l1<-lae(Prostate, ycol=9, kfold=10, my.formula=lpsa~.)
print(l1)
}
\keyword{Model Averaging}
\keyword{Shrinkage}
