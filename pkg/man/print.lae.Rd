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
data(Prostate)
l1<-lae(Prostate, ycol="lpsa", kfold=10)
print(l1)
}
\keyword{Model Averaging}
\keyword{Shrinkage}
