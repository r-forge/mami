\name{summary.lae}
\alias{summary.lae}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarizes results from Lasso Averaging Estimation
}
\description{
Summarizes the estimates and standard errors of Lasso (averaging) estimation and OLS estimation and reports variable importance.
}
\usage{
\method{summary}{lae}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
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
summary(l1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Shrinkage}
\keyword{Model Averaging}
