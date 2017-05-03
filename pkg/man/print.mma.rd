\name{print.mma}
\alias{print.mma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Printing results of Mallow's Model Averaging
}
\description{
Prints estimates of Mallow's Model Averaging.
}
\usage{
\method{print}{mma}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class \code{"mma"}.
}
  \item{\dots}{
Further arguments to be passed.
}
}
\author{
Michael Schomaker
}
\seealso{
\code{\link{mma}} for Mallow's Model Averaging. 
}
\examples{
data(Prostate)
m1 <- mma(Prostate,lpsa~.,ycol="lpsa")
m1
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Model Averaging}

