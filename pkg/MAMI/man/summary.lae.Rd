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
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Michael Schomaker
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

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
