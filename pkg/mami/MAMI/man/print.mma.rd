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
\code{\link{mma}} for Mallow's Model Averaging. 
}
\examples{
library(lasso2)
data(Prostate)
m1 <- mma(Prostate,modelformula=lpsa~.,ycol=9)
m1
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Model Averaging}

