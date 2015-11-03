\name{plot.lae}
\alias{plot.lae}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualize estimates from Lasso Averaging Estimation
}
\description{
The function plots (i) the weights for each candidate model relating to different candidate complexity parameters and (ii) the Lasso
averaging, Lasso, OLS estimation, and variable importance results.
}
\usage{
\method{plot}{lae}(x, xaxis = c("index", "realnumbers"), legend.place = "topright", 
         display.importance = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class \code{"lae"}.
}
  \item{xaxis}{
A string of either \code{"index"} or \code{"realnumbers"} to specify whether the x-axis should show the real values of the candidate complexity parameters or not.
}
  \item{legend.place}{
A string defining where the legend should be placed.
}
  \item{display.importance}{
A logical value defining whether the variable importance measure should be plotted or not.
}
  \item{\dots}{
Further graphical arguments to be passed.
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
l1 <- lae(Prostate, ycol=9, kfold=10, my.formula=lpsa~.)
plot(l1,xaxis="realnumbers")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Shrinkage}
\keyword{Model Averaging}% __ONLY ONE__ keyword per line
