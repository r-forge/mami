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
\author{
Michael Schomaker
}
\seealso{
\code{\link{lae}} for Lasso Averaging Estimation.
}
\examples{
data(Prostate)
l1 <- lae(Prostate, ycol="lpsa", kfold=10)
plot(l1, xaxis="realnumbers")
}
\keyword{Shrinkage}
\keyword{Model Averaging}
