\name{lae}
\alias{lae}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Lasso Averaging Estimation
}
\description{
Lasso (least absolute shrinkage and selection operator) estimation is performed and evaluated for different tuning parameter choices. To address
tuning parameter selection uncertainty a weighted average of these estimators is calculated. The weight vector is chosen such that a k-fold
cross validation criterion is minimized. 
}
\usage{
lae(X, ycol = 1, B.var = 50, nolambda = 100, kfold = 5, my.formula = NULL,
 standardize = TRUE, calc.variance = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A dataframe or matrix containing the data to be analyzed.
}
  \item{ycol}{
An integer specifying the column of the outcome variable.
}
  \item{B.var}{
An integer specifying the number of bootstrap replications to be used to estimate the standard error of the Lasso estimator.
}
  \item{nolambda}{
An integer specifying the number of candidate complexity parameters to consider.
}
  \item{kfold}{
An integer specifying the \code{kfold} cross validation criterion to (i) use for tuning parameter selection (ii) be minimized for Lasso averaging estimation.
}
  \item{my.formula}{
A formula specifying the full model.
}
  \item{standardize}{
A logical value speciying whether the covariate data should be standardized.
}
  \item{calc.variance}{
A logical value specifying whether the standard error of the estimates should be estimated at all (by means of bootstrapping).
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns an object of \code{class} `lae': 
\item{coefficients}{A matrix of coefficients and standard errors for Lasso averaging, Lasso selection, and OLS estimation.}
\item{variable.importance}{A matrix containing the relative importance of each variable based on model averaging weights.}
\item{sae.weights}{A vector containing the weights used for Lasso averaging.}
\item{sel.weights}{A vector indicating the complexity parameter that was chosen for Lasso estimation based on k-fold cross validation.}
\item{complexity.parameter}{A vector of the actual complexity parameter values used as candidate values for Lasso Averaging Estimation.}
}
\references{
Schomaker, M. (2012) \emph{Shrinkage Averaging Estimation}, Statistical Papers, 53:1015-1034
}
\author{
Michael Schomaker
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{plot.lae}} for visualizing the estimation results.
}
\examples{
library(lasso2)
data(Prostate)
lae(Prostate, ycol=9, kfold=10, my.formula=lpsa~.)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Model Averaging}
\keyword{Shrinkage}
