\name{lae}
\alias{lae}
\title{
Lasso Averaging Estimation
}
\description{
Lasso (least absolute shrinkage and selection operator) estimation is performed and evaluated for different tuning parameter choices. To address
tuning parameter selection uncertainty a weighted average of these estimators is calculated. The weight vector is chosen such that a k-fold
cross validation criterion is minimized. 
}
\usage{
lae(X, ycol = 1, kfold = 10, B.var = 100, calc.variance = FALSE,
 factor.variables = NULL, glm.family = "gaussian", tries = 10,
 standardize = TRUE, random = FALSE, pd = TRUE, ...)
}
\arguments{
  \item{X}{
A dataframe or matrix containing the data to be analyzed.
}
  \item{ycol}{
An integer or string specifying the column of the outcome variable.  The outcome for \code{glm.family="cox"} should be a two-column matrix with columns named 'time' and 'status'. The latter is a binary variable, with '1' indicating event, and '0' indicating right censored.
}
  \item{kfold}{
An integer specifying the \code{kfold} cross validation criterion to (i) use for tuning parameter selection (ii) be minimized for Lasso averaging estimation.
}
  \item{B.var}{
An integer specifying the number of bootstrap replications to be used to estimate the standard error of the Lasso estimator. 
}
  \item{calc.variance}{
A logical value specifying whether the standard error of the estimates should be estimated at all (by means of bootstrapping). See also details below.
}
  \item{factor.variables}{
A (vector of) string(s) specifying which variables should be treated as factors, i.e. recoded into dummy variables. Factor variables will automatically be recoded if not specified here.
}
  \item{glm.family}{
A character vector specifying one of the following families: \code{"gaussian"}, \code{"binomial"}, \code{"poisson"}.
}
  \item{tries}{
An integer for the number of tries in case \code{lae} fails; can be relevant if sub-datasets used during cross-validation lead to failure of lasso averaging. In this case, the sub-datasets are randomly re-defined.
}
  \item{standardize}{
A logical value speciying whether the covariate data should be standardized.
}
  \item{random}{
A logical value specifying whether creation of datasets for cross validation should be random or not.
}
  \item{pd}{
A logical value specifying whether messages should be printed or not.
}
  \item{\dots}{
Other arguments to be passed, i.e. to \code{cv.glmnet}. For example, pass \code{alpha=0} for the Ridge estimator and \code{alpha} smaller than 1 for Elastic Net.
}
}
\details{
Note that the candidate tuning parameters are selected automatically (by \code{cv.glmnet}). The bootstrap standard error for LASSO does not assume a fixed
tuning parameter, i.e. tuning parameter selection is done seperately in each bootstrap sample. Lasso averaging works on standard errors related to each tuning parameter,
but the variance between the different weighted estimates is taken into account. The importance measure based on the 
averaging weights could be interpreted as the importance of variables with respect to their predictive ability.
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
\seealso{
\code{\link{plot.lae}} for visualizing the estimation results.
}
\examples{
data(Prostate)
lae(Prostate, ycol="lpsa")
lae(Prostate, ycol="lpsa", factor.variables="gleason",calc.variance=TRUE, kfold=5)
}
\keyword{Model Averaging}
\keyword{Shrinkage}