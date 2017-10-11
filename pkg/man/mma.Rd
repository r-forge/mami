\name{mma}
\alias{mma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MMA: Mallow's Model Averaging 
}
\description{
Performs model averaging on a set of nested candidate models with the weights vector chosen such that a specific Mallow's criterion is minimized.
}
\usage{
mma(X, formula=NULL, ycol=1, variance=c("BA","boot"), bsa=200)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A dataframe or matrix of data.
}
  \item{formula}{
Formula of the full model.
}
  \item{ycol}{
Either a character vector or integer specifying the column with the outcome variable.
}
\item{variance}{A character vector specifying whether the variance is estimated due to the formula of Burnham
 and Anderson (\code{"BA"}) or based on bootstrapping (\code{"boot"}). Both options are pragmatic solutions.
 }
\item{bsa}{A positive integer specifying the number of bootstrap samples used if \code{variance} \code{=} \code{"boot"}.}
}
\details{
Mallow's Model Averaging (MMA) considers all nested submodels of the full model as candidate models, i.e if there are 7 variables there are 7 candidate models. 
The weight vector used to combine the respective estimates is chosen such that a certain Mallow's type criterion is minimized. The final weighted estimate
targets to minimize the mean squared prediciton error and is optimal in some sense, see Theorem 1 and Lemma 3 in Hansen, B. (2007, \emph{Least Squares Model Averaging}, Econometrica, 75:1175-1189).

Note however that the results of MMA depend on the ordering of the regresssors. 
}
\value{
Returns an object of \code{class} `mma': 
\item{coefficients}{A matrix of estimated coefficients and standard errors (and bootstrap standard errors if \code{variance} \code{=} \code{"boot"}).}
\item{averaging.weights}{A matrix containing the weights for each models used in the model averaging procedure.}
}
\references{
Hansen, B. (2007), \emph{Least Squares Model Averaging}, Econometrica, 75:1175-1189
}
\author{
Michael Schomaker
}
\examples{
data(Prostate)
mma(Prostate,formula=lpsa~.,ycol="lpsa")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Model Averaging}

