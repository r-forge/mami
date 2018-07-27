\name{jma}
\alias{jma}
\title{
JMA: Jackknife Model Averaging 
}
\description{
Performs model averaging on a set of (linear) candidate models with the weight vector chosen 
such that the leave-one-out cross validation error is minimized. 
}
\usage{
jma(y,x,method=c("JMA","MMA"),subset=c("nested","all"),factor.variables=NULL,pd=T,
                variance=c("none","boot"), bsa=100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
The response vector.
}
  \item{x}{
A matrix or dataframe containing the covariates.
}
  \item{method}{
A character vector specifiying whether jackknife model averaging (JMA) or
Mallow's model averaging (see \code{\link{mma}}) should be performed.
}
  \item{subset}{
A character vector specifiying whether all 2^p candidate models should be considered or the nested subset of p models. 
}
  \item{factor.variables}{
A (vector of) string(s) specifying which variables should be treated as factors, i.e. recoded into dummy variables. Factor variables will automatically be recoded if not specified here.
}
  \item{pd}{
A logical value specifying whether messages should be printed or not.
}
\item{variance}{A character vector specifying whether no variance should be estimated or based on bootstrapping (\code{"boot"}). 
 }
\item{bsa}{A positive integer specifying the number of bootstrap samples used if \code{variance} \code{=} \code{"boot"}.}
}
\details{
This function utilizes Jackknife model averaging as described in Hansen and Racine (2012), see reference below.

If \code{subset} \code{=} \code{"all"}, then 2^p candidate models are being evaluated. This means p can't be too large (say<20) and 
occasionally the residual matrix used in the quadtratic programming problem may not be positive definite. In the latter case, this matrix
is altered by \code{jma} such that it is positive definite, but results should be interpreted with care.
}
\value{
Returns an object of \code{class} `jma': 
\item{betahat}{estimates coefficients}
\item{se}{standarde error}
\item{lci}{lower confidence limit}
\item{uci}{upper confidence limit}
\item{weight}{JMA weight vector}
\item{yhat}{fitted values}
\item{ehat}{fitted residuals}
\item{y}{outcome variable}
\item{x}{matrix of covariates}
}
\references{
Hansen, B. and Racine, J. (2012), \emph{Jackknife Model Averaging}, Journal of Econometrics, 167:38-46
}
\author{
Michael Schomaker (based on the file of Bruce Hansen at \url{https://www.ssc.wisc.edu/~bhansen/progs/joe_12.html})
}
\examples{
data(Prostate)
jma(y=Prostate$lpsa,x=Prostate[,-9])
jma(y=Prostate$lpsa,x=Prostate[,-9], variance="boot", bsa=100)

# can also perform Mallow's Model Averaging as in mma()
jma(y=Prostate$lpsa,x=Prostate[,-9],method="MMA")
mma(Prostate, formula=lpsa~.,ycol="lpsa")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Model Averaging}

