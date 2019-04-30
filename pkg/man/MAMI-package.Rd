\name{mami-package}
\alias{mami-package}
\alias{MAMI}
\docType{package}
\title{
Model averaging (and model selection) after multiple imputation
}
\description{
Performs model selection/averaging on multiply imputed data sets (possibly created with \code{Amelia II} or \code{mice}) and combines the resulting estimates. The package also provides access to less frequently used 
model averaging techniques and offers integrated bootstrap estimation. The main function is \code{"mami"}.
}

\details{
\tabular{ll}{
Package: \tab MAMI\cr
Type: \tab Package\cr
Version: \tab 0.9.12.2\cr
Date: \tab 2019-04-30\cr
License: \tab GPL-2\cr
Depends: \tab R (>= 3.5.0)\cr
Imports: \tab quadprog, boot, MASS, survival, MuMIn, zoo, lme4, Amelia, BMA\cr
\tab glmnet, parallel, snow, doParallel, foreach\cr
Suggests: mice \cr
}

}
\author{
Michael Schomaker

Maintainer: Michael Schomaker <michael.schomaker@uct.ac.za>
}
\references{
Schomaker, M., Heumann, C. (2014) \emph{Model Selection and Model Averaging after Multiple Imputation},
Computational Statistics & Data Analysis, 71:758-770
}
\keyword{Model Averaging}

