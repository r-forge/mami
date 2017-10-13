\name{summary.mami}
\alias{summary.mami}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarizes results of model selection and model averaging after multiple imputation
}
\description{
Summarizes estimates, standard errors and confidence intervals for model selection/averaging estimators
after multiple imputation and reports variable importance if applicable.
}
\usage{
\method{summary}{mami}(object,string="--", ...)
}
\arguments{
  \item{object}{
An object of class \code{"mami"}.
}
  \item{string}{
A string specifying what is printed in the summary table for non-selected variables.
}
  \item{\dots}{
Further arguments to be passed.
}
}
\author{
Michael Schomaker
}
\seealso{
\code{\link{mami}} for model selection and model averaging after multiple imputation. 
}
\examples{
library(Amelia)
data(freetrade)
freetrade$pop <- log(freetrade$pop) # in line with original publication
freetrade_imp <- amelia(freetrade, p2s=0, ts = "year", cs = "country", noms="signed",
                        polytime = 2, intercs = TRUE, empri = 2)
m1 <- mami(freetrade_imp, method="MA.criterion", outcome="tariff",add.factor=c("country"))
summary(m1)
}
\keyword{Model Averaging}
\keyword{Model Selection}
\keyword{Multiple Imputation}