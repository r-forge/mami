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
\method{summary}{mami}(object,string="x", ...)
}
%- maybe also 'usage' for other objects documented here.
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
\code{\link{mami}} for model selection and model averaging after multiple imputation. 
}
\examples{
library(Amelia)
data(freetrade)
freetrade$pop <- log(freetrade$pop) # in line with original publication
freetrade_imp <- amelia(freetrade, p2s=0, ts = "year", cs = "country", noms="signed",
                        polytime = 2, intercs = TRUE, empri = 2)
m1 <- mami(freetrade_imp, method="criterion.average", ycol="tariff",add.factor=c("country"))
summary(m1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Model Averaging}
\keyword{Model Selection}
\keyword{Multiple Imputation}