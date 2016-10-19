\name{plot.mami}
\alias{plot.mami}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot bootstrap results when utilizing model selection/averaging after multiple imputation
}
\description{
Plots the bootstrap distribution of model selection/averaging estimators after multiple imputation. By default, the areas of the boostrap distribution that
are left and right of zero are shaded in different colors to better identify the importance of each variable.
}
\usage{
\method{plot}{mami}(x, shade.areas = TRUE, plots.p.page = c("1", "4", "9", "16"), 
color = c("lightgrey", "darkgrey"), adj.bw = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class \code{"mami"}.
}
  \item{shade.areas}{
A logical value indicating whether areas left and right of zero should be shaded or not.
}
  \item{plots.p.page}{
A string of either \code{"1"}, \code{"4"}, \code{"9"}, or \code{"16"} indicating the number of plots per page.
}
  \item{color}{
A string vector of length 2 indicating the colors used for the shaded areas.
}
  \item{adj.bw}{
An integer to adjust the bandwidth.
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(Amelia)
data(freetrade)
freetrade$pop <- log(freetrade$pop) # in line with original publication
freetrade_imp <- amelia(freetrade, ts = "year", cs = "country", noms="signed", polytime = 2,
                        intercs = TRUE, empri = 2)
m1 <- mami(freetrade_imp, method="criterion.selection", ycol="tariff",add.factor=c("country"),
           inference="+bootstrapping", B=25,X.org=freetrade,print.time=TRUE)
m1 # be patient with bootstrapping, increase B>=200 for better results
plot(m1) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Model Averaging}
\keyword{Model Selection}
\keyword{Multiple Imputation}
