\name{Prostate}
\alias{Prostate}
\docType{data}
\title{
Prostate data
}
\description{
These data come from a study that examined the correlation between the level of prostate specific antigen and a number of clinical
measures in men who were about to receive a radical prostatectomy. The data set has been 
copied from library \code{lasso2}. 
}
\usage{data(Prostate)}
\format{
  It is data frame with 97 rows and 9 columns.
  \describe{
    \item{\code{lcavol}}{log(cancer volume)}
    \item{\code{lweight}}{log(prostate weight)}
    \item{\code{age}}{age}
    \item{\code{lbph}}{log(benign prostatic hyperplasia amount)}
    \item{\code{svi}}{seminal vesicle invasion}
    \item{\code{lcp}}{log(capsular penetration)}
    \item{\code{gleason}}{Gleason score}
    \item{\code{pgg45}}{percentage Gleason scores 4 or 5}
    \item{\code{lpsa}}{log(prostate specific antigen)}
  }
}

\references{
Stamey, T. et al. (1989) \emph{Prostate specific antigen in the diagnosis and treatment of adenocarcinoma of the prostate: II. radical prostatectomy treated patients}, Journal of Urology, 141:1076-1083
}
\examples{
str(Prostate)
}
\keyword{Datasets}
