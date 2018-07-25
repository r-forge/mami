\name{listSLWrappers}
\alias{listSLWrappers}
\title{
Explains model averaging wrappers for super learning
}
\description{
Explains model averaging wrappers which can be used in conjunction with package \pkg{SuperLearner} (and \pkg{tmle} and \pkg{ltmle})
}
\usage{
listSLWrappers()
}
\author{
Michael Schomaker
}
\examples{
\dontrun{
listSLWrappers()
library(SuperLearner) # needs to be installed
SL.library <- c("SL.glm","SL.stepAIC", "SL.mean", "SL.step.interaction", "SL.mma.int", "SL.jma")
SL.library2 <- c("SL.glm","SL.stepAIC", "SL.mean", "SL.step.interaction", "SL.lae2")

data(Prostate)
P1 <- SuperLearner(Y=Prostate[,9], X=Prostate[,-9], SL.library = SL.library, verbose=T)
P2 <- SuperLearner(Y=Prostate[,5], X=Prostate[,-5], family="binomial", SL.library = SL.library2, verbose=T)
P2$coef
}
}
\keyword{super learning}
\keyword{wrapper}

