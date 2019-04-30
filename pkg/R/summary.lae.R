summary.lae <- function(object,...){

cat("Summary of LASSO averaging, LASSO and OLS estimation:\n")
cat("\n")
printCoefmat(object$coefficients)
cat("\n")
if(is.null(object$variable.importance)==FALSE){printCoefmat(object$variable.importance)}

}