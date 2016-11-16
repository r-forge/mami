print.lae<-function(x, ...){
cat("Coefficients:\n")
printCoefmat(x$coefficients)
cat("\n")
if(is.null(x$variable.importance)==FALSE){printCoefmat(x$variable.importance)}
cat("\n")
cat("LAE weights:\n")
print(x$sae.weights)
}
