print.mma<-function(x, ...){
cat("Coefficients:\n")
printCoefmat(round(x$coefficients, digits=3))
}
