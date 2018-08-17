print.jma<-function(x, ...){
cat("Coefficients:\n")
printCoefmat(round(cbind(x$betahat,t(x$se),t(x$lci),t(x$uci)), digits=5))
}
