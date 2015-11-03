summary.mami<-function(object,...){

cat("\n")

if(is.null(object$coefficients)==FALSE){
cat("Estimates for model averaging:\n\n")
print(object$coefficients)
cat("\n")}

if(is.null(object$coefficients.boot)==FALSE){
cat(paste("Bootstrap estimates for model averaging (based on",dim(object$boot.results[[2]])[1],"samples):\n\n",sep=" "))
print(object$coefficients.boot)
cat("\n")}

if(is.null(object$coefficients.s)==FALSE){
cat("Estimates for model selection:\n\n")
print(object$coefficients.s)
cat("\n")}

if(is.null(object$coefficients.boot.s)==FALSE){
cat(paste("Bootstrap estimates for model selection (based on",dim(object$boot.results[[1]])[1],"samples):\n\n",sep=" "))
print(object$coefficients.boot.s)
cat("\n")}

if(is.null(object$variable.importance)==FALSE){
cat("Variable importance (based on model averaging weights):\n\n")
print(object$variable.importance)}

cat("\n")

}