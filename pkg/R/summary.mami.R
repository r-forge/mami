summary.mami<-function(object,string="x",...){

blank <- function(mymatrix){
if(is.null(mymatrix)){mymatrix<-NULL} else {
mymatrix <- as.data.frame(mymatrix, row.names=make.names(rownames(mymatrix),unique=TRUE))
if(is.null(mymatrix)==TRUE | any(is.na(mymatrix))==TRUE){mymatrix<-mymatrix}else{ 
for(i in 1:dim(mymatrix)[1]){
if(all(mymatrix[i,]==1) | all(mymatrix[i,]==0)){mymatrix[i,]<-c(rep(string,length(mymatrix[i,])))}
}}
}
return(mymatrix)
}
object$coefficients.s <- blank(object$coefficients.s)
object$coefficients.boot.s <- blank(object$coefficients.boot.s)
object$coefficients.ma <- blank(object$coefficients.ma)
object$coefficients.ma.boot <- blank(object$coefficients.ma.boot)

cat("\n")

if(is.null(object$coefficients.ma)==FALSE){
cat("Estimates for model averaging:\n\n")
print(object$coefficients.ma)
cat("\n")}

if(is.null(object$coefficients.ma.boot)==FALSE){
cat(paste("Bootstrap estimates for model averaging (based on",dim(object$boot.results[[2]])[1],"bootstrap samples):\n\n",sep=" "))
print(object$coefficients.ma.boot)
cat("\n")}

if(is.null(object$coefficients.s)==FALSE){
cat("Estimates for model selection:\n\n")
print(object$coefficients.s)
cat("\n")}

if(is.null(object$coefficients.boot.s)==FALSE){
cat(paste("Bootstrap estimates for model selection (based on",dim(object$boot.results[[1]])[1],"bootstrap samples):\n\n",sep=" "))
print(object$coefficients.boot.s)
cat("\n")}

if(is.null(object$variable.importance)==FALSE){
cat("Variable importance (based on model averaging weights):\n\n")
print(object$variable.importance)}

cat("\n")


}