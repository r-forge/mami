print.mami<-function(x,string="x",...){

blank <- function(mymatrix){
kn<-rownames(mymatrix)
if(is.null(mymatrix)){mymatrix<-NULL} else {
mymatrix <- as.data.frame(mymatrix, row.names=make.names(rownames(mymatrix),unique=TRUE))
if(is.null(mymatrix)==TRUE | any(is.na(mymatrix))==TRUE){mymatrix<-mymatrix}else{ 
for(i in 1:dim(mymatrix)[1]){
if(all(mymatrix[i,]==1) | all(mymatrix[i,]==0)){mymatrix[i,]<-c(rep(string,length(mymatrix[i,])))}
}}
rownames(mymatrix)<-kn
}
return(mymatrix)
}
x$coefficients.s <- blank(x$coefficients.s)
x$coefficients.boot.s <- blank(x$coefficients.boot.s)
x$coefficients.ma <- blank(x$coefficients.ma)
x$coefficients.ma.boot <- blank(x$coefficients.ma.boot)




cat("\n")

if(is.null(x$coefficients.ma)==FALSE){
cat("Estimates for model averaging:\n\n")
print(x$coefficients.ma)
cat("\n")}

if(is.null(x$coefficients.ma.boot)==FALSE){
cat(paste("Bootstrap estimates for model averaging (based on",dim(x$boot.results[[2]])[1],"bootstrap samples):\n\n",sep=" "))
print(x$coefficients.ma.boot)
cat("\n")}

if(is.null(x$coefficients.s)==FALSE){
cat("Estimates for model selection:\n\n")
print(x$coefficients.s)
cat("\n")}

if(is.null(x$coefficients.boot.s)==FALSE){
cat(paste("Bootstrap estimates for model selection (based on",dim(x$boot.results[[1]])[1],"bootstrap samples):\n\n",sep=" "))
print(x$coefficients.boot.s)
cat("\n")}

if(is.null(x$variable.importance)==FALSE){
cat("Variable importance (based on model averaging weights):\n\n")
print(x$variable.importance)}

cat("\n")

}