print.mami<-function(x,string="--",...){

if(is.null(x$coefficients.s)==FALSE){x$coefficients.s <- round(x$coefficients.s,digits=6)}
if(is.null(x$coefficients.ma)==FALSE){x$coefficients.ma <- round(x$coefficients.ma,digits=6)}
if(is.null(x$coefficients.boot.s)==FALSE){x$coefficients.boot.s <- round(x$coefficients.boot.s,digits=6)}
if(is.null(x$coefficients.ma.boot)==FALSE){x$coefficients.ma.boot <- round(x$coefficients.ma.boot,digits=6)}

blank <- function(mymatrix){
kn<-rownames(mymatrix)
if(is.null(mymatrix)){mymatrix<-NULL} else {
mymatrix <- as.data.frame(mymatrix, row.names=make.names(rownames(mymatrix),unique=TRUE))
if(is.null(mymatrix)==TRUE){mymatrix<-mymatrix}else{ 
for(i in 1:dim(mymatrix)[1]){
if(all(mymatrix[i,]==1) | all(mymatrix[i,]==0) | all(mymatrix[i,] %in% c(NA,Inf,NaN,0,1))){mymatrix[i,]<-c(rep(string,length(mymatrix[i,])))}
}}
rownames(mymatrix)<-kn
}
return(mymatrix)
}
x$coefficients.s <- blank(x$coefficients.s)
x$coefficients.boot.s <- blank(x$coefficients.boot.s)
x$coefficients.ma <- blank(x$coefficients.ma)
x$coefficients.ma.boot <- blank(x$coefficients.ma.boot)
if(x$setup[[19]]==TRUE){
x$coefficients.ma <- x$coefficients.ma[,-2]
x$coefficients.ma <- x$coefficients.ma[,-5]
x$coefficients.s  <- x$coefficients.s[,-2]
x$coefficients.boot.s  <- x$coefficients.boot.s[,-5]
}

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
if(x$setup[[6]]!="BIC"& x$setup[[6]]!="BIC+"){cat("Variable importance (based on model averaging weights):\n\n")}else{
cat("Posterior effect probabilities:\n\n")
}
print(round(x$variable.importance,digits=2))}

if(is.null(x$setup[[18]])==FALSE){
cat("____________________\n")
cat(paste("The following variables have been excluded by screening:",paste(x$setup[[18]],collapse=" "),"\n"))
}

cat("\n")

}