summary.mami<-function(object,string="x",...){

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
object$coefficients.s <- blank(object$coefficients.s)
object$coefficients.boot.s <- blank(object$coefficients.boot.s)
object$coefficients.ma <- blank(object$coefficients.ma)
object$coefficients.ma.boot <- blank(object$coefficients.ma.boot)
factor.cleaning<-function(mystring){
cf<-NULL
for(i in 1:length(mystring)){
if(substr(mystring[i],1,10)=="as.factor("){cf<-c(cf, paste(substr(mystring[i],1,nchar(strsplit(mystring[i],")")[[1]][1])),")",sep=""))}else{
cf<-c(cf,mystring[i])}}
cf
}


cat("\n")
MA<-NULL
if(is.null(object$coefficients.ma)==FALSE){
MA<-as.data.frame(object$coefficients.ma[,c(1,3,4)])
MA.title <- "Estimates for model averaging: \n\n"}
if(is.null(object$coefficients.ma.boot)==FALSE){
MA<-cbind(MA,object$coefficients.ma.boot[,2:3])
MA.title <- paste("Estimates for model averaging (based on",dim(object$boot.results[[2]])[1],"bootstrap samples):\n\n",sep=" ")
colnames(MA)[2:5]<-c("LCI","UCI","Boot LCI","Boot UCI")
}
if(is.null(object$variable.importance)==FALSE){
MA$VI<-NA
subset1<- match(names(object$variable.importance),factor.cleaning(rownames(MA)))
MA$VI[subset1]<-object$variable.importance
MA[is.na(MA)]<-"--"
}

if(is.null(object$coefficients.ma)==FALSE){
cat(MA.title)
print(MA)
cat("\n")}
  
MS<-NULL
if(is.null(object$coefficients.s)==FALSE){
MS<-as.data.frame(object$coefficients.s[,c(1,3,4)])
MS.title <- "Estimates for model selection: \n\n"}
if(is.null(object$coefficients.boot.s)==FALSE){
MS<-cbind(MS,object$coefficients.boot.s[,2:3])
MS.title <- paste("Estimates for model selection (based on",dim(object$boot.results[[1]])[1],"bootstrap samples):\n\n",sep=" ")
colnames(MS)[2:5]<-c("LCI","UCI","Boot LCI","Boot UCI")
}
if(is.null(object$variable.importance)==FALSE){
MS$VI<-NA
subset2 <- match(names(object$variable.importance),factor.cleaning(rownames(MS)))
MS$VI[subset2]<-object$variable.importance
MS[is.na(MS)]<-"--"
}

if(is.null(object$coefficients.s)==FALSE){
cat(MS.title)
print(MS)
cat("\n")}

}