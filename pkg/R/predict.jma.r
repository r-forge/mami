predict.jma <- function(object,newdata=NULL,...){
    if(is.null(newdata)){newdata<-data.frame(object$x[,-1])}
    newdata <- newdata[,colnames(newdata)%in%colnames(object$x)]
    pred <- model.matrix(~.,data=newdata)%*%matrix(object$betahat,ncol=1)
    return(pred)
}