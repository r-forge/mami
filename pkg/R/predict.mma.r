predict.mma <- function(object,newdata=NULL,...){
    if(is.null(newdata)){newdata<-object$setup[[2]]}
    newdata <- newdata[,colnames(newdata)%in%colnames(object$setup[[2]])]
    mf <- paste("~",gsub(".*~","~",object$setup[[1]])[3])
    pred <- model.matrix(as.formula(mf),data=newdata)%*%matrix(object$coefficients[1,],ncol=1)
    return(pred)
}
