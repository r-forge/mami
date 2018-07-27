predict.lae <- function(object,newdata=NULL,clean=T,...){   
    if(is.null(newdata)){newdata<-data.frame(object$setup[[1]])}
    if(clean==T){newdata <- newdata[,colnames(newdata)%in%colnames(object$setup[[1]])]}
    pred <- model.matrix(~.,data=newdata)%*%matrix(object$coefficients[,1],ncol=1)
    if(object[[6]][[2]]=="binomial"){pred <- plogis(pred)}
    return(pred)
}