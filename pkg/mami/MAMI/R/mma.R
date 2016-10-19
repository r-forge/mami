mma <- function(somedata, modelformula=NULL, ycol=1){

if(is.null(modelformula)){stop("Please provide formula for full model.")}

avdata <- as.data.frame(somedata)
full.model   <- lm(modelformula, data=avdata)
n <- nrow(avdata)
p <- length(full.model$coefficients)
m <- length(attr(full.model$terms, "term.labels")) 
if(is.numeric(ycol[1])==FALSE){ycol <- match(ycol,colnames(somedata))}
myy <- colnames(avdata)[ycol]

coeffmat <- matrix(0,ncol=p,nrow=m+1)
null.model <-  lm(avdata[,ycol]~1)
coeffmat[1,1]<- null.model$coefficients
resmat <- matrix(NA,ncol=m+1,nrow=n)
resmat[,1]<- null.model$residuals 
variancemat <- matrix(0,ncol=p,nrow=m+1)
variancemat[1,1]<-((summary(null.model))$coefficients)[,2]

for(i in 1:m){
        mycovariables <- paste(attr(full.model$terms, "term.labels")[1:i])
        myindependent <- myy
        myformula     <- as.formula(paste(myindependent,"~", paste(mycovariables, collapse="+") ))
        mymodel       <- lm(myformula, data=avdata)
        coeffmat[i+1,c(1:length(mymodel$coefficients))]    <- mymodel$coefficients
        resmat[,i+1]                <- mymodel$residuals
        variancemat[i+1,c(1:length(mymodel$coefficients))] <- ((summary(mymodel))$coefficients)[,2]
}

E <- t(resmat)%*%resmat   # this corresponds to Hansens ebar^T ebar
K <- apply((apply(coeffmat,c(1,2),function(myv){myv!=0})),1,sum)+1 # this corresponds to Hansens K-Vektor
shat <- (1/(n-(p+1)))*sum((full.model$residuals)^2) # this corresponds to Hansens estimated Variance
Amat <- matrix(c(rep(1,m+1)),nrow=1,ncol=m+1)
Amat <- t(rbind(Amat,diag(1,m+1)))
bvec <-c(1,rep(0,m+1))
dvec <- -shat*K

weight <- solve.QP(E, dvec, Amat, bvec, meq=1)$solution
weight <- matrix(c(weight),nrow=1,ncol=length(weight))
colnames(weight)<- paste("M.", 1:length(weight))
rownames(weight)<- paste("MMA weights")

estimates <- weight%*%coeffmat
colnames(estimates)<- names(full.model$coefficients)
rownames(estimates)<- paste("averaged est.")

avestimate <- matrix(c(rep(estimates,each=m+1)),ncol=p,nrow=m+1)
se <- (weight%*%sqrt((coeffmat-avestimate)^2+variancemat^2))
rownames(se) <- paste("standard error")
colnames(se)<- names(full.model$coefficients)


res <- list(coefficients=rbind(estimates,se),
            averaging.weights = weight)
class(res) <- "mma"
res
 

}



