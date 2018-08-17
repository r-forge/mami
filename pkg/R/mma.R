mma <- function(X, formula=NULL, ycol=1, variance=c("BA","boot"), bsa=200){

if(is.null(formula)){stop("Please provide formula for full model.")}
variance            <- match.arg(variance)
avdata <- as.data.frame(X)
full.model   <- lm(formula=formula, data=avdata)
n <- nrow(avdata)
p <- length(full.model$coefficients)
m <- length(attr(full.model$terms, "term.labels")) 
if(is.numeric(ycol[1])==FALSE){ycol <- match(ycol,colnames(X))}
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

weight <- quadprog::solve.QP(E, dvec, Amat, bvec, meq=1)$solution
weight <- matrix(c(weight),nrow=1,ncol=length(weight))
colnames(weight)<- paste("M.", 1:length(weight))
rownames(weight)<- paste("MMA weights")

estimates <- weight%*%coeffmat
colnames(estimates)<- names(full.model$coefficients)
rownames(estimates)<- paste("averaged est.")

avestimate <- matrix(c(rep(estimates,each=m+1)),ncol=p,nrow=m+1)

lce <- NULL
uce <- NULL

if(variance=="BA"){se <- (weight%*%sqrt((coeffmat-avestimate)^2+variancemat^2))}
if(variance=="boot"){
mma.se.boot <- function(mydata,indices){
mydata <- mydata[indices,]
mydata <- as.data.frame(mydata)
full.model.b   <- lm(formula=formula, data=mydata)
n.b <- nrow(mydata)
p.b <- length(full.model.b$coefficients)
m.b <- length(attr(full.model.b$terms, "term.labels")) 
if(is.numeric(ycol[1])==FALSE){ycol <- match(ycol,colnames(mydata))}
myy.b <- colnames(mydata)[ycol]

coeffmat.b <- matrix(0,ncol=p.b,nrow=m.b+1)
null.model.b <-  lm(mydata[,ycol]~1)
coeffmat.b[1,1]<- null.model.b$coefficients
resmat.b <- matrix(NA,ncol=m.b+1,nrow=n.b)
resmat.b[,1]<- null.model.b$residuals 
variancemat.b <- matrix(0,ncol=p.b,nrow=m.b+1)
variancemat.b[1,1]<-((summary(null.model.b))$coefficients)[,2]

for(j in 1:m.b){
        mycovariables.b <- paste(attr(full.model.b$terms, "term.labels")[1:j])
        myindependent.b <- myy.b
        myformula.b     <- as.formula(paste(myindependent.b,"~", paste(mycovariables.b, collapse="+") ))
        mymodel.b       <- lm(myformula.b, data=mydata)
        coeffmat.b[j+1,c(1:length(mymodel.b$coefficients))]    <- mymodel.b$coefficients
        resmat.b[,j+1]                <- mymodel.b$residuals
        variancemat.b[j+1,c(1:length(mymodel.b$coefficients))] <- ((summary(mymodel.b))$coefficients)[,2]
}

E.b <- t(resmat.b)%*%resmat.b   # this corresponds to Hansens ebar^T ebar
K.b <- apply((apply(coeffmat.b,c(1,2),function(myv){myv!=0})),1,sum)+1 # this corresponds to Hansens K-Vektor
shat.b <- (1/(n.b-(p.b+1)))*sum((full.model.b$residuals)^2) # this corresponds to Hansens estimated Variance
Amat.b <- matrix(c(rep(1,m.b+1)),nrow=1,ncol=m.b+1)
Amat.b <- t(rbind(Amat.b,diag(1,m.b+1)))
bvec.b <-c(1,rep(0,m.b+1))
dvec.b <- -shat.b*K.b

weight.b <- quadprog::solve.QP(E.b, dvec.b, Amat.b, bvec.b, meq=1)$solution
weight.b <- matrix(c(weight.b),nrow=1,ncol=length(weight.b))
estimates.b <- weight.b%*%coeffmat.b
c(estimates.b)
}
result.boot <- try(boot(X,mma.se.boot,bsa),silent=TRUE)
mysd <- function(myvalues){sd(na.omit(myvalues))}
se <- matrix(apply(result.boot$t,2,mysd),nrow=1)
lower95 <- function(xx){quantile(na.omit(xx),probs=0.025)}
upper95 <- function(xx){quantile(na.omit(xx),probs=0.975)}
lce <- matrix(apply(result.boot$t,2,lower95),nrow=1)
uce <- matrix(apply(result.boot$t,2,upper95),nrow=1)
rownames(lce)<- paste("95% CI (lower)")
rownames(uce)<- paste("95% CI (upper)")
}

rownames(se) <- paste("standard error")
colnames(se)<- names(full.model$coefficients)


res <- list(coefficients=rbind(estimates,se,lce,uce),
            averaging.weights = weight,setup=list(formula,X[,-ycol]))
class(res) <- "mma"
res
 

}



