###################################
# Lasso Averaging Estimaton (LAE) #
###################################

lae <- function(X, ycol=1, kfold=10, B.var=100,  calc.variance=FALSE, factor.variables=NULL, glm.family="gaussian",
 tries=10, standardize=TRUE, random=FALSE, pd=TRUE, ...){


# Preparation of datamatrix
if(glm.family%in%c("gaussian","binomial","poisson")==FALSE){stop("Choose glm.family as any of the following: 'gaussian','binomial','poisson'")}
if(pd==T & (glm.family=="binomial" | glm.family=="poisson")){cat("Note that the deviance is used as the loss function to select the tuning parameter for the LASSO estimate. \n However, the weights for Lasso averaging have to be determined based on asquared loss function.")}
X <- data.matrix(X)
keepnames <- colnames(X)
if(is.numeric(ycol[1])==FALSE){ycol <- match(ycol,colnames(X))}

tf <- function(mymatrix){colnames(mymatrix)[sapply(mymatrix,is.factor)]}
myindependent <- keepnames[ycol]
covariates <- keepnames[-ycol]
if(length(tf(X))>0){
factor.variables <- c(factor.variables,tf(X))
factor.variables <- intersect(factor.variables,factor.variables)
factor.variables <- intersect(colnames(X),factor.variables)
}
if(is.null(factor.variables)==FALSE){
if(pd==T){cat(paste("Note: The following variables are treated as factors and are recoded into dummies:", paste(factor.variables,collapse=" ")," \n\n"))}
X <- as.data.frame(X)
for(i in match(factor.variables,covariates)){
X[,covariates[i]] <- as.factor(X[,covariates[i]])
}
myformula <- as.formula(paste(myindependent,"~", paste(covariates, collapse="+")))
X <- cbind(X[ycol],model.matrix(myformula,data=as.data.frame(X))[,-1])
colnames(X)[1]<- keepnames[ycol]
ycol=1
X<- as.matrix(X)
}
if(standardize==TRUE){X[,-ycol] <- scale(X[,-ycol])}
n <- nrow(X)
m <- ncol(X)  


# these are our final estimates
shrink_sel     <- NULL       # Shrinkage estimation based on tuning parameter SELECTION
var_shrink_sel <- NULL
shrink_ave     <- NULL       # Shrinkage averaging estimation incorporating tuning parameter SELECTION UNCERTAINTY
var_shrink_ave <- NULL
# ...and this information will also be returned at the end of the function
sae_weights    <-NULL
sel_weights    <-NULL
complexity_parameter   <-NULL
variable_importance    <-NULL



#############
### Lasso ###
#############
cvIndex <- rep(1:kfold,trunc(n/kfold)+1)[1:n]
if(random==TRUE){cvIndex <- sample(cvIndex)}
lasso <- cv.glmnet(X[,-ycol],X[,ycol],foldid=cvIndex,nfolds=kfold,family=glm.family)
lasso_coeff <- glmnet::coef.glmnet(lasso, s = c(lasso$lambda))
complexity_parameter   <- lasso$lambda

# Shrinkage with tuning parameter selection based on CV_k 
        nolambda <- length(lasso$lambda)
        cvRMSE <- lasso$cvm  
        weights_minCV <-  c(rep(0, nolambda))
        myindex <- function(myvector){
        for(i in 1:length(myvector)){
        if(min(myvector)==myvector[i]){return(i)}}}
        weights_minCV[myindex(cvRMSE)] <- 1
        
# Shrinkage averaging estimation = OCV_k 
        
        # Shrinkage averaging estimation = OCV_k
        SAE<- function(avdata,mynewlambda,mykfold=kfold,ycol2=ycol,reshuffle=FALSE){

        nr <- nrow(avdata)
        newlambda<-mynewlambda
        MSPE <- 0 
        cvIndex3 <- rep(1:mykfold,trunc(nrow(avdata)/mykfold)+1)[1:nrow(avdata)]
        if(reshuffle==TRUE | random==TRUE){cvIndex3<-sample(cvIndex3)}
        
            for(j in 1:mykfold){    
        myavdata<- avdata[cvIndex3!=j,]
        cvIndex4 <- rep(1:mykfold,trunc(nrow(myavdata)/mykfold)+1)[1:nrow(myavdata)]
        myavlasso <-  cv.glmnet(myavdata[,-ycol2],myavdata[,ycol2],lambda=newlambda,foldid=cvIndex4,family=glm.family)
        Preavdata <- avdata[cvIndex3==j,]
        myavPred<-predict(myavlasso,newx=Preavdata[,-ycol2],s=newlambda)
        MSPE2 <- apply((avdata[cvIndex3==j,ycol2] - myavPred)^2,2,mean)
        MSPE <- MSPE + MSPE2 
        if(j==1){resmat <- matrix(NA,ncol=nrow(avdata),nrow=length(newlambda))}
        resmat[,cvIndex3==j] <- avdata[cvIndex3==j,ycol2] - myavPred
        }       

        E <- resmat%*%t(resmat)
        E <- make.positive.definite(E) # is normally not needed...just in case there is a problem
        Amat <- matrix(c(rep(1,length(newlambda))),nrow=1,ncol=length(newlambda))
        Amat <- t(rbind(Amat,diag(1,length(newlambda))))
        bvec <-c(1,rep(0,length(newlambda)))
        dvec <- c(rep(0,length(newlambda)))

        weight <- solve.QP(E, dvec, Amat, bvec, meq=1)$solution
        weight <- round(weight,digits=4)
        weight <- matrix(c(weight),nrow=1,ncol=length(newlambda))
        colnames(weight)<- paste("w",1:(length(newlambda)), sep="")
        rownames(weight)<- paste("OCV weights")
        return(weight)
        }

        sae_weights <- NULL
        tt<-1
        sae_weights <- suppressWarnings(try(SAE(X,lasso$lambda,kfold)))
        if(class(sae_weights)=="try-error"){
        while(tt<tries){
        cat(paste("Encountered problem when calculating the weights: re-shuffling cross validation sets now for the ",tt,". time... \n", sep=""))
        sae_weights <- try(SAE(X,lasso$lambda,kfold,reshuffle=TRUE))
        print(tt)
        if(class(sae_weights)!="try-error"){tt<-tries}else{tt<-tt+1}
        }}
        sel_weights <- weights_minCV
        # Calaculate estimates and their standard error
        shrink_sel  <- as.vector(weights_minCV%*%t(as.matrix(lasso_coeff)))
        shrink_ave  <- as.vector(sae_weights%*%t(as.matrix(lasso_coeff)))        
          
        # Calculate Variable Importance Measure
        index_imp <- matrix(rep(0,ncol(lasso_coeff)*nrow(lasso_coeff)),ncol=ncol(lasso_coeff),nrow=nrow(lasso_coeff))
        for(j in 1:ncol(lasso_coeff)){
        for(i in 1:nrow(lasso_coeff)){
        if(abs(lasso_coeff[i,j])>1e-05){index_imp[i,j]<-1}
        }}
        variable_importance<- index_imp%*%matrix(sae_weights,ncol=1, nrow=length(sae_weights))
        rownames(variable_importance) <- names(lasso_coeff)
        colnames(variable_importance) <- "Importance"

# Variance of lasso estimate via bootstrapping
if(calc.variance==TRUE){
boot.lasso <- function(mydata,indices){
mydata <- mydata[indices,]
myboot <- cv.glmnet(mydata[,-1],mydata[,1],nfolds=kfold,lambda=lasso$lambda,family=glm.family)
as.vector(cbind(glmnet::coef.glmnet(myboot, s = lasso$lambda),glmnet::coef.glmnet(myboot, s = myboot$lambda.min)))
}
mylasso_boot<- boot(cbind(X[,ycol],X[,-ycol]),boot.lasso,B.var)
var_lasso <-matrix(apply(mylasso_boot$t,2,sd),ncol=nolambda+1,nrow=nrow(lasso_coeff))
var_lasso1 <- var_lasso[,-dim(var_lasso)[2]]
var_lasso2 <- matrix(var_lasso[,dim(var_lasso)[2]],ncol=nolambda,nrow(lasso_coeff))
}

if(calc.variance==TRUE){
          av_1 <- matrix(c(rep(shrink_sel,each=nolambda)),ncol=length(shrink_sel),nrow=nolambda)
          av_3 <- matrix(c(rep(shrink_ave,each=nolambda)),ncol=length(shrink_ave),nrow=nolambda)
          var_shrink_sel  <- (weights_minCV%*%sqrt((t(as.matrix(lasso_coeff))-av_1)^2+t(var_lasso2)^2))
          var_shrink_ave   <- (sae_weights%*%sqrt((t(as.matrix(lasso_coeff))-av_3)^2+t(var_lasso1)^2))
          }else{
           var_shrink_sel <- rep(NA,length(shrink_sel))
           var_shrink_ave <- rep(NA,length(shrink_ave))
          }
        
# OLS estimate
mycovariables <- colnames(X)[-ycol]
myindependent <- colnames(X)[ycol]
my.formula     <- as.formula(paste(myindependent, "~", paste(mycovariables, collapse="+")))
ols           <- coefficients(glm(my.formula,data=as.data.frame(X),family=glm.family))
var_ols       <- summary(glm(my.formula, data=as.data.frame(X),family=glm.family))[[12]][,2]  

################################################################################

# Final results
mycoefficients <- t(rbind(shrink_ave,var_shrink_ave,shrink_sel,var_shrink_sel,ols,var_ols))
rownames(mycoefficients)<-rownames(lasso_coeff)
colnames(mycoefficients)<- c("LAE est","se","LASSO","se","(G)LM","se")
rownames(variable_importance)<-rownames(lasso_coeff)

################################################################################

# Return our results

res= list(coefficients=mycoefficients,
          variable.importance = variable_importance,
          sae.weights = sae_weights,
          sel.weights = sel_weights,
          complexity.parameter = complexity_parameter    
)

class(res) <- "lae"
res


}