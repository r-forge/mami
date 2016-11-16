#################################
# Lasso Averaging Estimaton #
#################################

lae <- function(X,ycol=1,B.var=50,nolambda=100,kfold=5, my.formula=NULL, standardize=TRUE, calc.variance=TRUE){


# Preparation of datamatrix
X <- as.data.frame(X)
keepnames <- colnames(X)
if(is.numeric(ycol)==FALSE){ycol<-sum(as.numeric(names(X)==ycol)*seq(1:length(names(X))))}
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

################################################################################
myindex <- function(myvector){
for(i in 1:length(myvector)){
if(min(myvector)==myvector[i]){return(i)}
}}

################################################################################



#############
### Lasso ###
#############

mybound <- (1:nolambda)/(nolambda)
mycovariables <- names(X)[-ycol]
myindependent <- names(X)[ycol]
if(is.null(my.formula)){my.formula     <- as.formula(paste(myindependent, "~", paste(mycovariables, collapse="+")))}
lasso <- suppressWarnings(l1ce(my.formula, X, bound=mybound))
lasso_coeff <- coefficients(lasso)
lasso_res <- residuals(lasso)
gcv_inf <- gcv(lasso, type="Tibshirani")
complexity_parameter <- mybound

# Variance of lasso estimate via bootstrapping
if(calc.variance==TRUE){
boot.lasso <- function(mydata,indices){
mydata <- mydata[indices,]
myboot <- suppressWarnings(l1ce(my.formula, data=mydata, bound=mybound))
coefficients(myboot)
}
mylasso_boot<- boot(as.data.frame(X),boot.lasso,B.var)
var_lasso <-matrix(apply(mylasso_boot$t,2,sd),ncol=ncol(lasso_coeff),nrow=nrow(lasso_coeff))
mean_lasso <-matrix(apply(mylasso_boot$t,2,mean),ncol=ncol(lasso_coeff),nrow=nrow(lasso_coeff))
}

# OLS estimate
ols           <- coefficients(lm(my.formula, data=X))
var_ols       <- summary(lm(my.formula, data=X))[[4]][,2]

                  
                  
        ##################
        # OCV estimators #
        ##################
        
        
        #  OCV Shrinkage averaging estimation
        
        ## Cross Validation
        candidateBounds <- gcv_inf[,1]
        cvRMSE <- rep(0,nolambda)

        cvIndex <- rep(1:kfold,trunc(n/kfold)+1)[1:n]
        set.seed(666)
        cvIndex <- sample(cvIndex)
          for (j in 1:kfold) {
            for (k in (1:nolambda)) {
        mylassodata<- X[cvIndex!=j,]
        mylassodata <- as.data.frame(mylassodata)
        colnames(mylassodata)<- names(X)
        mylasso<-suppressWarnings(l1ce(my.formula,data=mylassodata,bound=candidateBounds[k]))
        Predata <- X[cvIndex==j,]
        colnames(Predata) <- colnames(mylassodata)
        myPred<-predict(mylasso,newdata=Predata)
        cvRMSE[k] <- cvRMSE[k] + mean((X[cvIndex==j,ycol] - myPred)^2)
          }}

        cvRMSE <- cvRMSE/kfold
        # Shrinkage with tuning parameter selection based on CV_k 
        weights_minCV <-  c(rep(0, nolambda))
        weights_minCV[myindex(cvRMSE)] <- 1
        
        # Shrinkage averaging estimation = OCV_k
        jma<- function(avdata,mynewlambda,mykfold,myformula2=my.formula){

        n <- nrow(avdata)
        resmat <- matrix(NA,ncol=nrow(avdata),nrow=mynewlambda)
        mysupernewbound <- (1:nolambda)/(nolambda)

          for(k in 1:mynewlambda){
            for(j in 1:mykfold){

        cvRMSE <- rep(0,mynewlambda)
        cvIndex <- rep(1:mykfold,trunc(n/mykfold)+1)[1:n]
        mylassodata<- avdata[cvIndex!=j,]
        mylasso<-suppressWarnings(l1ce(myformula2,data=mylassodata,bound=mysupernewbound[k]))
        Predata <- avdata[cvIndex==j,]
        colnames(Predata) <- colnames(mylassodata)
        myPred<-predict(mylasso,newdata=Predata)
        cvRMSE[k] <- cvRMSE[k] + mean((avdata[cvIndex==j,ycol] - myPred)^2)
        resmat[k,cvIndex==j] <- avdata[cvIndex==j,ycol] - myPred
        }}

        E <- resmat%*%t(resmat)
        E <- make.positive.definite(E) # is normally not needed...just in case there is a problem
        Amat <- matrix(c(rep(1,mynewlambda)),nrow=1,ncol=mynewlambda)
        Amat <- t(rbind(Amat,diag(1,mynewlambda)))
        bvec <-c(1,rep(0,mynewlambda))
        dvec <- c(rep(0,mynewlambda))

        weight <- solve.QP(E, dvec, Amat, bvec, meq=1)$solution
        weight <- round(weight,digits=4)
        weight <- matrix(c(weight),nrow=1,ncol=mynewlambda)
        colnames(weight)<- paste("w",1:(mynewlambda), sep="")
        rownames(weight)<- paste("OCV weights")
        return(weight)
        }

        sae_weights <- jma(X,nolambda,kfold)
        sel_weights <- weights_minCV
        # Calaculate estimates and their standard error
        shrink_sel  <- weights_minCV%*%lasso_coeff
        shrink_ave  <- sae_weights%*%lasso_coeff
        
        if(calc.variance==TRUE){
          av_1 <- matrix(c(rep(shrink_sel,each=nolambda)),ncol=length(shrink_sel),nrow=nolambda)
          av_3 <- matrix(c(rep(shrink_ave,each=nolambda)),ncol=length(shrink_ave),nrow=nolambda)
          var_shrink_sel  <- (weights_minCV%*%sqrt((lasso_coeff-av_1)^2+var_lasso^2))
          var_shrink_ave   <- (sae_weights%*%sqrt((lasso_coeff-av_3)^2+var_lasso^2))
          }else{
           var_shrink_sel <- rep(NA,length(shrink_sel))
           var_shrink_ave <- rep(NA,length(shrink_ave))
          }
          
        # Calculate Variable Importance Measure
        index_imp <- matrix(rep(0,ncol(lasso_coeff)*nrow(lasso_coeff)),ncol=ncol(lasso_coeff),nrow=nrow(lasso_coeff))
        for(j in 1:ncol(lasso_coeff)){
        for(i in 1:nrow(lasso_coeff)){
        if(abs(lasso_coeff[i,j])>1e-05){index_imp[i,j]<-1}
        }}
        variable_importance<- t(index_imp)%*%matrix(sae_weights,ncol=1, nrow=length(sae_weights))
        rownames(variable_importance) <- names(ols)
        colnames(variable_importance) <- "Importance"
        
   

################################################################################

# Final results
mycoefficients <- t(rbind(shrink_ave,var_shrink_ave,shrink_sel,var_shrink_sel,ols,var_ols))
rownames(mycoefficients)<-names(ols)
colnames(mycoefficients)<- c("LAE est","se","LASSO","se","OLS","se")
rownames(variable_importance)[1]<-"Interc."

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