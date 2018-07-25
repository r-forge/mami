# on loading the package
.onAttach <- function(libname = find.package("MAMI"), pkgname = "MAMI") {
packageStartupMessage("The package manual can be found at http://mami.r-forge.r-project.org/")
#r <- unclass(lsf.str(envir = asNamespace("MAMI"), all = T))
#r <- r[grep("^SL", r)]
#for(name in r){eval(parse(text=paste0(name, '<-MAMI:::', name)))}
}


# HELPER functions for mami()
ma <- function(maobject){
ctb <- summary(maobject)[[9]]
estimate <- ctb[,1][c(1,order(rownames(ctb))[-1])]
variance <- ctb[,2][c(1,order(rownames(ctb))[-1])]
importance <- maobject$importance[1:(length(estimate)-1)]
importance <- importance[order(names(importance))]
return(list(estimate,variance,importance))
}

ms.h <- function(msobject,fullnames,myfamily,keepnames,oc,add.parameters){
s.estimate <- msobject[,1]
s.se <- msobject[,2+as.numeric(myfamily=="cox")]
if((length(keepnames)-length(oc)+1+add.parameters)!=length(fullnames)){fullnames<-c("(Intercept)",fullnames)}
coefficientmat <- matrix(NA,nrow=1,ncol=length(fullnames))
variancemat    <- matrix(NA,nrow=1,ncol=length(fullnames))
colnames(coefficientmat) <- fullnames
colnames(variancemat) <- fullnames
for(i in 1:dim(coefficientmat)[2]){for(j in 1:length(s.estimate)){
if(colnames(coefficientmat)[i]==rownames(msobject)[j]){
coefficientmat[,i] <- s.estimate[j]
variancemat[,i] <- s.se[j]
}}}
for(i in 1:dim(variancemat)[1]){for(j in 1:dim(variancemat)[2]){
if(is.na(coefficientmat[i,j]==T)){coefficientmat[i,j]<-0}
if(is.na(variancemat[i,j]==T)){variancemat[i,j]<-0}}}
return(list(coefficientmat,variancemat))
}

ma.coxph.h <- function(maobject,id){
ctb <- summary(maobject)[[9]]
if(is.null(id)){
estimate <- c(0,ctb[,1][c(1,order(rownames(ctb))[-1])])
variance <- c(0,ctb[,2][c(1,order(rownames(ctb))[-1])])}else{
estimate <- c(0,ctb[,1][c(order(rownames(ctb)))])
variance <- c(0,ctb[,2][c(order(rownames(ctb)))])
}
names(estimate)[1] <- names(variance)[1] <- "(Intercept)"
importance <- (maobject$importance[substr(names(maobject$importance),1,8)!="cluster(" & substr(names(maobject$importance),1,7)!="strata("])#[1:(length(estimate)-1)]
importance <- c(importance[order(names(importance))])
return(list(estimate,variance,importance))
}

CV.cs <- function(mymodel,kf,cvr){
if(sum(mymodel$prior.weights!=1)!=0){stop("Cross validation should not be used with weights")}
cvIndex <- rep(1:kf,trunc(length(mymodel$y)/kf)+1)[1:length(mymodel$y)]
    if(cvr==TRUE){cvIndex <- sample(cvIndex)}
    cv <-0
          for (j in 1:kf){
                mymodelcv  <- glm(mymodel$formula, data=as.data.frame(mymodel$data[cvIndex!=j,]), family=mymodel$family)
                mypred   <- predict(mymodelcv,newdata=mymodel$data[cvIndex==j,])
                cv <- cv + mean((mymodel$data[cvIndex==j,]$y - mypred)^2)
                }
          sqrt(cv)
          }
          
CV.l <- function(mymodel,kf,cvr,X){
if(sum(weights(mymodel)!=1)!=0){stop("Cross validation should not be used with weights")}
    cvIndex <- as.numeric(factor(levels(getME(mymodel, name="flist")[[1]])))
    if(cvr==TRUE){cvIndex <- sample(cvIndex)}
    cvIndex2 <- cut_number(cvIndex,n=kf,labels = FALSE)
    cv <-0
    md1 <- as.data.frame(cbind(getME(mymodel, name="y"),getME(mymodel, name="X"),as.numeric(factor(getME(mymodel, name="flist")[[1]]))))
    names(md1)[1] <- names(X[[1]])[1]
    names(md1)[dim(md1)[2]] <- names(getME(mymodel, name="flist"))
          for (j in 1:kf){
                md2   <- md1[as.numeric(factor((getME(mymodel, name="flist")[[1]])))%in%cvIndex[cvIndex2!=j],]
                cvformula <- as.formula(formula(mymodel))
                mymodelcv  <- suppressWarnings(glmer(as.formula(cvformula, env=as.environment(as.data.frame(md2))), data=as.data.frame(md2), family=family(mymodel)$family))
                mypred   <- predict(mymodelcv,newdata=md1[as.numeric(factor((getME(mymodel, name="flist")[[1]])))%in%cvIndex[cvIndex2==j],],allow.new.levels = TRUE)
                cv <- cv + mean((getME(mymodel, name="y")[as.numeric(factor((getME(mymodel, name="flist")[[1]])))%in%cvIndex[cvIndex2==j]] - mypred)^2)
                }
          sqrt(cv)
          }

       
ma.bma.h <- function(bmaobject,model){
estimate <- bmaobject$postmean
variance <- bmaobject$postsd
importance <- bmaobject$probne0/100
if(model!="cox"){importance <- c(importance[order(names(importance))])
nms <- NULL
    for (i in 1:length(bmaobject$output.names)) {
        if (is.na(bmaobject$output.names[i][1]))
            nms <- c(nms, names(bmaobject$output.names[i]))
        else nms <- c(nms, paste(names(bmaobject$output.names[i]), unlist(bmaobject$output.names[i])[-1],
            sep = "."))
    }
}
if(model=="cox"){
nms <- NULL
    for (i in 1:length(bmaobject$output.names)) {
        if (is.na(bmaobject$output.names[i][1]))
            nms <- c(nms, names(bmaobject$output.names[i]))
        else nms <- c(nms, paste(names(bmaobject$output.names[i]), unlist(bmaobject$output.names[i])[-1],
            sep = "."))
    }
names(importance) <- bmaobject[[15]]
names(estimate) <-  nms
names(variance) <- nms
}
est2 <- bmaobject$mle[1,]
var2 <- bmaobject$se[1,]
names(est2) <-  nms
names(var2) <- nms
if(any(grepl("strata",nms))){
estimate <- estimate[!grepl("strata",nms)]
variance <- variance[!grepl("strata",nms)]
importance <- importance[!grepl("strata",names(importance))]
est2 <-  est2[!grepl("strata",nms)]
var2 <-  var2[!grepl("strata",nms)]
}
return(list(estimate,variance,importance,est2,var2))
}


# HELPER functions from other packages
round_any <- function(x, accuracy, f = round)
{
    f(x/accuracy)*accuracy
}

zero_range <- function(x, tol = 1000 * .Machine$double.eps)
{
    if (length(x) == 1)
        return(TRUE)
    if (length(x) != 2)
        stop("x must be length 1 or 2")
    if (any(is.na(x)))
        return(NA)
    if (x[1] == x[2])
        return(TRUE)
    if (all(is.infinite(x)))
        return(FALSE)
    m <- min(abs(x))
    if (m == 0)
        return(FALSE)
    abs((x[1] - x[2])/m) < tol
}

fullseq <- function(range, size, ..., pad = FALSE)
{
    if(zero_range(range))
        return(range + size*c(-1, 1)/2)
    x <- seq(round_any(range[1], size, floor), round_any(range[2],
        size, ceiling), by = size)
    if(pad){
        c(min(x) - size, x, max(x) + size)
    }
    else {
        x
    }
}

breaks <- function (x, equal, nbins = NULL, binwidth = NULL)
{
    equal <- match.arg(equal, c("numbers", "width"))
    if ((!is.null(nbins) && !is.null(binwidth)) || (is.null(nbins) &&
        is.null(binwidth))) {
        stop("Specify exactly one of n and width")
    }
    rng <- range(x, na.rm = TRUE, finite = TRUE)
    if(equal == "width"){
        if(!is.null(binwidth)) {
            fullseq(rng, binwidth)
        }
        else{
            seq(rng[1], rng[2], length.out = nbins + 1)
        }
    }
    else{
        if(!is.null(binwidth)){
            probs <- seq(0, 1, by = binwidth)
        }
        else{
            probs <- seq(0, 1, length.out = nbins + 1)
        }
        stats::quantile(x, probs, na.rm = TRUE)
    }
}

cut_number <- function(x, n = NULL, ...)
{
    brk <- breaks(x, "n", n)
    if (anyDuplicated(brk))
        stop("Insufficient data values to produce ", n, " bins.",
            call. = FALSE)
    cut(x, brk, include.lowest = TRUE, ...)
}

make.positive.definite <- function(m, tol)
{
    if (!is.matrix(m))
        m = as.matrix(m)
    d = dim(m)[1]
    if (dim(m)[2] != d)
        stop("Input matrix is not square!")
    es = eigen(m, symmetric = TRUE)
    esv = es$values
    if (missing(tol))
        tol = d * max(abs(esv)) * .Machine$double.eps
    delta = 2 * tol
    tau = pmax(0, delta - esv)
    dm = es$vectors %*% diag(tau, d) %*% t(es$vectors)
    return(m + dm)
}

###### SUPER LEARNER #######

SL.mma <- function(Y, X, newX, family,...){
if(all(Y%in%c(0,1))){stop("binomial outcome data not allowed")}
if(as.numeric(all(X[,1]==1))){X <- X[,-1]}
if(as.numeric(all(newX[,1]==1))){newX <- newX[,-1]}
mm <- mma(cbind(Y,X), formula = Y ~.)
mo <- mm$coefficients
newX <- newX[,colnames(newX)%in%colnames(X)]
pred <- predict.mma(mm,newdata=newX)
fit <- list(object = mm)
class(fit) <- 'SL.mma'
out <- list(pred = pred, fit = fit)
return(out)
}

predict.SL.mma <- function(object, newdata, ...){
pred <- predict.mma(object=object$object,newdata=newdata)
return(pred)
}

SL.mma.int <- function(Y, X, newX, family,...){
if(all(Y%in%c(0,1))){stop("binomial outcome data not allowed")}
if(as.numeric(all(X[,1]==1))){X <- X[,-1]}
if(as.numeric(all(newX[,1]==1))){newX <- newX[,-1]}
mm <- mma(cbind(Y,X), formula = Y ~.^2)
mo <- mm$coefficients
newX <- newX[,colnames(newX)%in%colnames(X)]
pred <- model.matrix(~.^2, data = newX)%*%matrix(mo[1,],ncol=1)
fit <- list(object = mm)
class(fit) <- 'SL.mma.int'
out <- list(pred = pred, fit = fit)
return(out)
}

predict.SL.mma.int <- function(object, newdata, ...){
pred <- predict.mma(object=object$object,newdata=newdata)
return(pred)
}

SL.mma2 <- function(Y, X, newX, family, cts.num=10, ...){
if(all(Y%in%c(0,1))){stop("binomial outcome data not allowed")}
if(as.numeric(all(X[,1]==1))){X <- X[,-1]}
if(as.numeric(all(newX[,1]==1))){newX <- newX[,-1]}
cts.x <- apply(X, 2, function(x) (length(unique(x)) > cts.num))
    if (sum(!cts.x) > 0) {
        mma.model <- paste("Y~",paste(colnames(X),
            collapse = "+"), "+", paste(paste("I(",
            colnames(X[, cts.x, drop = FALSE]),
            "^2)", sep = ""), collapse = "+"))
    }else{
        mma.model <- paste("Y~",paste(colnames(X),
            collapse = "+"), "+" ,paste(paste("I(",
            colnames(X[, cts.x, drop = FALSE]),
            "^2)", sep = ""), collapse = "+"))
    }
    if (sum(!cts.x) == length(cts.x)){
        mma.model <- paste("Y~", paste(colnames(X),
            collapse = "+"), sep = "")
    }
mm <- mma(cbind(Y,X), formula = as.formula(mma.model))
mo <- mm$coefficients
newX <- newX[,colnames(newX)%in%colnames(X)]
pred <- model.matrix(as.formula(gsub("Y~","~",mma.model)), data = newX)%*%matrix(mo[1,],ncol=1)
fit <- list(object = mm)
class(fit) <- 'SL.mma2'
out <- list(pred = pred, fit = fit)
return(out)
}

predict.SL.mma2 <- function(object, newdata, ...){
pred <- predict.mma(object=object$object,newdata=newdata)
return(pred)
}

SL.mma.all <- function(Y, X, newX, family,...){
if(all(Y%in%c(0,1))){stop("binomial outcome data not allowed")}
if(as.numeric(all(X[,1]==1))){X <- X[,-1]}
if(as.numeric(all(newX[,1]==1))){newX <- newX[,-1]}
jm <- jma(y=Y, x=X, method="MMA", subset="all")
mo <- jm$betahat
newX <- newX[,colnames(newX)%in%colnames(X)]
predic <- model.matrix(~.,data=newX)%*%matrix(mo,ncol=1)
fit <- list(object = jm)
class(fit) <- 'SL.mma.all'
out <- list(pred = predic, fit = fit)
return(out)
}

predict.SL.mma.all <- function(object, newdata, ...){
pred <- predict.jma(object=object$object,newdata=newdata)
return(pred)
}

SL.jma <- function(Y, X, newX, family,...){
if(all(Y%in%c(0,1))){stop("binomial outcome data not allowed")}
if(as.numeric(all(X[,1]==1))){X <- X[,-1]}
if(as.numeric(all(newX[,1]==1))){newX <- newX[,-1]}
jm <- jma(y=Y, x=X, method="JMA", subset="nested")
mo <- jm$betahat
newX <- newX[,colnames(newX)%in%colnames(X)]
predic <- model.matrix(~.,data=newX)%*%matrix(mo,ncol=1)
fit <- list(object = jm)
class(fit) <- 'SL.jma'
out <- list(pred = predic, fit = fit)
return(out)
}

predict.SL.jma <- function(object, newdata, ...){
pred <- predict.jma(object=object$object,newdata=newdata)
return(pred)
}

SL.jma.all <- function(Y, X, newX, family,...){
if(all(Y%in%c(0,1))){stop("binomial outcome data not allowed")}
if(as.numeric(all(X[,1]==1))){X <- X[,-1]}
if(as.numeric(all(newX[,1]==1))){newX <- newX[,-1]}
jm <- jma(y=Y, x=X, method="JMA", subset="all")
mo <- jm$betahat
newX <- newX[,colnames(newX)%in%colnames(X)]
predic <- model.matrix(~.,data=newX)%*%matrix(mo,ncol=1)
fit <- list(object = jm)
class(fit) <- 'SL.jma.all'
out <- list(pred = predic, fit = fit)
return(out)
}

predict.SL.jma.all <- function(object, newdata, ...){
pred <- predict.jma(object=object$object,newdata=newdata)
return(pred)
}

SL.jma.int <- function(Y, X, newX, family,...){
if(all(Y%in%c(0,1))){stop("binomial outcome data not allowed")}
if(as.numeric(all(X[,1]==1))){X <- X[,-1]}
if(as.numeric(all(newX[,1]==1))){newX <- newX[,-1]}
for(i in 1:dim(X)[2]){if(is.numeric(X[,i])){if(max(abs(X[,i]))>500){X[,i]<-scale(X[,i])}}}
for(i in 1:dim(newX)[2]){if(is.numeric(newX[,i])){if(max(abs(newX[,i]))>500){newX[,i]<-scale(newX[,i])}}}
jm <- jma(y=Y, x=model.matrix(~.^2, data = X)[,-1], method="JMA", subset="nested")
mo <- jm$betahat
newX <- newX[,colnames(newX)%in%colnames(X)]
predic <- model.matrix(~.^2, data = newX)%*%matrix(mo,ncol=1)
fit <- list(object = jm)
class(fit) <- 'SL.jma.int'
out <- list(pred = predic, fit = fit)
return(out)
}

predict.SL.jma.int <- function(object, newdata, ...){
pred <- predict.jma(object=object$object,newdata=newdata)
return(pred)
}

SL.jma2 <- function(Y, X, newX, family, cts.num=10,...){
if(all(Y%in%c(0,1))){stop("binomial outcome data not allowed")}
if(as.numeric(all(X[,1]==1))){X <- X[,-1]}
if(as.numeric(all(newX[,1]==1))){newX <- newX[,-1]}
for(i in 1:dim(X)[2]){if(is.numeric(X[,i])){if(max(abs(X[,i]))>500){X[,i]<-scale(X[,i])}}}
for(i in 1:dim(newX)[2]){if(is.numeric(newX[,i])){if(max(abs(newX[,i]))>500){newX[,i]<-scale(newX[,i])}}}
cts.x <- apply(X, 2, function(x) (length(unique(x)) > cts.num))
    if(sum(!cts.x) > 0) {
        jma.model <- paste("Y~",paste(colnames(X),
            collapse = "+"), "+", paste(paste("I(",
            colnames(X[, cts.x, drop = FALSE]),
            "^2)", sep = ""), collapse = "+"))
    }else{
        jma.model <- paste("Y~",paste(colnames(X),
            collapse = "+"), "+" ,paste(paste("I(",
            colnames(X[, cts.x, drop = FALSE]),
            "^2)", sep = ""), collapse = "+"))
    }
    if(sum(!cts.x) == length(cts.x)){
        jma.model <- paste("Y~", paste(colnames(X),
            collapse = "+"), sep = "")
    }
jm <- jma(y=Y, x=model.matrix(as.formula(jma.model), data = X)[,-1], method="JMA", subset="nested")
mo <- jm$betahat
newX <- newX[,colnames(newX)%in%colnames(X)]
predic <- model.matrix(as.formula(gsub("Y~","~",jma.model)), data = newX)%*%matrix(mo,ncol=1)
fit <- list(object = jm)
class(fit) <- 'SL.jma2'
out <- list(pred = predic, fit = fit)
return(out)
}

predict.SL.jma2 <- function(object, newdata, ...){
pred <- predict.jma(object=object$object,newdata=newdata)
return(pred)
}


SL.lae  <- function(Y, X, newX, family, ...){
if(family$family=="binomial" & !all(Y%in%c(0,1))){family$family="gaussian"}
if(as.numeric(all(X[,1]==1))){X <- X[,-1]}
if(as.numeric(all(newX[,1]==1))){newX <- newX[,-1]}
if(family$family=="gaussian"){
lam <- lae(cbind(Y,model.matrix(~.,data=X)[,-1]), kfold=10, glm.family = "gaussian", pd=F,factor.variables=NULL)
mo <- lam$coefficients[,1]
newX <- newX[,colnames(newX)%in%colnames(X)]
predic <- c(model.matrix(~.,data=newX)%*%matrix(mo,ncol=1))
fit <- list(object = lam)
}
if(family$family=="binomial"){
lam <- lae(cbind(Y,model.matrix(~.,data=X)[,-1]), kfold=10, glm.family = "binomial", pd=F,factor.variables=NULL)
mo <- lam$coefficients[,1]
newX <- newX[,colnames(newX)%in%colnames(X)]
predic <- c(plogis(model.matrix(~.,data=newX)%*%matrix(mo,ncol=1)))
fit <- list(object = lam)
}
class(fit) <- 'SL.lae'
out <- list(pred = predic, fit = fit)
return(out)
}

predict.SL.lae <- function(object, newdata, ...){
pred <- predict.lae(object=object$object,newdata=newdata)
return(pred)
}

SL.lae.int  <- function(Y, X, newX, family, ...){
if(family$family=="binomial" & !all(Y%in%c(0,1))){family$family="gaussian"}
if(as.numeric(all(X[,1]==1))){X <- X[,-1]}
if(as.numeric(all(newX[,1]==1))){newX <- newX[,-1]}
lam  <- lae(cbind(Y,model.matrix(~.^2, data = X)[,-1]), kfold=10, glm.family = family$family, pd=F, factor.variables = NULL)
mo <- lam$coefficients[,1]
newX <- newX[,colnames(newX)%in%colnames(X)]
predic <- model.matrix(~.^2, data = newX)%*%matrix(mo,ncol=1)
if(family$family=="binomial"){predic <- plogis(predic)}
fit <- list(object = lam)
class(fit) <- 'SL.lae.int'
out <- list(pred = predic, fit = fit)
return(out)
}

predict.SL.lae.int <- function(object, newdata,...){
pdata <- data.frame(model.matrix(as.formula(object$mf), data = data.frame(newdata))[,-1])
pred <- try(predict.lae(object=object$object,newdata=pdata),silent=TRUE)
if(class(pred)=="try-error"){pred <- try(predict.lae(object=object$object,newdata=pdata,clean=F))}
return(pred)
}

SL.lae2 <- function(Y, X, newX, family, cts.num=10,...){
if(family$family=="binomial" & !all(Y%in%c(0,1))){family$family="gaussian"}
if(as.numeric(all(X[,1]==1))){X <- X[,-1]}
if(as.numeric(all(newX[,1]==1))){newX <- newX[,-1]}
cts.x <- apply(X, 2, function(x) (length(unique(x)) > cts.num))
    if (sum(!cts.x) > 0) {
        lae.model <- paste("Y~",paste(colnames(X),
            collapse = "+"), "+", paste(paste("I(",
            colnames(X[, cts.x, drop = FALSE]),
            "^2)", sep = ""), collapse = "+"))
    }else{
        lae.model <- paste("Y~",paste(colnames(X),
            collapse = "+"), "+" ,paste(paste("I(",
            colnames(X[, cts.x, drop = FALSE]),
            "^2)", sep = ""), collapse = "+"))
    }
    if (sum(!cts.x) == length(cts.x)){
        lae.model <- paste("Y~", paste(colnames(X),
            collapse = "+"), sep = "")
    }
lam <- lae(cbind(Y,model.matrix(as.formula(lae.model), data = X)[,-1]), kfold=10, glm.family = family$family, pd=F, factor.variables = NULL)
mo  <- lam$coefficients[,1]
newX <- newX[,colnames(newX)%in%colnames(X)]
predic <- model.matrix(as.formula(gsub("Y~","~",lae.model)), data = newX)%*%matrix(mo,ncol=1)
if(family$family=="binomial"){predic <- plogis(predic)}
fit <- list(object = lam, mf = lae.model)
class(fit) <- 'SL.lae2'
out <- list(pred = predic, fit = fit)
return(out)
}

predict.SL.lae2 <- function(object, newdata,...){
pdata <- data.frame(model.matrix(as.formula(object$mf), data = data.frame(newdata))[,-1])
pred <- try(predict.lae(object=object$object,newdata=pdata),silent=TRUE)
if(class(pred)=="try-error"){pred <- try(predict.lae(object=object$object,newdata=pdata,clean=F))}
return(pred)
}
