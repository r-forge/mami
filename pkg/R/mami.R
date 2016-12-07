###################################################################
# Model selection and model averaging (after multiple imputation) #
###################################################################


mami <- function(X,method=c("criterion.average","criterion.selection","MMA","LASSO/LAE"), 
  criterion=c("AIC","BIC","BIC+","CV","GCV"), B=20, X.org=NULL,inference=c("standard","+bootstrapping"),
  missing.data=c("imputed","none","CC"), var.remove=NULL, user.weights=NULL, candidate.models=c("all",
  "restricted","very restricted"), model.family=c("gaussian","binomial","poisson","coxph"),
  add.factor=NULL, add.interaction=NULL, add.transformation=NULL, ycol=1, CI=0.95, kfold=5, id=NULL,
  use.stratum=NULL, report.exp=FALSE, print.time=FALSE, print.warnings=TRUE, ...){

#library(boot, quietly=TRUE, warn.conflicts=FALSE)
#library(quadprog, quietly=TRUE, warn.conflicts=FALSE)
#library(survival, quietly=TRUE, warn.conflicts=FALSE)
#library(MuMIn, quietly=TRUE, warn.conflicts=FALSE)
#library(MASS, quietly=TRUE, warn.conflicts=FALSE)
#library(zoo, quietly=TRUE, warn.conflicts=FALSE)
#library(lme4, quietly=TRUE, warn.conflicts=FALSE)
#library(lasso2, quietly=TRUE, warn.conflicts=FALSE)
#library(corpcor, quietly=TRUE, warn.conflicts=FALSE) 
#library(BMA, quietly=TRUE, warn.conflicts=FALSE)      

ptm <- proc.time()

method            <- match.arg(method)
criterion         <- match.arg(criterion)
inference         <- match.arg(inference)
missing.data      <- match.arg(missing.data)
model.family      <- match.arg(model.family)
candidate.models  <- match.arg(candidate.models)

# User information
if(method=="MMA" & model.family!="gaussian"){stop("Mallows model averaging only works for the linear model.")}
if(print.warnings==TRUE & method=="MMA"){cat("Note: The results of MMA depend on the ordering of the regressors. \n")}
if(print.warnings==TRUE & missing.data=="CC"){cat("Note: you should use complete cases only if you are sure that this is appropriate! \n")}
if(model.family=="quasipoisson"){stop("Quasi-likelihood methods are not supported yet.")}
if(inference=="+bootstrapping" & is.null(X.org) & missing.data=="imputed"){stop("Please also provide unimputed data under `X.org' when using bootstrapping.")}
if(print.warnings==TRUE & method=="criterion.selection" & inference=="standard"){cat("Note: Bootstrap confidence intervals are preferred after model selection. Standard confidence intervals after model selection may be overoptimistic. \n")}
if(method=="criterion.selection" & is.null(id)==FALSE & model.family!="coxph" & (criterion!="CV"&criterion!="GCV")){stop("Model selection (with stepAIC) for mixed models not possible.\n However, you can use method=`criterion.average' instead which also produces results for model selection,\n though definition of AIC in mixed models is not aways clear.")}
if(is.null(user.weights)==FALSE & inference=="+bootstrapping"){stop("Bootstrapping not allowed when specifying weights.")}
if(method=="criterion.selection" & criterion=="BIC" & model.family=="coxph"){stop("BIC not allowed for Cox model. `N' is not clearly defined.")}
if(method=="LASSO/LAE" & is.null(add.factor)==FALSE){stop("Lasso Averaging Estimation for factor variables not implemented yet.")}
if(method=="LASSO/LAE" & model.family!="gaussian"){stop("Only linear model is allowed for LASSO estimation.")}
if(method=="LASSO/LAE" & is.null(id)==FALSE){stop("LASSO estimation not available for longitudinal data.")}
if(report.exp==TRUE & model.family=="gaussian"){warning("Exponentiated coefficients are not useful for linear models.")}
if(is.null(use.stratum)==FALSE & model.family!="coxph"){stop("You can only specify strata in a Cox Model.")}
if(missing.data=="CC" & any(is.na(X))==FALSE & print.warnings==TRUE){stop("You specified a complete case analysis but there is no missing data.")}
if(model.family=="coxph" & (criterion=="CV" | criterion=="GCV")){stop("Cross Validation not allowed in the Cox PH model.")}
if((criterion=="CV" | criterion=="GCV") & method=="criterion.average"){stop("Use method='criterion.selection' for model selection with (G)CV")}
if(criterion=="CV" & kfold==1){stop("Choose kfold>1")}
if(criterion=="CV" & print.warnings==TRUE){cat("Note: a much faster option is leave-one-out cross validation approximated by `criterion=GCV' \n")} 
if(criterion=="CV" & is.null(id)==FALSE & is.null(add.transformation)==FALSE){stop("It is currently not possible to add transformations to mixed models when criterion='CV'")}
if(print.warnings==TRUE & criterion=="CV" & is.null(id)==FALSE){cat("Note that using cross validation in mixed models does not make full use of the random intercept estimates. \n")}
if(is.null(add.interaction)==FALSE & is.list(add.interaction)==FALSE){stop("Use `list()' to add interactions")}
if(print.warnings==TRUE & (is.null(add.transformation)==FALSE | is.null(add.interaction)==FALSE) & (length(add.interaction)+length(add.transformation)>=4)){cat("Note: Many interactions and/or transformation increase model complexity and computation time. \n If the full model can not be fit because of overfitting, mami may crash.\n")}
if(print.warnings==TRUE & criterion=="BIC+" & method=="criterion.selection"){stop("You can use `criterion=BIC' for model selection. `BIC+' is needed for model averaging, particularly when analysing a big dataset.")}
if(print.warnings==TRUE & criterion=="BIC+" & (candidate.models=="restricted" | candidate.models=="very restricted")){cat("Note: candidate models are restricted by the leaps algorithm of the `BMA' package. \n")}
if(criterion=="BIC+" & is.null(id)==FALSE){stop("BIC+ can only be used with cross-sectional data, but you specified an `id'.")}
if(is.null(id)==FALSE & inference=="+bootstrapping"){cat("Note: You utilized the longitudinal bootstrap. \n")}

# Preparation of datamatrix
M = NULL          # Number of imputations
if(class(X)=="amelia"){
M=length(X$imputations)
X_ <- X$imputations
am <- X$arguments
}
if(class(X)=="data.frame"){
M=1
X_ <- list(X)
}
if(class(X)=="list"){
if(inference=="+bootstrapping"){stop("To utilize bootstrapping please impute with Amelia II.")}
M= length(X)
X_ <- X
}
if(class(X)=="mids"){stop("mice not implemented yet")}
if(is.null(M)){
M=length(X)
X_ <- X
}
if(is.null(M)){stop("Data format not recognized. Choose between dataframe or multiple imputation via Amelia II  or a list of user-specified imputed datasets.")}
X <- X_
if(is.numeric(ycol[1])==FALSE){ycol <- match(ycol,colnames(X[[1]]))}
if(is.numeric(var.remove)==FALSE){var.remove <- match(var.remove,colnames(X[[1]]))}
if(any(!(ycol %in% 1))==TRUE | is.null(var.remove)==FALSE){for(m in 1:M){X[[m]] <- subset(X[[m]],select=c(names(X[[m]])[ycol],names(X[[m]])[-c(ycol,var.remove)])) }}
keepnames <- colnames(X[[1]])
if((candidate.models=="restricted" & (dim(X[[1]])[2]-length(ycol))<5) |  (candidate.models=="very restricted" & (dim(X[[1]])[2]-length(ycol))<10)){stop("You should not restrict your candidate models for this analysis")}
if(missing.data=="CC"){X <- list(na.omit(as.data.frame(X)))}
if(missing.data!="imputed"){X.org<-X[[1]]}    
if(missing.data=="none"     & any(is.na(X[[1]]))==TRUE){stop("There is still missing data but you specified there is none.")}
if(missing.data=="imputed"  & any(is.na(X[[1]]))==TRUE){stop("There is still missing data but you specified the data is imputed.")}
if(missing.data=="none" & length(X)>1){stop("You specified there is no missing data but you provide multiple datasets as `X'.")}
#options(na.action = "na.fail")
          
# these are our final model averaging estimates we want to report or keep
est                     <- NULL         # final estimator(s)
se_est                  <- NULL         # se of final estimator
lower                   <- NULL         # lower confidence limit
upper                   <- NULL         # upper confidence limit
ests_boot               <- NULL         # all bootstrap estimates
est_boot                <- NULL         # final bootstrap estimator(s) based on re-imputed data
est_boot2               <- NULL         # final bootstrap estimator(s) based on mean of bootstrap samples
se_est_boot             <- NULL         # se of final bootstrap estimator
lower_boot              <- NULL         # lower bootstrap confidence limit
upper_boot              <- NULL         # upper bootstrap confidence limit
variable_importance     <- NULL         # based on model averaging weights (if applicable)
# these are our final model selection estimates we want to report or keep
est.s                     <- NULL         # final estimator(s)
se_est.s                  <- NULL         # se of final estimator
lower.s                   <- NULL         # lower confidence limit
upper.s                   <- NULL         # upper confidence limit
ests_boot.s               <- NULL         # all bootstrap estimates
est_boot.s                <- NULL         # final bootstrap estimator(s) based on re-imputed data
est_boot2.s               <- NULL         # final bootstrap estimator(s) based on mean of bootstrap samples
se_est_boot.s             <- NULL         # se of final bootstrap estimator
lower_boot.s              <- NULL         # lower bootstrap confidence limit
upper_boot.s              <- NULL         # upper bootstrap confidence limit
#
result.boot.s             <- NULL
result.boot.ma            <- NULL

##############################################
## Inference for model selection/averaging ###
##############################################

# Regression formula
tf <- function(mymatrix){colnames(mymatrix)[sapply(mymatrix,is.factor)]}
if(length(tf(X[[1]]))>0){
add.factor <- c(add.factor,tf(X[[1]]))
add.factor <- intersect(add.factor,add.factor)
add.factor <- intersect(colnames(X[[1]]),add.factor)
}
if(print.warnings==TRUE & is.null(add.factor)==FALSE){cat(paste("Note: The following variables are treated as factors:", paste(add.factor,collapse=" "),"\n"))}
myfamily <- model.family
oc <- c(1:length(ycol))
myindependent <- keepnames[oc]
usenames <- keepnames
fn <- function(myname){paste("as.factor(",myname,")",sep="")}
add.parameters <- 0
if(is.null(add.factor)==FALSE){for(i in match(add.factor,usenames)){
add.parameters <- add.parameters + length(levels(as.factor(X[[1]][,usenames[i]]))) - 2
usenames[i] <- fn(usenames[i])
}}
param.transform <- NULL
if(is.null(add.transformation)==FALSE){for(i in 1:length(add.transformation)){
param.transform <- paste("I(",add.transformation,")",collapse="+",sep="")}
add.parameters <- add.parameters+length(add.transformation)}
param.interaction <- NULL
if(is.null(add.factor)){add.factor<-NA}
if(is.null(add.interaction)==FALSE){for(j in 1:length(add.interaction)){
  add.factor.param <- 1
  if(is.na(any(match(add.factor,add.interaction[[j]])))){add.parameters <- add.parameters+1}
  if(is.na(any(match(add.factor,add.interaction[[j]])))==FALSE & length(na.omit(match(add.factor,add.interaction[[j]])))==1){add.parameters <- add.parameters + length(levels(as.factor((X[[1]])[,(add.interaction[[j]][na.omit(match(add.factor,add.interaction[[j]]))])])))-1} 
  if(is.na(any(match(add.factor,add.interaction[[j]])))==FALSE & length(na.omit(match(add.factor,add.interaction[[j]])))> 1){for(k in 1:length(na.omit(match(add.factor,add.interaction[[j]])))){add.factor.param <- add.factor.param*(length(levels(as.factor((X[[1]])[,(add.interaction[[j]][na.omit(match(add.factor,add.interaction[[j]]))[k]])])))-1)} 
  add.parameters   <- add.parameters + add.factor.param} 
}}
interact <- function(a){paste(a[1],":",a[2],sep="")}
if(is.null(add.interaction)==FALSE){
for(j in 1:length(add.interaction)){if(is.na(any(match(add.factor,add.interaction[[j]])))==FALSE){
  if(length(na.omit(match(add.factor,add.interaction[[j]])))==1){add.interaction[[j]][na.omit(match(add.factor,add.interaction[[j]]))] <- fn(add.interaction[[j]][na.omit(match(add.factor,add.interaction[[j]]))])}
  if(length(na.omit(match(add.factor,add.interaction[[j]])))> 1){for(k in 1:length(na.omit(match(add.factor,add.interaction[[j]])))){add.interaction[[j]][k] <- fn(add.interaction[[j]][k])}}}
param.interaction <- paste(lapply(add.interaction,interact),collapse="+",sep="")
}}
saveid <- id
if(is.null(id)==FALSE){
if(model.family!="coxph"){
usenames <- usenames[-match(id,usenames)]
id <- paste("(1|",id,")",sep="")
add.parameters <- add.parameters-1} else
{usenames <- usenames[-match(id,usenames)]
id <- paste("cluster(",id,")",sep="")
add.parameters <- add.parameters-1} 
}
if(is.null(use.stratum)==FALSE){usenames <- usenames[-match(use.stratum,usenames)]
use.stratum <- paste("strata(",use.stratum,")",sep="")
add.parameters <- add.parameters-1
}
if(is.null(user.weights)==FALSE){usenames <- usenames[-match(user.weights,usenames)]
add.parameters <- add.parameters-1
} else{user.weights<-""} 
mycovariables <- paste(c(usenames[-oc],param.transform,param.interaction,id),collapse="+")
if(length(ycol) > 1){
myindependent <- paste("Surv(",paste(keepnames[oc],collapse=","),")",sep="")
mycovariables <- paste(c(usenames[-oc],param.transform,param.interaction,use.stratum,id),collapse="+")
}
myformula <- paste(myindependent,"~", mycovariables)


# We need these matrices later on...
mycoefflist <- matrix(NA,ncol=M,nrow=(length(keepnames)-length(ycol)+1+add.parameters))
myvarlist   <- matrix(NA,ncol=M,nrow=(length(keepnames)-length(ycol)+1+add.parameters))
mycoefflist.s <- matrix(NA,ncol=M,nrow=(length(keepnames)-length(ycol)+1+add.parameters))
myvarlist.s <- matrix(NA,ncol=M,nrow=(length(keepnames)-length(ycol)+1+add.parameters))
myvarimplist <- matrix(NA,ncol=M,nrow=length(c(usenames[-oc],param.transform,param.interaction)))  
if(criterion=="BIC+" & model.family=="coxph"){
mycoefflist <- matrix(NA,ncol=M,nrow=(length(keepnames)-length(ycol)+1+add.parameters-1))
myvarlist   <- matrix(NA,ncol=M,nrow=(length(keepnames)-length(ycol)+1+add.parameters-1))
mycoefflist.s <- matrix(NA,ncol=M,nrow=(length(keepnames)-length(ycol)+1+add.parameters-1))
myvarlist.s <- matrix(NA,ncol=M,nrow=(length(keepnames)-length(ycol)+1+add.parameters-1))
myvarimplist <- matrix(NA,ncol=M,nrow=length(c(usenames[-oc],param.transform,param.interaction))-1)  
}

###########################
# Method: Model selection #
###########################

if(method=="criterion.selection" & (criterion!="CV" & criterion!="GCV")){

ptm1.5 <- proc.time()
if(model.family=="coxph"){mycoefflist.s <- mycoefflist.s[-1,]
                          myvarlist.s <- myvarlist.s[-1,]}        
for(m in 1:M)try({
if(print.warnings==TRUE & missing.data=="imputed"){cat("Busy with model selection in imputed dataset no.",m, "\n")}
if(model.family!="coxph"){full.model <- glm(as.formula(myformula), data = as.data.frame(X[[m]]),   weights=X[[m]]$user.weights, family=model.family)} else{
  if(model.family=="coxph"){full.model <- coxph(as.formula(myformula), data = as.data.frame(X[[m]]),   weights=X[[m]]$user.weights)}}
sel_m <- stepAIC(full.model, trace=0, direction="both", k=2+(as.numeric(criterion=="BIC")*(-2+log(dim(X[[m]])[1]))))
if(model.family!="coxph"){sel_coef <- summary(sel_m)[[12]]} else{
   if(model.family=="coxph"){sel_coef <- summary(sel_m)[[7]]}}
mycoef   <- matrix(c(rep(0,length(coefficients(full.model)))),ncol=1,nrow=length(coefficients(full.model)),dimnames=list(names(coefficients(full.model)),"coeff."))
mystd    <- matrix(c(rep(0,length(coefficients(full.model)))),ncol=1,nrow=length(coefficients(full.model)),dimnames=list(names(coefficients(full.model)),"coeff."))
for(i in 1:length(mycoef)){for(j in 1:length(sel_coef[,1])){
if(rownames(mycoef)[i]==rownames(sel_coef)[j]){
mycoef[i] <- sel_coef[j,1]
if(model.family!="coxph"){mystd[i] <- sel_coef[j,2]} else{
if(model.family=="coxph"){mystd[i] <- sel_coef[j,3]}}
}}}
mycoefflist.s[,m] <- mycoef
myvarlist.s[,m]   <- mystd
})
rownames(mycoefflist.s) <- names(coefficients(full.model)) 
rownames(myvarlist.s) <- names(coefficients(full.model)) 


ptm2 <-  proc.time()
inbetweentime <- round((((ptm2-ptm1.5)/60))[1],digits=2)
if(print.time==TRUE & inference=="+bootstrapping"){cat(paste("The analysis will run for about another", inbetweentime*(B+1), "minute(s) \n"))}

if(inference=="+bootstrapping"){
if(print.warnings==TRUE){cat("Now start with bootstrap estimation... \n" )}

mydata1=NULL
mydata2=NULL
if(is.null(id)==FALSE){ # Prepare longitudinal bootstrap if needed!
flag1 <- rep(1,dim(X.org)[1])
for(i in 2:length(flag1)){
if((X.org[,saveid])[i]==(X.org[,saveid])[i-1]){flag1[i]<-0}
}
flag1 <- as.logical(flag1)
mydata1 <- X.org[flag1,]
mydata2 <- X.org
}

msmi.boot <- function(mydata,indices,mydata1=mydata1,mydata2=mydata2){
mydata <- mydata[indices,]

  if(is.null(mydata1)==FALSE){     # Longitudinal Bootstrap if needed!
  our.ids      <- mydata1[,saveid]
  b.sample     <- sample(our.ids, length(our.ids), replace = TRUE)
  startsample  <- rep(NA,dim(mydata2)[2])
  b.loop       <- 1
  b.samples    <- as.numeric(rownames(table(b.sample)))
  while(b.loop < max(table(b.sample))+5){
  startsample  <- rbind(startsample,mydata2[(c(mydata2[,saveid]) %in% c(b.samples)),])
  b.samples    <- as.numeric(rownames(table(b.sample)[table(b.sample) > b.loop]))
  b.loop       <- b.loop+1
  }
  startsample<-startsample[-1,]
  flag2 <- rep(1,dim(startsample)[1])
  for(i in 2:length(startsample[,saveid])){
  if((startsample[,saveid])[i]==(startsample[,saveid])[i-1]){flag2[i]<-0}
  }
  flag2 <- cumsum(flag2)
  startsample[,saveid] <- flag2
  mydata <- startsample
  rownames(mydata) <- seq(1:dim(mydata)[1])
  }

    if(missing.data=="imputed"){mydata.out <-  amelia(mydata, m = M, p2s = as.numeric(print.warnings==TRUE),frontend = FALSE, idvars = am$idvars,
                                  ts = am$ts, cs = am$cs, polytime = am$polytime, splinetime = am$splinetime, intercs = am$intercs,
                                  lags = am$lags, leads = am$leads, startvals = am$startvals, tolerance = am$tolerance,
                                  logs = am$logs, sqrts = am$sqrts, lgstc = am$lgstc, noms = am$noms, ords = am$ords,
                                  incheck = TRUE, collect = FALSE, arglist = am$arglist, empri = am$empri,
                                  priors = am$priors, autopri = am$autopri, emburn = c(0,0), bounds = am$bounds,
                                  max.resample = am$max.resample, overimp = am$overimp)} 
    if(missing.data!="imputed"){mydata.out <- NULL
                                mydata.out$imputations[[1]] <- mydata} 
if(any(!(ycol %in% 1))==TRUE | is.null(var.remove)==FALSE){for(m in 1:M){mydata.out$imputations[[m]] <- try(subset(mydata.out$imputations[[m]],select=c(names(mydata.out$imputations[[m]])[ycol],names(mydata.out$imputations[[m]])[-c(ycol,var.remove)])))}}
mycoefflist_b <- matrix(NA,ncol=M,nrow=(length(keepnames)-length(ycol)+1+add.parameters))
if(model.family=="coxph"){mycoefflist_b <- mycoefflist_b[-1,]} 
for(m in 1:M)try({ 
if(print.warnings==TRUE & missing.data=="imputed"){cat("Busy with model selection in imputed dataset no.",m, "\n")}
if(model.family!="coxph"){full.model.b <- glm(as.formula(myformula), data = as.data.frame(mydata.out$imputations[[m]]), family=model.family)} else{
  if(model.family=="coxph"){full.model.b <- coxph(as.formula(myformula), data = as.data.frame(mydata.out$imputations[[m]]))}}
sel_m.b <- try(stepAIC(full.model.b, trace=0, direction="both", k=2+(as.numeric(criterion=="BIC")*(-2+log(dim(X[[m]])[1])))))
if(model.family!="coxph"){sel_coef.b <- summary(sel_m.b)[[12]]} else{
   if(model.family=="coxph"){sel_coef.b <- summary(sel_m.b)[[7]]}}
mycoef.b   <- matrix(c(rep(0,length(coefficients(full.model)))),ncol=1,nrow=length(coefficients(full.model)),dimnames=list(names(coefficients(full.model)),"coeff."))
for(i in 1:length(mycoef)){for(j in 1:length(sel_coef.b[,1])){
if(rownames(mycoef)[i]==rownames(sel_coef.b)[j]){
mycoef.b[i] <- sel_coef.b[j,1]
}}}
if(length(mycoefflist_b[,m])==length(mycoef.b)){mycoefflist_b[,m] <- mycoef.b}
})
c(apply(mycoefflist_b,1,mean))
}
result.boot <- try(boot(X.org,msmi.boot,B,mydata1=mydata1,mydata2=mydata2))
if(any(is.na(result.boot$t))){warning("Sorry: Some bootstrap samples have been removed as model selection did not work there well enough.")}
result.boot.s <- result.boot$t[apply(apply(result.boot$t,1,is.na),2,any)==F,]
ests_boot.s <- result.boot.s
est_boot.s <-  result.boot$t0
est_boot2.s <- apply(result.boot.s,2,mean)
se_est_boot.s <- apply(result.boot.s,2,sd)
nquantile <- function(myvector){quantile(myvector,probs=c(((1-CI)/2),1-((1-CI)/2)))}
boot.quantiles <- apply(result.boot.s,2,nquantile)
lower_boot.s <- boot.quantiles[1,]
upper_boot.s<- boot.quantiles[2,]
}

}

###############
# Method: MMA #
###############

if(method=="MMA"){

ptm1.5 <- proc.time()
for(m in 1:M){
if(print.warnings==TRUE & missing.data=="imputed"){cat("Busy with Mallow's model averaging in imputed dataset no.",m, "\n")}
av_m <- mma(X[[m]],myformula,...)
mycoefflist[,m] <- av_m$coefficients[1,]
myvarlist[,m]   <- av_m$coefficients[2,]
}
rownames(mycoefflist) <- names(av_m$coefficients[1,])
rownames(myvarlist)     <- names(av_m$coefficients[1,])


ptm2 <-  proc.time()
inbetweentime <- round((((ptm2-ptm1.5)/60))[1])
if(print.time==TRUE & inference=="+bootstrapping"){cat(paste("The analysis will run for about another", inbetweentime*(B+1), "minute(s) \n"))}

if(inference=="+bootstrapping"){
if(print.warnings==TRUE){cat("Now start with bootstrap estimation... \n" )}
mami.boot <- function(mydata,indices){
mydata <- mydata[indices,]
  if(missing.data=="imputed"){mydata.out <-  amelia(mydata, m = M, p2s = as.numeric(print.warnings==TRUE),frontend = FALSE, idvars = am$idvars,
                                  ts = am$ts, cs = am$cs, polytime = am$polytime, splinetime = am$splinetime, intercs = am$intercs,
                                  lags = am$lags, leads = am$leads, startvals = am$startvals, tolerance = am$tolerance,
                                  logs = am$logs, sqrts = am$sqrts, lgstc = am$lgstc, noms = am$noms, ords = am$ords,
                                  incheck = TRUE, collect = FALSE, arglist = am$arglist, empri = am$empri,
                                  priors = am$priors, autopri = am$autopri, emburn = c(0,0), bounds = am$bounds,
                                  max.resample = am$max.resample, overimp = am$overimp)} 
  if(missing.data!="imputed"){mydata.out <- NULL
                              mydata.out$imputations[[1]] <- mydata} 
if(any(!(ycol %in% 1))==TRUE | is.null(var.remove)==FALSE){for(m in 1:M){mydata.out$imputations[[m]] <- try(subset(mydata.out$imputations[[m]],select=c(names(mydata.out$imputations[[m]])[ycol],names(mydata.out$imputations[[m]])[-c(ycol,var.remove)])))}}
mycoefflist_b <- matrix(NA,ncol=M,nrow=(length(keepnames)-length(ycol)+1+add.parameters))
for(m in 1:M)try({  
av_m.b <- try(mma(mydata.out$imputations[[m]],myformula,...))
if(length(mycoefflist_b[,m]) == length(av_m.b$coefficients[1,])){mycoefflist_b[,m] <- av_m.b$coefficients[1,]}
})
c(apply(mycoefflist_b,1,mean))
}
result.boot <- try(boot(X.org,mami.boot,B))
if(any(is.na(result.boot$t))){warning("Sorry: Some bootstrap samples have been removed as MMA did not work there.")}
result.boot.ma <- result.boot$t[apply(apply(result.boot$t,1,is.na),2,any)==F,]
ests_boot <- result.boot.ma
est_boot <- result.boot$t0
est_boot2 <- apply(result.boot.ma,2,mean)
se_est_boot <- apply(result.boot.ma,2,sd)
nquantile <- function(myvector){quantile(myvector,probs=c(((1-CI)/2),1-((1-CI)/2)))}
boot.quantiles <- apply(result.boot.ma,2,nquantile)
lower_boot <- boot.quantiles[1,]
upper_boot <- boot.quantiles[2,]
}

}

#################
# Method: Lasso #
#################

if(method=="LASSO/LAE"){

ptm1.5 <- proc.time()

for(m in 1:M){
if(print.warnings==TRUE & missing.data=="imputed"){cat("Busy with LASSO (averaging) estimation in imputed dataset no.",m, "\n")}
av_m <- lae(X[[m]],my.formula=as.formula(myformula),kfold=kfold,...)
mycoefflist[,m] <- av_m$coefficients[,1]
myvarlist[,m]   <- av_m$coefficients[,2]
mycoefflist.s[,m] <-  av_m$coefficients[,3]
myvarlist.s[,m]   <-  av_m$coefficients[,4]
myvarimplist[,m] <- av_m$variable.importance[-1]
}
rownames(mycoefflist) <- names(av_m$coefficients[,1])
rownames(myvarlist)     <- names(av_m$coefficients[,1])
rownames(mycoefflist.s) <- names(av_m$coefficients[,1])
rownames(myvarlist.s)     <- names(av_m$coefficients[,1])
rownames(myvarimplist)     <- names(av_m$coefficients[,1])[-1]

ptm2 <-  proc.time()
inbetweentime <- round((((ptm2-ptm1.5)/60))[1],digits=2)
if(print.time==TRUE & inference=="+bootstrapping"){cat(paste("The analysis will run for about another", inbetweentime*(B+1), "minute(s) \n"))}

if(inference=="+bootstrapping"){
if(print.warnings==TRUE){cat("Now start with bootstrap estimation... \n" )}
mami.boot <- function(mydata,indices){
mydata <- mydata[indices,]
  if(missing.data=="imputed"){mydata.out <-  amelia(mydata, m = M, p2s = as.numeric(print.warnings==TRUE),frontend = FALSE, idvars = am$idvars,
                                  ts = am$ts, cs = am$cs, polytime = am$polytime, splinetime = am$splinetime, intercs = am$intercs,
                                  lags = am$lags, leads = am$leads, startvals = am$startvals, tolerance = am$tolerance,
                                  logs = am$logs, sqrts = am$sqrts, lgstc = am$lgstc, noms = am$noms, ords = am$ords,
                                  incheck = TRUE, collect = FALSE, arglist = am$arglist, empri = am$empri,
                                  priors = am$priors, autopri = am$autopri, emburn = c(0,0), bounds = am$bounds,
                                  max.resample = am$max.resample, overimp = am$overimp)} 
  if(missing.data!="imputed"){mydata.out <- NULL
                              mydata.out$imputations[[1]] <- mydata} 
if(any(!(ycol %in% 1))==TRUE | is.null(var.remove)==FALSE){for(m in 1:M){mydata.out$imputations[[m]] <- try(subset(mydata.out$imputations[[m]],select=c(names(mydata.out$imputations[[m]])[ycol],names(mydata.out$imputations[[m]])[-c(ycol,var.remove)])))}}
mycoefflist_b <- matrix(NA,ncol=M,nrow=(length(keepnames)-length(ycol)+1+add.parameters))
mycoefflist_b.s <- matrix(NA,ncol=M,nrow=(length(keepnames)-length(ycol)+1+add.parameters))
for(m in 1:M)try({ 
if(print.warnings==TRUE & missing.data=="imputed"){cat("Busy with LASSO (averaging) estimation in imputed dataset no.",m, "\n")} 
av_m.b <- lae(mydata.out$imputations[[m]],my.formula=as.formula(myformula),kfold=kfold,...) 
if(length(mycoefflist_b[,m]) == length(av_m.b$coefficients[,1])){
mycoefflist_b[,m] <- av_m.b$coefficients[,1]
mycoefflist_b.s[,m] <- av_m.b$coefficients[,3]
}})
c(apply(mycoefflist_b.s,1,mean),apply(mycoefflist_b,1,mean))
}
result.boot <- try(boot(X.org,mami.boot,B))
if(any(apply(apply(result.boot$t,1,is.na),2,any))==TRUE & print.warnings==TRUE){cat("There occurred fitting problems in one or many bootstrap samples. The relevant samples have been removed. \n")}
result.boot.s <- result.boot$t[apply(apply(result.boot$t[,1:(length(keepnames)-length(ycol)+1+add.parameters)],1,is.na),2,any)==F,1:(length(keepnames)-length(ycol)+1+add.parameters)]    
ests_boot.s <- result.boot.s
est_boot.s <- result.boot$t0[1:(length(keepnames)-length(ycol)+1+add.parameters)]
est_boot2.s <- apply(result.boot.s,2,mean)
se_est_boot.s <- apply(result.boot.s,2,sd)
nquantile <- function(myvector){quantile(myvector,probs=c(((1-CI)/2),1-((1-CI)/2)))}
boot.quantiles.s <- apply(result.boot.s,2,nquantile)
lower_boot.s <- boot.quantiles.s[1,]
upper_boot.s <- boot.quantiles.s[2,]
result.boot.ma <- result.boot$t[apply(apply(result.boot$t[,c((length(keepnames)-length(ycol)+1+1+add.parameters):dim(result.boot$t)[2])],1,is.na),2,any)==F,c((length(keepnames)-length(ycol)+1+1+add.parameters):dim(result.boot$t)[2])]
ests_boot <- result.boot.ma
est_boot <- result.boot$t0[c((length(keepnames)-length(ycol)+1+1+add.parameters):dim(result.boot$t)[2])]
est_boot2 <- apply(result.boot.ma,2,mean)
se_est_boot <- apply(result.boot.ma,2,sd)
boot.quantiles <- apply(result.boot.ma,2,nquantile)
lower_boot <- boot.quantiles[1,]
upper_boot <- boot.quantiles[2,]
} 

}


#############################################
# Method:  Model averaging, Akaike weights  #
#############################################


if(method=="criterion.average" | (method=="criterion.selection" & (criterion=="CV"|criterion=="GCV"))){

mymax <- round(length(keepnames[-ycol])/(as.numeric(candidate.models==c("all"))+ 2*as.numeric(candidate.models==c("restricted"))+4*as.numeric(candidate.models==c("very restricted"))), digits=0)
#mydelta <- round(25/(as.numeric(candidate.models==c("all"))+ 2*as.numeric(candidate.models==c("restricted"))+5*as.numeric(candidate.models==c("very restricted"))), digits=0)
if((candidate.models=="restricted" | candidate.models=="very restricted") & print.warnings==TRUE & criterion!="BIC+"){cat(paste("Note: Due to your specified restriction the maximum amount of variables allowed in the candidate models is ",mymax,".\n If your restriction is too tight `mami' may fail. \n",sep=""))} 
 

ma <- function(maobject){
ctb <- summary(maobject)[[9]]
estimate <- ctb[,1][c(1,order(rownames(ctb))[-1])]
variance <- ctb[,2][c(1,order(rownames(ctb))[-1])]
importance <- maobject$importance[1:(length(estimate)-1)]
importance <- importance[order(names(importance))]
return(list(estimate,variance,importance))
}

ms <- function(msobject,fullnames,myfamily){
s.estimate <- msobject[,1]
s.se <- msobject[,2+as.numeric(myfamily=="coxph")] 
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

if(model.family=="coxph"){
ma.coxph <- function(maobject){
ctb <- summary(maobject)[[9]]
estimate <- c(0,ctb[,1][order(rownames(ctb))])
variance <- c(0,ctb[,2][order(rownames(ctb))])
importance <- (maobject$importance[substr(names(maobject$importance),1,8)!="cluster(" & substr(names(maobject$importance),1,7)!="strata("])[1:(length(estimate)-1)]
importance <- c(importance[order(names(importance))])
return(list(estimate,variance,importance))
}
}

if(criterion=="CV" & is.null(id)==TRUE){
CV <- function(mymodel,kf=kfold){
if(sum(mymodel$prior.weights!=1)!=0){stop("Cross validation should not be used with weights")}
cvIndex <- rep(1:kf,trunc(length(mymodel$y)/kf)+1)[1:length(mymodel$y)]
    cvIndex <- sample(cvIndex)
    cv <-0
          for (j in 1:kf){
                mymodelcv  <- glm(mymodel$formula, data=as.data.frame(mymodel$data[cvIndex!=j,]), family=mymodel$family)
                mypred   <- predict(mymodelcv,newdata=mymodel$data[cvIndex==j,])
                cv <- cv + mean((mymodel$data[cvIndex==j,]$y - mypred)^2)
                }
          sqrt(cv)
          }
}
if(criterion=="CV" & is.null(id)==FALSE){
CV <- function(mymodel,kf=kfold){
if(sum(weights(mymodel)!=1)!=0){stop("Cross validation should not be used with weights")}
    cvIndex <- as.numeric(factor(levels(getME(mymodel, name="flist")[[1]])))
    cvIndex <- sample(cvIndex)
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
}
if(criterion=="GCV"){
GCV <- function(gmodel){
N<- length(residuals(gmodel))
if(is.null(id)==TRUE){myy<-gmodel$y}
if(is.null(id)==FALSE){myy<-getME(mymodel, name="y")}
1/N*sum(((myy-predict(gmodel))/(1-sum(hatvalues(gmodel))))^2)
}
}


if(criterion=="BIC+"){
ma.bma <- function(bmaobject){
estimate <- bmaobject$postmean
variance <- bmaobject$postsd
importance <- bmaobject$probne0/100
if(model.family!="coxph"){importance <- c(importance[order(names(importance))])}
if(model.family=="coxph"){
names(importance) <-  bmaobject[[15]]
names(estimate) <-bmaobject[[15]]
names(variance) <-bmaobject[[15]]
}
return(list(estimate,variance,importance))
}
ma.bms <- function(bmaobject){
estimate <- bmaobject$mle[1,]
variance <- bmaobject$se[1,]
return(list(estimate,variance))
}
}


ptm1.5 <- proc.time()
myil <- rep(list(NA),M)

for(m in 1:M)try({
if(print.warnings==TRUE & missing.data=="imputed" & (criterion!="CV" & criterion!="GCV")){cat("Busy with model averaging in imputed dataset no.",m, "\n")} 
if(print.warnings==TRUE & missing.data=="imputed" & (criterion=="CV" | criterion=="GCV")){cat("Busy with model selection in imputed dataset no.",m, "\n")}
if(criterion!="BIC+"){
if(model.family!="coxph" & is.null(id)){mymodel <- glm(as.formula(myformula), data = as.data.frame(X[[m]]), na.action=na.fail, weights=X[[m]]$user.weights, family=model.family)} 
if(model.family!="coxph" & is.null(id)==FALSE){mymodel <- suppressWarnings(glmer(as.formula(myformula), na.action=na.fail, data = as.data.frame(X[[m]]),   weights=X[[m]]$user.weights, family=model.family))}  # suppress warnings since we use glmer for gaussian family
if(model.family=="coxph"){mymodel <- coxph(as.formula(myformula), data = as.data.frame(X[[m]]), na.action=na.fail, weights=X[[m]]$user.weights)}
mymo <- suppressMessages(suppressWarnings(dredge(mymodel, m.lim=c(0,mymax),rank=get(criterion),fixed=c(id,use.stratum),...)))
if(print.warnings==TRUE & (criterion!="CV" & criterion!="GCV") & m==M){cat(paste("Model averaging was based on",dim(mymo)[1],"candidate models \n"))}
if(print.warnings==TRUE & (criterion=="CV" | criterion=="GCV") & m==M){cat(paste("Model selection was based on",dim(mymo)[1],"candidate models \n"))}
if(any(is.na(mymo[,"delta"]))==FALSE | model.family=="coxph"){mym  <- suppressWarnings(get.models(mymo,subset=NA))} else {stop("Problems in model averaging function `dredge'. Possibly try `criterion='BIC+' '.")}
if(model.family!="coxph"){av_m <- suppressWarnings(ma(model.avg(mym)))} else{
       if(model.family=="coxph"){av_m <- try(suppressWarnings(ma.coxph(model.avg(mym))))}}
if(!(model.family!="coxph" & is.null(id)==FALSE)){sel_m <- ms(coefficients(summary(mym[[1]])),names(coefficients(mymodel)),model.family)} else{
                                                  sel_m <- ms(coefficients(summary(mym[[1]])),rownames(summary(mymodel)[[10]]),model.family)             }
  } else {
  if(criterion=="BIC+" & model.family!="coxph"){mymodel   <- bic.glm(as.formula(myformula), data = as.data.frame(X[[m]]), na.action=na.fail, wt=X[[m]]$user.weights, glm.family=model.family,...)}
  if(criterion=="BIC+" & model.family=="coxph"){mymodel   <- bic.surv(as.formula(myformula), data = as.data.frame(X[[m]]), na.action=na.fail, wt=X[[m]]$user.weights,...)}
  av_m  <- ma.bma(mymodel)
  sel_m <- ma.bms(mymodel)
  if(print.warnings==TRUE){cat(paste("Model averaging was based on",length(mymodel$bic),"candidate models \n"))}
  }
  
mycoefflist[,m] <- av_m[[1]]
myvarlist[,m]   <- av_m[[2]]
mycoefflist.s[,m] <- sel_m[[1]]
myvarlist.s[,m]   <- sel_m[[2]]

myvarimp <- matrix(av_m[[3]],nrow=1,ncol=length(av_m[[3]]),dimnames=list(NULL,names(av_m[[3]])))
myvarimp <- myvarimp[,sort(colnames(myvarimp))]
myil[[m]] <- myvarimp
})

nimp <- max(unlist(lapply(myil, length)))
ii <- 0
for(ii in 1:M){
ind <- as.numeric(length(myil[[ii]])==nimp)
if(ind==TRUE){cn <- names(myil[[ii]])}}
myvarimplist <- matrix(0,ncol=M,nrow=nimp,dimnames=list(cn,NULL)) 
for(i in 1:M){myvarimplist[(is.na(match(rownames(myvarimplist),names(myil[[i]])))==FALSE),i]<-myil[[i]]}

ptm2 <-  proc.time()
inbetweentime <- round((((ptm2-ptm1.5)/60))[1],digits=2)
if(print.time==TRUE & inference=="+bootstrapping"){cat(paste("The analysis will run for about another", inbetweentime*(B+1), "minute(s) \n"))}
if(!(model.family!="coxph" & is.null(id)==FALSE)){mynames <- names(mymodel$coefficients)} else{
mynames<-rownames(summary(mymodel)[[10]])}
if((length(keepnames)-length(oc)+1+add.parameters)!=length(mynames)){mynames<-c("(Intercept)",mynames)}

if(inference=="+bootstrapping"){
if(print.warnings==TRUE){cat("Now start with bootstrap estimation... \n" )}

mydata1=NULL
mydata2=NULL
if(is.null(id)==FALSE){ # Prepare longitudinal bootstrap if needed!
flag1 <- rep(1,dim(X.org)[1])
for(i in 2:length(flag1)){
if((X.org[,saveid])[i]==(X.org[,saveid])[i-1]){flag1[i]<-0}
}
flag1 <- as.logical(flag1)
mydata1 <- X.org[flag1,]
mydata2 <- X.org
}

mami.boot <- function(mydata,indices){
mydata <- mydata[indices,]
  if(is.null(mydata1)==FALSE){     # Longitudinal Bootstrap if needed!
  our.ids      <- mydata1[,saveid]
  b.sample     <- sample(our.ids, length(our.ids), replace = TRUE)
  startsample  <- rep(NA,dim(mydata2)[2])
  b.loop       <- 1
  b.samples    <- as.numeric(rownames(table(b.sample)))
  while(b.loop < max(table(b.sample))+5){
  startsample  <- rbind(startsample,mydata2[(c(mydata2[,saveid]) %in% c(b.samples)),])
  b.samples    <- as.numeric(rownames(table(b.sample)[table(b.sample) > b.loop]))
  b.loop       <- b.loop+1
  }
  startsample<-startsample[-1,]
  flag2 <- rep(1,dim(startsample)[1])
  for(i in 2:length(startsample[,saveid])){
  if((startsample[,saveid])[i]==(startsample[,saveid])[i-1]){flag2[i]<-0}
  }
  flag2 <- cumsum(flag2)
  startsample[,saveid] <- flag2
  mydata <- startsample
  rownames(mydata) <- seq(1:dim(mydata)[1])
  }
  if(missing.data=="imputed"){mydata.out <-  amelia(mydata, m = M, p2s = as.numeric(print.warnings==TRUE),frontend = FALSE, idvars = am$idvars,
                                  ts = am$ts, cs = am$cs, polytime = am$polytime, splinetime = am$splinetime, intercs = am$intercs,
                                  lags = am$lags, leads = am$leads, startvals = am$startvals, tolerance = am$tolerance,
                                  logs = am$logs, sqrts = am$sqrts, lgstc = am$lgstc, noms = am$noms, ords = am$ords,
                                  incheck = TRUE, collect = FALSE, arglist = am$arglist, empri = am$empri,
                                  priors = am$priors, autopri = am$autopri, emburn = c(0,0), bounds = am$bounds,
                                  max.resample = am$max.resample, overimp = am$overimp, boot.type = am$boot.type)} 
  if(missing.data!="imputed"){mydata.out <- NULL
                              mydata.out$imputations[[1]] <- mydata} 
if(any(!(ycol %in% 1))==TRUE | is.null(var.remove)==FALSE){for(m in 1:M){mydata.out$imputations[[m]] <- try(subset(mydata.out$imputations[[m]],select=c(names(mydata.out$imputations[[m]])[ycol],names(mydata.out$imputations[[m]])[-c(ycol,var.remove)])))}}
if(criterion!="BIC+"){
mycoefflist_b   <- matrix(NA,ncol=M,nrow=(length(keepnames)-length(ycol)+1+add.parameters))
mycoefflist_b.s <- matrix(NA,ncol=M,nrow=(length(keepnames)-length(ycol)+1+add.parameters))} else {
if(model.family!="coxph"){mycoefflist_b      <- matrix(NA,ncol=M,nrow=(length(keepnames)-length(ycol)+1+add.parameters))
                          mycoefflist_b.s    <- matrix(NA,ncol=M,nrow=(length(keepnames)-length(ycol)+1+add.parameters)) }
if(model.family=="coxph"){mycoefflist_b      <- matrix(NA,ncol=M,nrow=(length(keepnames)-length(ycol)  +add.parameters))
                          mycoefflist_b.s    <- matrix(NA,ncol=M,nrow=(length(keepnames)-length(ycol)  +add.parameters)) }
}
for(m in 1:M)try({  # try() is used here to let the bootstrap run even when a single bootstrap sample causes the specified model to crash. Will later be removed.
if(print.warnings==TRUE & missing.data=="imputed" & (criterion!="CV" & criterion!="GCV")){cat("Busy with model averaging in imputed dataset no.",m, "\n")} 
if(print.warnings==TRUE & missing.data=="imputed" & (criterion=="CV" | criterion=="GCV")){cat("Busy with model selection in imputed dataset no.",m, "\n")}
if(criterion!="BIC+"){
if(model.family!="coxph" & is.null(id)){       mymodel.b <- try(glm(as.formula(myformula), na.action=na.fail, data = as.data.frame(mydata.out$imputations[[m]]), family=model.family))} 
if(model.family!="coxph" & is.null(id)==FALSE){mymodel.b <- suppressWarnings(glmer(as.formula(myformula), na.action=na.fail, data = as.data.frame(mydata.out$imputations[[m]]), family=model.family))}  # suppress warnings since we use glmer for gaussian family
if(model.family=="coxph"){                     mymodel.b <- try(coxph(as.formula(myformula), na.action=na.fail, data = as.data.frame(mydata.out$imputations[[m]])))}
mym.b <- suppressMessages(suppressWarnings(get.models(dredge(mymodel.b, m.lim=c(0,mymax),rank=get(criterion),fixed=c(id,use.stratum),...),subset=NA)))
if(model.family!="coxph"){av_m.b <- try(suppressWarnings(ma(model.avg(mym.b))))} else {
       if(model.family=="coxph"){av_m.b <- try(suppressWarnings(ma.coxph(model.avg(mym.b))))}}
if(!(model.family!="coxph" & is.null(id)==FALSE)){sel_mb <- try(ms(coefficients(summary(mym.b[[1]])),mynames,model.family))} else{
                                                  sel_mb <- try(ms(coefficients(summary(mym.b[[1]])),mynames,model.family))}
  } else {
  if(criterion=="BIC+" & model.family!="coxph"){mymodel.b   <- try(bic.glm(as.formula(myformula), data = as.data.frame(mydata.out$imputations[[m]]), na.action=na.fail, glm.family=model.family,...)) }
  if(criterion=="BIC+" & model.family=="coxph"){mymodel.b   <- try(bic.surv(as.formula(myformula), data = as.data.frame(mydata.out$imputations[[m]]), na.action=na.fail,...))}
  av_m.b  <- try(ma.bma(mymodel.b))
  sel_mb <- try(ma.bms(mymodel.b))
  }
if(length(mycoefflist[,m])==length(av_m[[1]])){mycoefflist[,m] <- av_m[[1]]}
if(length(mycoefflist_b[,m])==length(av_m.b[[1]])){mycoefflist_b[,m] <- av_m.b[[1]]}
if(length(mycoefflist_b.s[,m])==length(sel_mb[[1]])){mycoefflist_b.s[,m] <- sel_mb[[1]]}
})
c(apply(mycoefflist_b.s,1,mean),apply(mycoefflist_b,1,mean))
}
result.boot <- try(boot(X.org,mami.boot,B))
if(any(apply(apply(result.boot$t,1,is.na),2,any))==TRUE & print.warnings==TRUE & method!="criterion.selection"){cat("There occurred fitting problems in one or many bootstrap samples. The relevant samples have been removed. \n")}
if(criterion=="BIC+" & model.family=="coxph"){ result.boot.s <- result.boot$t[apply(apply(result.boot$t[,1:(length(keepnames)-length(ycol)  +add.parameters)],1,is.na),2,any)==F,1:(length(keepnames)-length(ycol)  +add.parameters)] } else {  
                                               result.boot.s <- result.boot$t[apply(apply(result.boot$t[,1:(length(keepnames)-length(ycol)+1+add.parameters)],1,is.na),2,any)==F,1:(length(keepnames)-length(ycol)+1+add.parameters)]    }
ests_boot.s <- result.boot.s
if(criterion=="BIC+" & model.family=="coxph"){ est_boot.s <- result.boot$t0[1:(length(keepnames)-length(ycol)  +add.parameters)]} else {
                                               est_boot.s <- result.boot$t0[1:(length(keepnames)-length(ycol)+1+add.parameters)]}
est_boot2.s <- apply(result.boot.s,2,mean)
se_est_boot.s <- apply(result.boot.s,2,sd)
nquantile <- function(myvector){quantile(myvector,probs=c(((1-CI)/2),1-((1-CI)/2)))}
boot.quantiles.s <- apply(result.boot.s,2,nquantile)
lower_boot.s <- boot.quantiles.s[1,]
upper_boot.s <- boot.quantiles.s[2,]
if(criterion=="BIC+" & model.family=="coxph"){result.boot.ma <- result.boot$t[apply(apply(result.boot$t[,c((length(keepnames)-length(ycol)+1  +add.parameters):dim(result.boot$t)[2])],1,is.na),2,any)==F,c((length(keepnames)-length(ycol)  +1+add.parameters):dim(result.boot$t)[2])]} else {
                                              result.boot.ma <- result.boot$t[apply(apply(result.boot$t[,c((length(keepnames)-length(ycol)+1+1+add.parameters):dim(result.boot$t)[2])],1,is.na),2,any)==F,c((length(keepnames)-length(ycol)+1+1+add.parameters):dim(result.boot$t)[2])]}
ests_boot <- result.boot.ma
if(criterion=="BIC+" & model.family=="coxph"){ est_boot <- result.boot$t0[c((length(keepnames)-length(ycol)+1  +add.parameters):dim(result.boot$t)[2])]} else {
                                               est_boot <- result.boot$t0[c((length(keepnames)-length(ycol)+1+1+add.parameters):dim(result.boot$t)[2])]}
est_boot2 <- apply(result.boot.ma,2,mean)
se_est_boot <- apply(result.boot.ma,2,sd)
boot.quantiles <- apply(result.boot.ma,2,nquantile)
lower_boot <- boot.quantiles[1,]
upper_boot <- boot.quantiles[2,]
}   


if(criterion=="BIC+"){mynames <- names(av_m[[1]])}
rownames(mycoefflist)   <- mynames 
rownames(myvarlist)     <- mynames 
rownames(mycoefflist.s) <- mynames 
rownames(myvarlist.s)   <- mynames                     

}


#######################################
# Inference after multiple imputation #
#######################################

if(method!="criterion.selection"){
est <- apply(mycoefflist,1,mean)
var_within <- apply(myvarlist^2,1,mean)
coeffdiff<-(matrix(cbind(rep(est,M)),ncol=M,nrow=length(est))-mycoefflist)^2
var_between <-  apply(coeffdiff,1,sum)
variance <- var_within + ((M+1)/(M*(M-1+1e-100)))*var_between
se_est <- round(sqrt(variance),digits=5)
if(M>1){mydf <- (M-1)*(1+((M*var_within)/(((M+1)/(M-1+0.000001))*var_between)))^2}
if(M==1){mydf <- Inf}  # only coarse approx., since df not clear 
mytq <- qt(1-((1-CI)/2), df = mydf)
upper <- est + mytq*se_est
lower <- est - mytq*se_est
}

est.s <- apply(mycoefflist.s,1,mean)
var_within.s <- apply(myvarlist.s^2,1,mean)
coeffdiff.s<-(matrix(cbind(rep(est.s,M)),ncol=M,nrow=length(est.s))-mycoefflist.s)^2
var_between.s <-  apply(coeffdiff.s,1,sum)
variance.s <- var_within.s + ((M+1)/(M*(M-1+1e-100)))*var_between.s
se_est.s <- round(sqrt(variance.s),digits=5)
if(M>1){mydf.s <- (M-1)*(1+((M*var_within.s)/(((M+1)/(M-1+0.000001))*var_between.s)))^2}
if(M==1){mydf.s <- Inf}  # only coarse approx., since df not clear 
mytq.s <- qt(1-((1-CI)/2), df = mydf.s)
upper.s <- est.s + mytq.s*se_est.s
lower.s <- est.s - mytq.s*se_est.s

if(method=="LASSO/LAE" | method=="criterion.average"){
variable_importance <- sort(round(apply(myvarimplist,1,mean),digits=2),decreasing=TRUE)
if(is.null(id)==FALSE & model.family=="coxph"){variable_importance<-variable_importance[((substr(names(variable_importance),1,8)=="cluster(")==FALSE) & ((substr(names(variable_importance),1,7)=="strata(")==FALSE)]}
}

################################################################################

# Final results
mycoefficients=NULL
if(method!="criterion.selection"){
mycoefficients <- round(as.matrix(cbind(est,se_est,lower,upper)),digits=5)
colnames(mycoefficients)<- c("Estimate","Std.Error","Lower CI","Upper CI")}
mycoefficients.b = NULL
if(inference=="+bootstrapping" & method!="criterion.selection"){mycoefficients.b <- round(as.matrix(cbind(est,lower_boot,upper_boot,est_boot2,se_est_boot)),digits=5)}
if(inference=="+bootstrapping"& method!="criterion.selection"){
colnames(mycoefficients.b)<- c("Estimate","Lower CI","Upper CI","(Boot. mean)","(Boot. SE)")
rownames(mycoefficients.b)<-rownames(mycoefficients)
}
mycoefficients.s <- round(as.matrix(cbind(est.s,se_est.s,lower.s,upper.s)),digits=5)
mycoefficients.s[is.na(mycoefficients.s)]<-0
colnames(mycoefficients.s)<- c("Estimate","Std.Error","Lower CI","Upper CI")
mycoefficients.b.s = NULL
if(inference=="+bootstrapping" & method!="MMA" & method!="LASSO/SAE"){mycoefficients.b.s <- round(as.matrix(cbind(est.s,lower_boot.s,upper_boot.s,est_boot2.s,se_est_boot.s)),digits=5)}
if(inference=="+bootstrapping"  & method!="MMA" & method!="LASSO/SAE"){
colnames(mycoefficients.b.s)<- c("Estimate","Lower CI","Upper CI","(Boot. mean)","(Boot. SE)")
rownames(mycoefficients.b.s)<-rownames(mycoefficients.s)
}
if(inference=="+bootstrapping"){colnames(result.boot.s) <- rownames(mycoefficients.s)}
if(inference=="+bootstrapping"){colnames(result.boot.ma) <- rownames(mycoefficients)}
if(model.family=="coxph" & criterion!="BIC+"){
if(method!="criterion.selection"){mycoefficients <- mycoefficients[-1,]}
if(method!="criterion.selection"){mycoefficients.b <- mycoefficients.b[-1,]}
if(method!="criterion.selection"){mycoefficients.s <- mycoefficients.s[-1,]}
if(method!="criterion.selection"){mycoefficients.b.s <- mycoefficients.b.s[-1,]}
}
if(report.exp==TRUE){
if(is.null(mycoefficients)==FALSE){mycoefficients <- round(exp(mycoefficients),digits=4)}
if(is.null(mycoefficients.b)==FALSE){mycoefficients.b <- round(exp(mycoefficients.b),digits=4)}
if(is.null(mycoefficients.s)==FALSE){mycoefficients.s <- round(exp(mycoefficients.s),digits=4)}
if(is.null(mycoefficients.b.s)==FALSE){mycoefficients.b.s <- round(exp(mycoefficients.b.s),digits=4)}
}
if(method=="MMA"){mycoefficients.s=NULL}


################################################################################

ptm3 <-  proc.time()
mytime <- round((((ptm3-ptm)/60))[1],digits=2)
if(print.time==TRUE){cat(paste("The analysis time was", mytime, "minute(s) \n"))}

################################################################################

# Return our results

res= list(coefficients.ma=mycoefficients,
          coefficients.ma.boot = mycoefficients.b,
          coefficients.s = mycoefficients.s,
          coefficients.boot.s = mycoefficients.b.s,
          variable.importance = variable_importance,
          boot.results = list(result.boot.s,result.boot.ma)
)

class(res) <- "mami"
res



##########################

}





