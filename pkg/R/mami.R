###################################################################
# Model selection and model averaging (after multiple imputation) #
###################################################################


mami <- function(X, missing.data=c("imputed","none","CC"),  
  model=c("gaussian","binomial","poisson","cox"), outcome=1, id=NULL, 
  method=c("MA.criterion","MS.criterion","LASSO","LAE","MMA"), 
  criterion=c("AIC","BIC","BIC+","CV","GCV"), kfold=5, cvr=F,
  inference=c("standard","+boot"), B=20, X.org=NULL, CI=0.95,
  var.remove=NULL, aweights=NULL, add.factor=NULL, add.interaction=NULL, 
  add.transformation=NULL, add.stratum=NULL,    
  candidate.models=c("all","restricted","very restricted"), screen=0,
  report.exp=FALSE, print.time=FALSE, print.warnings=TRUE, ...){

ptm <- proc.time()

method            <- match.arg(method)
criterion         <- match.arg(criterion)
inference         <- match.arg(inference)
missing.data      <- match.arg(missing.data)
model      <- match.arg(model)
candidate.models  <- match.arg(candidate.models)

# Checks
if(is.null(aweights)==FALSE & inference=="+boot"){stop("Bootstrapping not allowed when specifying weights.")}
if(is.null(add.interaction)==FALSE & is.list(add.interaction)==FALSE){stop("Use `list()' to add interactions")}
if(is.null(add.stratum)==FALSE & model!="cox"){stop("You can only specify strata in a Cox Model.")}
if(missing.data=="CC" & any(is.na(X))==FALSE & print.warnings==TRUE){stop("You specified a complete case analysis but there is no missing data.")}
if(method=="MS.criterion" & is.null(id)==FALSE & model!="cox" & (criterion!="CV"&criterion!="GCV")){stop("Model selection (with stepAIC) for mixed models not possible.\n However, you can use method=`MA.criterion' instead which also produces results for model selection,\n though definition of AIC in mixed models is not aways clear.")}
if(method=="MS.criterion" & criterion=="BIC" & model=="cox"){stop("BIC not allowed for Cox model. Possibly use `BIC+'")}
if(inference=="+boot" & is.null(X.org) & missing.data=="imputed"){stop("Please also provide unimputed data under `X.org' when using bootstrapping.")}
if(model=="cox" & (criterion=="CV" | criterion=="GCV")){stop("Cross Validation not allowed in the Cox PH model.")}
if(method=="MMA" & model!="gaussian"){stop("Mallows model averaging only works for the linear model.")}
if(criterion=="BIC+"){
if(method=="MS.criterion"){stop("You can use `criterion=BIC' for model selection. `BIC+' is appropriate for model averaging, see manual.")}
if(is.null(id)==FALSE){stop("BIC+ can only be used with cross-sectional data, but you specified an `id'.")}
}
if(method=="LAE"){
if(model=="cox"){stop("Lasso averaging not possible for model='cox'. However, try method='LASSO' instead.")}
if(method=="MA.criterion" & model%in%c("gaussian","binomial","poisson")==FALSE){stop("Only the following families are allowed for LASSO averaging estimation: 'gaussian','binomial','poisson'")}
if((is.null(add.transformation)==FALSE | is.null(add.interaction)==FALSE)){stop("Don't use 'add.transformation' and 'add.interaction' with method='LAE'. \n Prepare your dataset such that transformations and interactions are already contained in it.")}
}
if((criterion=="CV" | criterion=="GCV") & method=="MA.criterion"){stop("Use method='MS.criterion' for model selection with (G)CV.")}
if(criterion=="CV"){
if(kfold==1){stop("Choose kfold>1")}
if(is.null(id)==FALSE & is.null(add.transformation)==FALSE){stop("It is currently not possible to add transformations to mixed models when criterion='CV'")}
}
if((method=="LASSO" | method=="LAE") & is.null(id)==FALSE){stop("LASSO estimation not available for longitudinal data.")}


cores=1
#if(detectCores()<cores){stop(paste("Your maximum amount of available cores is",detectCores(),". Choose cores <= ",detectCores(),"\n",sep=""))}
if(print.warnings==TRUE){
if(method=="MMA"){cat("Note: The results of MMA depend on the ordering of the regressors. \n")}
if(missing.data=="CC"){cat("Note: you should use complete cases only if you are sure that this is appropriate! \n")}
if((method=="MS.criterion" | method=="LASSO") & inference=="standard"){cat("Note: Bootstrap confidence intervals are preferred after model selection. Standard confidence intervals after model selection may be overoptimistic. \n")}
if((method=="LASSO" | method=="LAE")){cat("Note: LASSO CIs based on the bootstrap SE, as reported below, are not bias corrected and not randomization valid. Interpret with caution. \n")} 
if(report.exp==TRUE & model=="gaussian"){warning("Exponentiated coefficients are typically not meaningful for the linear model.")}
if(criterion=="CV"){cat("Note: a much faster option is leave-one-out cross validation approximated by `criterion=GCV'. \n")} 
if(criterion=="CV" & is.null(id)==FALSE){cat("Note that using cross validation in mixed models does not make full use of the random intercept estimates. \n")}
if((is.null(add.transformation)==FALSE | is.null(add.interaction)==FALSE) & (length(add.interaction)+length(add.transformation)>=4)){cat("Note: Many interactions and/or transformation increase model complexity and computation time. \n If the full model can not be fit because of overfitting, mami may crash.\n")}
if(criterion=="BIC+" & (candidate.models=="restricted" | candidate.models=="very restricted")){cat("Note: candidate models are restricted by the leaps algorithm of the `BMA' package. \n")}
if(cores>1 & ((method=="MMA" | method=="LAE" | method=="LASSO") | (method=="MS.criterion" & (criterion!="CV" & criterion!="GCV")))){cat("Note: you benefit from cores>1 only when \n 1) using 'method=MA.criterion' or \n 2) 'method=MS.criterion' and 'criterion=(G)CV' \n")}
}

# Preparation of datamatrix
M = NULL          # Number of imputations
imp.method <- "none"
if(class(X)=="amelia"){
M=length(X$imputations)
X_ <- X$imputations
am <- X$arguments
imp.method <- "amelia"
}
if(class(X)=="data.frame"){
M=1
X_ <- list(X)
}
if(class(X)=="list"){
if(inference=="+boot"){stop("To utilize bootstrapping please impute with Amelia II or mice.")}
M= length(X)
X_ <- X
}
if(class(X)=="mids"){
X_ <- complete(X, action="long")
M =max(as.numeric(X_$.imp))
X_ <- split(X_[,-c(1:2)],X_$.imp)
imp.method <- "mice"
mp <- list(X$method, X$predictorMatrix, X$visitSequence, X$post, X$form)
}
if(is.null(M)){
M=length(X)
X_ <- X
}
if(is.null(M)){stop("Data format not recognized.\n Choose between a data frame; or multiple imputation via Amelia II or mice; or a list of user-specified imputed data sets.")}
X <- X_  
if(is.numeric(outcome[1])==FALSE){outcome <- match(outcome,colnames(X[[1]]))}
if(is.null(add.interaction)==FALSE){add.interaction<-lapply(add.interaction,sort)}

exclude=NULL
if(screen>0){
if(is.null(id)==F){stop("Screening with LASSO not allowed for longitudinal data")}
if(screen>dim(X[[1]])[2]-3-length(outcome)){stop(paste("You are allowed to only screen a maximum of",dim(X[[1]])[2]-2-length(outcome),"variables in this data set"))}
cvIndex <- rep(1:kfold,trunc(dim(X[[1]])[1]/kfold)+1)[1:dim(X[[1]])[1]]
if(cvr==T){cvIndex<-sample(cvIndex)}
vri=NULL
if(is.numeric(var.remove)==FALSE & is.null(var.remove)==FALSE){vri<- match(var.remove,colnames(X[[1]]))}
if(model!="cox"){sm <- cv.glmnet(as.matrix(X[[1]][,-c(outcome,vri)]),as.matrix(X[[1]][,outcome]),foldid=cvIndex,nfolds=kfold,family=model)}
if(model=="cox"){
asi<-add.stratum
if(is.numeric(add.stratum)==FALSE & is.null(add.stratum)==FALSE){asi<- match(add.stratum,colnames(X[[1]]))}
s.surv.outcome <- matrix(cbind(as.numeric(X[[1]][,outcome[1]]),X[[1]][,outcome[2]]),ncol=2,dimnames=list(NULL,c("time","status")))
sm <- cv.glmnet(as.matrix(X[[1]][,-c(outcome,asi,vri)]),s.surv.outcome,foldid=cvIndex,nfolds=kfold,family="cox")
}
si <- abs((max(sm$nzero)-screen)-sm$nzero)
exclude <- rownames(glmnet::coef.glmnet(sm, s = sm$lambda[si==min(si)][1]))[(glmnet::coef.glmnet(sm, s = sm$lambda[si==min(si)][1]))[,1]==0]
var.remove <- c(var.remove,exclude)
if(print.warnings==T){cat(paste("You specified that (about)",screen,"variables should be removed in a screening step.\n", "Therefore, the following variable(s) have been removed (using LASSO):",paste(exclude,collapse=" "),"\n"))}
if(is.null(add.transformation)==F){
sit <- apply(matrix(sapply(exclude, function(exclude) grepl(exclude, add.transformation)),nrow=length(add.transformation)),1,sum)>0
add.transformation <- add.transformation[!sit]
if(print.warnings==T){if(any(sit)){
if(length(add.transformation)==0){add.transformation<-NULL}
cat(paste("Due to screening the list of suggested transformations has been adjusted. \n"))}}
}
if(is.null(add.interaction)==F){
sii <- apply(matrix(sapply(exclude, function(exclude) grepl(exclude, add.interaction)),nrow=length(add.interaction)),1,sum)>0
add.interaction <- add.interaction[!sii]
if(print.warnings==T){if(any(sii)){
if(length(add.interaction)==0){add.interaction<-NULL}
cat(paste("Due to screening the list of suggested interactions has been adjusted. \n"))}}
}
if(is.null(add.factor)==F){
sif <- apply(matrix(sapply(exclude, function(exclude) grepl(exclude, add.factor)),nrow=length(add.factor)),1,sum)>0
add.factor <- add.factor[!sif]
if(print.warnings==T){if(any(sif)){
if(length(add.factor)==0){add.factor<-NULL}
cat(paste("Due to screening the list of suggested factors has been adjusted. \n"))}}
}
}

if(is.numeric(var.remove)==FALSE){var.remove <- match(var.remove,colnames(X[[1]]))}
if(model=="poisson" & is.factor(X[[1]][,outcome])==TRUE){for(m in 1:M){X[[m]][,outcome] <- as.numeric(as.character(X[[m]][,outcome]))}}
if(model=="cox" | model=="poisson"){for(m in 1:M){if(any(sign(X[[m]][,outcome])==-1)){stop("Cox/Poisson model requires non-negative outcome variables")}}}
if(any(!(outcome %in% 1))==TRUE | is.null(var.remove)==FALSE){for(m in 1:M){X[[m]] <- subset(X[[m]],select=c(names(X[[m]])[outcome],names(X[[m]])[-c(outcome,var.remove)])) }}
keepnames <- colnames(X[[1]])
if((candidate.models=="restricted" & (dim(X[[1]])[2]-length(outcome))<5) |  (candidate.models=="very restricted" & (dim(X[[1]])[2]-length(outcome))<10)){stop("You should not restrict your candidate models for this analysis")}
if(missing.data=="CC"){
if(print.warnings==T){cat(paste("Complete case analysis: ",round(dim(na.omit(as.data.frame(X)))[1]/dim(as.data.frame(X))[1],digits=3)*100,"% of data is used. \n", sep=""))}
X <- list(na.omit(as.data.frame(X)))}
if(missing.data!="imputed"){X.org<-X[[1]]}    
if(missing.data=="none"     & any(is.na(X[[1]]))==TRUE){stop("There is still missing data but you specified there is none.")}
if(missing.data=="imputed"  & any(is.na(X[[1]]))==TRUE){stop("There is still missing data but you specified the data is imputed.")}
if(missing.data=="none" & length(X)>1){stop("You specified there is no missing data but you provide multiple datasets as `X'.")}

          
# matrices for storing
est <- se_est <- lower <- upper <- ests_boot <- est_boot <- est_boot2 <- se_est_boot <-lower_boot <-upper_boot <- variable_importance <- NULL 
est.s <- se_est.s <- lower.s <- upper.s <- ests_boot.s <- est_boot.s <- est_boot2.s <- se_est_boot.s <- lower_boot.s <- upper_boot.s  <- result.boot.s <- result.boot.ma  <- NULL


# Regression formula
tf <- function(mymatrix){colnames(mymatrix)[sapply(mymatrix,is.factor)]}
if(length(tf(X[[1]]))>0){
add.factor <- c(add.factor,tf(X[[1]]))
add.factor <- intersect(add.factor,add.factor)
add.factor <- intersect(colnames(X[[1]]),add.factor)
}
save.factor <-add.factor
if(print.warnings==TRUE & is.null(add.factor)==FALSE){cat(paste("Note: The following variables are treated as factors:", paste(add.factor,collapse=" "),"\n"))}
myfamily <- model
oc <- c(1:length(outcome))
myindependent <- keepnames[oc]
usenames <- keepnames
fn <- function(myname){paste("factor(",myname,")",sep="")}
add.parameters <- 0
if(is.null(add.factor)==FALSE){for(i in match(add.factor,usenames)){
add.parameters <- add.parameters + length(levels(factor(X[[1]][,usenames[i]]))) - 2
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
  if(is.na(any(match(add.factor,add.interaction[[j]])))==FALSE & length(na.omit(match(add.factor,add.interaction[[j]])))==1){add.parameters <- add.parameters + length(levels(factor((X[[1]])[,(add.interaction[[j]][na.omit(match(add.factor,add.interaction[[j]]))])])))-1} 
  if(is.na(any(match(add.factor,add.interaction[[j]])))==FALSE & length(na.omit(match(add.factor,add.interaction[[j]])))> 1){for(k in 1:length(na.omit(match(add.factor,add.interaction[[j]])))){add.factor.param <- add.factor.param*(length(levels(factor((X[[1]])[,(add.interaction[[j]][na.omit(match(add.factor,add.interaction[[j]]))[k]])])))-1)} 
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
if(model!="cox"){
usenames <- usenames[-match(id,usenames)]
id <- paste("(1|",id,")",sep="")
add.parameters <- add.parameters-1} else
{usenames <- usenames[-match(id,usenames)]
id <- paste("cluster(",id,")",sep="")
add.parameters <- add.parameters-1} 
}
if(is.null(add.stratum)==FALSE){usenames <- usenames[-match(add.stratum,usenames)]
add.stratum <- paste("strata(",add.stratum,")",sep="")
add.parameters <- add.parameters-1
}
if(is.null(aweights)==FALSE){usenames <- usenames[-match(aweights,usenames)]
add.parameters <- add.parameters-1
} else{aweights<-""} 
mycovariables <- paste(c(usenames[-oc],param.transform,param.interaction,id),collapse="+")
if(length(outcome) > 1){
myindependent <- paste("Surv(",paste(keepnames[oc],collapse=","),")",sep="")
mycovariables <- paste(c(usenames[-oc],param.transform,param.interaction,add.stratum,id),collapse="+")
}
myformula <- paste(myindependent,"~", mycovariables)


# We need these matrices later on...
mycoefflist <- myvarlist <- mycoefflist.s <- myvarlist.s <- matrix(NA,ncol=M,nrow=(length(keepnames)-length(outcome)+1+add.parameters))
if(method!="LASSO" & method!="LAE"){myvarimplist <- matrix(NA,ncol=M,nrow=length(c(usenames[-oc],param.transform,param.interaction)))}else{myvarimplist <- matrix(NA,ncol=M,nrow=(length(keepnames)-length(outcome)+1+add.parameters-1))}
if((criterion=="BIC+" | method=="LASSO") & model=="cox"){
mycoefflist <- myvarlist <- mycoefflist.s <- myvarlist.s <- matrix(NA,ncol=M,nrow=(length(keepnames)-length(outcome)+1+add.parameters-1))
myvarimplist <- matrix(NA,ncol=M,nrow=length(c(usenames[-oc],param.transform,param.interaction))-1)  
}

###########################
# Method: Model selection #
###########################

if(method=="MS.criterion" & (criterion!="CV" & criterion!="GCV")){

ptm1.5 <- proc.time()
if(model=="cox"){mycoefflist.s <- mycoefflist.s[-1,]
                          myvarlist.s <- myvarlist.s[-1,]}        
for(m in 1:M)try({
if(print.warnings==TRUE & missing.data=="imputed"){cat("Busy with model selection in imputed dataset no.",m, "\n")}
if(model!="cox"){full.model <- glm(as.formula(myformula), data = as.data.frame(X[[m]]),   weights=X[[m]]$aweights, family=model)} else{
  if(model=="cox"){full.model <- coxph(as.formula(myformula), data = as.data.frame(X[[m]]),   weights=X[[m]]$aweights)}}
sel_m <- suppressWarnings(stepAIC(full.model, trace=0, direction="both", k=2+(as.numeric(criterion=="BIC")*(-2+log(dim(X[[m]])[1])))))
if(model!="cox"){sel_coef <- summary(sel_m)[[12]]}else{ 
if(model=="cox" & class(sel_m)[1]!="coxph.null"){sel_coef <- summary(sel_m)[[7]]}
if(model=="cox" & class(sel_m)[1]=="coxph.null"){sel_coef<-matrix(rep(0,dim(mycoefflist.s)[1]),ncol=1,dimnames=list(letters[1:dim(mycoefflist.s)[1]],NULL))}}
mycoef   <- matrix(c(rep(0,length(coefficients(full.model)))),ncol=1,nrow=length(coefficients(full.model)),dimnames=list(names(coefficients(full.model)),"coeff."))
mystd    <- matrix(c(rep(0,length(coefficients(full.model)))),ncol=1,nrow=length(coefficients(full.model)),dimnames=list(names(coefficients(full.model)),"coeff."))
for(i in 1:length(mycoef)){for(j in 1:dim(sel_coef)[1]){
if(rownames(mycoef)[i]==rownames(sel_coef)[j]){
mycoef[i] <- sel_coef[j,1]
if(model!="cox"){mystd[i] <- sel_coef[j,2]}else{ 
if(model=="cox" & class(sel_m)[1]!="coxph.null"){mystd[i] <- sel_coef[j,3]}
if(model=="cox" & class(sel_m)[1]=="coxph.null"){mystd[i] <- 0}}
}}}
mycoefflist.s[,m] <- mycoef
myvarlist.s[,m]   <- mystd
})
rownames(mycoefflist.s) <- names(coefficients(full.model)) 
rownames(myvarlist.s) <- names(coefficients(full.model)) 


ptm2 <-  proc.time()
inbetweentime <- round((((ptm2-ptm1.5)/60))[1],digits=2)
if(print.time==TRUE & inference=="+boot"){cat(paste("The analysis will run for about another", inbetweentime*(B+1), "minute(s) \n"))}

if(inference=="+boot"){
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

    if(missing.data=="imputed" & imp.method=="amelia"){mydata.out <-  amelia(mydata, m = M, p2s = as.numeric(print.warnings==TRUE),frontend = FALSE, idvars = am$idvars,
                                  ts = am$ts, cs = am$cs, polytime = am$polytime, splinetime = am$splinetime, intercs = am$intercs,
                                  lags = am$lags, leads = am$leads, startvals = am$startvals, tolerance = am$tolerance,
                                  logs = am$logs, sqrts = am$sqrts, lgstc = am$lgstc, noms = am$noms, ords = am$ords,
                                  incheck = TRUE, collect = FALSE, arglist = am$arglist, empri = am$empri,
                                  priors = am$priors, autopri = am$autopri, emburn = c(0,0), bounds = am$bounds,
                                  max.resample = am$max.resample, overimp = am$overimp)}
    if(missing.data=="imputed" & imp.method=="mice"){mydata.out <-  mice(mydata, m = M, method=mp[[1]], predictorMatrix=mp[[2]],
                                  visitSequence=mp[[3]], post= mp[[4]], printFlag=print.warnings, form=mp[[5]])
                                  mydata.out <- complete(mydata.out, action="long")
                                  mydata.out <- split(mydata.out[,-c(1:2)],mydata.out$.imp)
                                  mydata.out$imputations <- mydata.out
                                  }  
    if(missing.data!="imputed"){mydata.out <- NULL
                                mydata.out$imputations[[1]] <- mydata} 
if(any(!(outcome %in% 1))==TRUE | is.null(var.remove)==FALSE){for(m in 1:M){mydata.out$imputations[[m]] <- try(subset(mydata.out$imputations[[m]],select=c(names(mydata.out$imputations[[m]])[outcome],names(mydata.out$imputations[[m]])[-c(outcome,var.remove)])))}}
mycoefflist_b <- matrix(NA,ncol=M,nrow=(length(keepnames)-length(outcome)+1+add.parameters))
if(model=="cox"){mycoefflist_b <- mycoefflist_b[-1,]} 
for(m in 1:M)try({ 
if(print.warnings==TRUE & missing.data=="imputed"){cat("Busy with model selection in imputed dataset no.",m, "\n")}
if(model!="cox"){full.model.b <- glm(as.formula(myformula), data = as.data.frame(mydata.out$imputations[[m]]), family=model)} else{
  if(model=="cox"){full.model.b <- coxph(as.formula(myformula), data = as.data.frame(mydata.out$imputations[[m]]))} }
sel_m.b <- try(suppressWarnings(stepAIC(full.model.b, trace=0, direction="both", k=2+(as.numeric(criterion=="BIC")*(-2+log(dim(X[[m]])[1]))))))
if(model!="cox"){sel_coef.b <- summary(sel_m.b)[[12]]} else{
   if(model=="cox" & class(sel_m)[1]!="coxph.null"){sel_coef.b <- summary(sel_m.b)[[7]]}
   if(model=="cox" & class(sel_m)[1]=="coxph.null"){sel_coef.b <- matrix(rep(0,dim(mycoefflist_b)[1]),ncol=1,dimnames=list(letters[1:dim(mycoefflist_b)[1]],NULL))}
   }
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
if(print.time==TRUE & inference=="+boot"){cat(paste("The analysis will run for about another", inbetweentime*(B+1), "minute(s) \n"))}

if(inference=="+boot"){
if(print.warnings==TRUE){cat("Now start with bootstrap estimation... \n" )}
mami.boot <- function(mydata,indices){
mydata <- mydata[indices,]
  if(missing.data=="imputed" & imp.method=="amelia"){mydata.out <-  amelia(mydata, m = M, p2s = as.numeric(print.warnings==TRUE),frontend = FALSE, idvars = am$idvars,
                                  ts = am$ts, cs = am$cs, polytime = am$polytime, splinetime = am$splinetime, intercs = am$intercs,
                                  lags = am$lags, leads = am$leads, startvals = am$startvals, tolerance = am$tolerance,
                                  logs = am$logs, sqrts = am$sqrts, lgstc = am$lgstc, noms = am$noms, ords = am$ords,
                                  incheck = TRUE, collect = FALSE, arglist = am$arglist, empri = am$empri,
                                  priors = am$priors, autopri = am$autopri, emburn = c(0,0), bounds = am$bounds,
                                  max.resample = am$max.resample, overimp = am$overimp)} 
  if(missing.data=="imputed" & imp.method=="mice"){mydata.out <-  mice(mydata, m = M, method=mp[[1]], predictorMatrix=mp[[2]],
                                  visitSequence=mp[[3]], post= mp[[4]], printFlag=print.warnings, form=mp[[5]])
                                  mydata.out <- complete(mydata.out, action="long")
                                  mydata.out <- split(mydata.out[,-c(1:2)],mydata.out$.imp)
                                  mydata.out$imputations <- mydata.out
                                  } 
  if(missing.data!="imputed"){mydata.out <- NULL
                              mydata.out$imputations[[1]] <- mydata} 
if(any(!(outcome %in% 1))==TRUE | is.null(var.remove)==FALSE){for(m in 1:M){mydata.out$imputations[[m]] <- try(subset(mydata.out$imputations[[m]],select=c(names(mydata.out$imputations[[m]])[outcome],names(mydata.out$imputations[[m]])[-c(outcome,var.remove)])))}}
mycoefflist_b <- matrix(NA,ncol=M,nrow=(length(keepnames)-length(outcome)+1+add.parameters))
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
# Method: LAE   #
#################

if(method=="LAE"){

ptm1.5 <- proc.time()

for(m in 1:M){
if(print.warnings==TRUE & missing.data=="imputed"){cat("Busy with LASSO (averaging) estimation in imputed dataset no.",m, "\n")}
av_m <- lae(X[[m]],ycol=1,kfold=kfold, factor.variables=save.factor, calc.variance=TRUE, glm.family=model,pd=F,random=cvr,...)
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
if(print.time==TRUE & inference=="+boot"){cat(paste("The analysis will run for about another", inbetweentime*(B+1), "minute(s) \n"))}

if(inference=="+boot"){
if(print.warnings==TRUE){cat("Now start with bootstrap estimation... \n" )}
mami.boot <- function(mydata,indices){
mydata <- mydata[indices,]
  if(missing.data=="imputed" & imp.method=="amelia"){mydata.out <-  amelia(mydata, m = M, p2s = as.numeric(print.warnings==TRUE),frontend = FALSE, idvars = am$idvars,
                                  ts = am$ts, cs = am$cs, polytime = am$polytime, splinetime = am$splinetime, intercs = am$intercs,
                                  lags = am$lags, leads = am$leads, startvals = am$startvals, tolerance = am$tolerance,
                                  logs = am$logs, sqrts = am$sqrts, lgstc = am$lgstc, noms = am$noms, ords = am$ords,
                                  incheck = TRUE, collect = FALSE, arglist = am$arglist, empri = am$empri,
                                  priors = am$priors, autopri = am$autopri, emburn = c(0,0), bounds = am$bounds,
                                  max.resample = am$max.resample, overimp = am$overimp)} 
  if(missing.data=="imputed" & imp.method=="mice"){mydata.out <-  mice(mydata, m = M, method=mp[[1]], predictorMatrix=mp[[2]],
                                  visitSequence=mp[[3]], post= mp[[4]], printFlag=print.warnings, form=mp[[5]])
                                  mydata.out <- complete(mydata.out, action="long")
                                  mydata.out <- split(mydata.out[,-c(1:2)],mydata.out$.imp)
                                  mydata.out$imputations <- mydata.out
                                  } 
  if(missing.data!="imputed"){mydata.out <- NULL
                              mydata.out$imputations[[1]] <- mydata} 
if(any(!(outcome %in% 1))==TRUE | is.null(var.remove)==FALSE){for(m in 1:M){mydata.out$imputations[[m]] <- try(subset(mydata.out$imputations[[m]],select=c(names(mydata.out$imputations[[m]])[outcome],names(mydata.out$imputations[[m]])[-c(outcome,var.remove)])))}}
mycoefflist_b <- matrix(NA,ncol=M,nrow=(length(keepnames)-length(outcome)+1+add.parameters))
mycoefflist_b.s <- matrix(NA,ncol=M,nrow=(length(keepnames)-length(outcome)+1+add.parameters))
for(m in 1:M)try({ 
if(print.warnings==TRUE & missing.data=="imputed"){cat("Busy with LASSO (averaging) estimation in imputed dataset no.",m, "\n")} 
av_m.b <- lae(mydata.out$imputations[[m]],ycol=1,kfold=kfold, factor.variables=NULL, calc.variance=TRUE, glm.family=model,pd=F,random=cvr,...) 
if(length(mycoefflist_b[,m]) == length(av_m.b$coefficients[,1])){
mycoefflist_b[,m] <- av_m.b$coefficients[,1]
mycoefflist_b.s[,m] <- av_m.b$coefficients[,3]
}})
c(apply(mycoefflist_b.s,1,mean),apply(mycoefflist_b,1,mean))
}
result.boot <- try(boot(X.org,mami.boot,B))
if(any(apply(apply(result.boot$t,1,is.na),2,any))==TRUE & print.warnings==TRUE){cat("There occurred fitting problems in one or many bootstrap samples. The relevant samples have been removed. \n")}
result.boot.s <- result.boot$t[apply(apply(result.boot$t[,1:(length(keepnames)-length(outcome)+1+add.parameters)],1,is.na),2,any)==F,1:(length(keepnames)-length(outcome)+1+add.parameters)]    
ests_boot.s <- result.boot.s
est_boot.s <- result.boot$t0[1:(length(keepnames)-length(outcome)+1+add.parameters)]
est_boot2.s <- apply(result.boot.s,2,mean)
se_est_boot.s <- apply(result.boot.s,2,sd)
nquantile <- function(myvector){quantile(myvector,probs=c(((1-CI)/2),1-((1-CI)/2)))}
boot.quantiles.s <- apply(result.boot.s,2,nquantile)
lower_boot.s <- boot.quantiles.s[1,]
upper_boot.s <- boot.quantiles.s[2,]
result.boot.ma <- result.boot$t[apply(apply(result.boot$t[,c((length(keepnames)-length(outcome)+1+1+add.parameters):dim(result.boot$t)[2])],1,is.na),2,any)==F,c((length(keepnames)-length(outcome)+1+1+add.parameters):dim(result.boot$t)[2])]
ests_boot <- result.boot.ma
est_boot <- result.boot$t0[c((length(keepnames)-length(outcome)+1+1+add.parameters):dim(result.boot$t)[2])]
est_boot2 <- apply(result.boot.ma,2,mean)
se_est_boot <- apply(result.boot.ma,2,sd)
boot.quantiles <- apply(result.boot.ma,2,nquantile)
lower_boot <- boot.quantiles[1,]
upper_boot <- boot.quantiles[2,]
} 
}

#################
# Method: LASSO #
#################

if(method=="LASSO"){

ptm1.5 <- proc.time()           

cvIndex <- rep(1:kfold,trunc(dim(X[[1]])[1]/kfold)+1)[1:dim(X[[1]])[1]]
if(cvr==T){cvIndex<-sample(cvIndex)}
keepnamesL <- colnames(X[[1]])
if(model!="cox"){
boot.lasso <- function(mydata,indices){
mydata <- mydata[indices,]
myboot <- cv.glmnet(mydata[,-1],mydata[,1],foldid=cvIndex,nfolds=kfold,family=model,...)
as.vector(glmnet::coef.glmnet(myboot, s = myboot$lambda.min))
}}else{
boot.cox.lasso <- function(mydata,indices){
mydata <- mydata[indices,]
myboot <- cv.glmnet(mydata[,-c(1:2)],mydata[,1:2],foldid=cvIndex,nfolds=kfold,family="cox",...)
as.vector(glmnet::coef.glmnet(myboot, s = myboot$lambda.min))
}
}

for(m in 1:M){
if(print.warnings==TRUE & m==1){cat(paste("Note: The LASSO standard error is calculated based on", max(100,B), "bootstrap samples. \n"))}
if(print.warnings==TRUE & missing.data=="imputed"){cat("Busy with LASSO model selection in imputed dataset no.",m, "\n")}
if(model=="cox"){surv.outcome <- matrix(cbind(as.numeric(X[[m]][,1]),X[[m]][,2]),ncol=2,dimnames=list(NULL,c("time","status")))}
X[[m]] <- cbind(X[[m]][,1],model.matrix(as.formula(myformula),data=as.data.frame(X[[m]]))[,-1])
colnames(X[[m]])[1]<-keepnamesL[1]
if(model!="cox"){lasso  <- cv.glmnet(X[[m]][,-1],X[[m]][,1],foldid=cvIndex,nfolds=kfold,family=model,...)}
if(model=="cox"){lasso  <- cv.glmnet(as.matrix(X[[m]][,-1]),surv.outcome,foldid=cvIndex,nfolds=kfold,family="cox",...)}
lasso_coeff <- glmnet::coef.glmnet(lasso, s = c(lasso$lambda.min))
if(model!="cox"){mylasso_boot<- boot(X[[m]],boot.lasso,max(100,B))}
if(model=="cox"){mylasso_boot<- boot(cbind(surv.outcome,as.matrix(X[[m]][,-1])),boot.cox.lasso,max(100,B))}
var_lasso <-matrix(apply(mylasso_boot$t,2,sd),ncol=1,nrow=nrow(lasso_coeff))
mycoefflist.s[,m] <-  as.vector(lasso_coeff)
myvarlist.s[,m]   <-  var_lasso 
}
if(model!="cox"){
rownames(mycoefflist.s) <- c("Intercept",colnames(X[[m]])[-1])
rownames(myvarlist.s)     <- c("Intercept",colnames(X[[m]])[-1])}else{
rownames(mycoefflist.s) <- c(colnames(X[[m]])[-1])
rownames(myvarlist.s)     <- c(colnames(X[[m]])[-1])
}


ptm2 <-  proc.time()
inbetweentime <- round((((ptm2-ptm1.5)/60))[1],digits=2)
if(print.time==TRUE & inference=="+boot"){cat(paste("The analysis will run for about another", inbetweentime*(B+1), "minute(s) \n"))}

if(inference=="+boot"){
if(print.warnings==TRUE){cat("Now start with bootstrap estimation... \n" )}
mami.boot <- function(mydata,indices){
mydata <- mydata[indices,]
  if(missing.data=="imputed" & imp.method=="amelia"){mydata.out <-  amelia(mydata, m = M, p2s = as.numeric(print.warnings==TRUE),frontend = FALSE, idvars = am$idvars,
                                  ts = am$ts, cs = am$cs, polytime = am$polytime, splinetime = am$splinetime, intercs = am$intercs,
                                  lags = am$lags, leads = am$leads, startvals = am$startvals, tolerance = am$tolerance,
                                  logs = am$logs, sqrts = am$sqrts, lgstc = am$lgstc, noms = am$noms, ords = am$ords,
                                  incheck = TRUE, collect = FALSE, arglist = am$arglist, empri = am$empri,
                                  priors = am$priors, autopri = am$autopri, emburn = c(0,0), bounds = am$bounds,
                                  max.resample = am$max.resample, overimp = am$overimp)} 
  if(missing.data=="imputed" & imp.method=="mice"){mydata.out <-  mice(mydata, m = M, method=mp[[1]], predictorMatrix=mp[[2]],
                                  visitSequence=mp[[3]], post= mp[[4]], printFlag=print.warnings, form=mp[[5]])
                                  mydata.out <- complete(mydata.out, action="long")
                                  mydata.out <- split(mydata.out[,-c(1:2)],mydata.out$.imp)
                                  mydata.out$imputations <- mydata.out
                                  } 
  if(missing.data!="imputed"){mydata.out <- NULL
                              mydata.out$imputations[[1]] <- mydata} 
if(any(!(outcome %in% 1))==TRUE | is.null(var.remove)==FALSE){for(m in 1:M){mydata.out$imputations[[m]] <- try(subset(mydata.out$imputations[[m]],select=c(names(mydata.out$imputations[[m]])[outcome],names(mydata.out$imputations[[m]])[-c(outcome,var.remove)])))}}
if(model!="cox"){mycoefflist_b.s <- matrix(NA,ncol=M,nrow=(length(keepnames)-length(outcome)+1+add.parameters))}else{
mycoefflist_b.s <- matrix(NA,ncol=M,nrow=(length(keepnames)-length(outcome)+1+add.parameters-1))}
for(m in 1:M)try({ 
if(print.warnings==TRUE & missing.data=="imputed"){cat("Busy with LASSO model selection in imputed dataset no.",m, "\n")} 
if(model=="cox"){surv.outcome.b <- matrix(cbind(as.numeric(mydata.out$imputations[[m]][,1]),mydata.out$imputations[[m]][,2]),ncol=2,dimnames=list(NULL,c("time","status")))}
mydata.out$imputations[[m]] <- cbind(mydata.out$imputations[[m]][,1],model.matrix(as.formula(paste("~", mycovariables)),data=as.data.frame(mydata.out$imputations[[m]][,-1]))[,-1])
colnames(mydata.out$imputations[[m]])[1]<-keepnamesL[1]
if(model!="cox"){lasso.b  <- cv.glmnet(mydata.out$imputations[[m]][,-1],mydata.out$imputations[[m]][,1],foldid=cvIndex,nfolds=kfold,family=model,...)}
if(model=="cox"){lasso.b  <- cv.glmnet(as.matrix(mydata.out$imputations[[m]][,-1]),surv.outcome.b,foldid=cvIndex,nfolds=kfold,family="cox",...)}
lasso_coeff.b <- glmnet::coef.glmnet(lasso.b, s = c(lasso.b$lambda.min))
if(model!="cox"){mylasso_boot.b<- boot(mydata.out$imputations[[m]],boot.lasso,max(100,B))}
if(model=="cox"){mylasso_boot.b<- boot(cbind(surv.outcome.b,as.matrix(mydata.out$imputations[[m]][,-1])),boot.cox.lasso,max(100,B))}
var_lasso.b <-matrix(apply(mylasso_boot.b$t,2,sd),ncol=1,nrow=nrow(lasso_coeff))
mycoefflist_b.s[,m] <-  as.vector(lasso_coeff.b)
})
c(apply(mycoefflist_b.s,1,mean))
}
result.boot <- try(boot(X.org,mami.boot,B))
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


#############################################
# Method:  Model averaging, Akaike weights  #
#############################################


if(method=="MA.criterion" | (method=="MS.criterion" & (criterion=="CV"|criterion=="GCV"))){

if(cores==1){clust=NULL}else{
if(print.warnings==TRUE){cat(paste("Note: You initialized parallel computation using",cores,"cores...initializing cluster now..."))}
clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
clust <- try(makeCluster(getOption("cl.cores", cores), type = clusterType))
if(print.warnings==TRUE){
print(clust)
cat("\n")
}}

mymax <- round((dim(mycoefflist)[1]-1*((criterion=="BIC+" & model=="cox")==FALSE))/(as.numeric(candidate.models==c("all"))+ 2*as.numeric(candidate.models==c("restricted"))+4*as.numeric(candidate.models==c("very restricted"))), digits=0)
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

if(model=="cox"){
ma.coxph <- function(maobject){
ctb <- summary(maobject)[[9]]
estimate <- c(0,ctb[,1][c(1,order(rownames(ctb))[-1])])
variance <- c(0,ctb[,2][c(1,order(rownames(ctb))[-1])])
names(estimate)[1] <- names(variance)[1] <- "(Intercept)" 
importance <- (maobject$importance[substr(names(maobject$importance),1,8)!="cluster(" & substr(names(maobject$importance),1,7)!="strata("])#[1:(length(estimate)-1)]
importance <- c(importance[order(names(importance))])
return(list(estimate,variance,importance))
}
}

if(criterion=="CV" & is.null(id)==TRUE){
CV <- function(mymodel,kf=kfold){
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
}
if(criterion=="CV" & is.null(id)==FALSE){
CV <- function(mymodel,kf=kfold){
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
}



ptm1.5 <- proc.time()
myil <- rep(list(NA),M)

for(m in 1:M)try({
ddata <- as.data.frame(X[[m]]) 
if(print.warnings==TRUE & missing.data=="imputed" & (criterion!="CV" & criterion!="GCV")){cat("Busy with model averaging in imputed dataset no.",m, "\n")} 
if(print.warnings==TRUE & missing.data=="imputed" & (criterion=="CV" | criterion=="GCV")){cat("Busy with model selection in imputed dataset no.",m, "\n")}
if(criterion!="BIC+"){                           
if(model!="cox" & is.null(id)){mymodel <- glm(as.formula(myformula), data = ddata, na.action=na.fail, weights=X[[m]]$aweights, family=model)} 
if(model!="cox" & is.null(id)==FALSE){mymodel <- suppressWarnings(glmer(as.formula(myformula), na.action=na.fail, data = ddata, weights=X[[m]]$aweights, family=model))}  # suppress warnings since we use glmer for gaussian family
if(model=="cox"){mymodel <- coxph(as.formula(myformula), data = ddata, na.action=na.fail, weights=X[[m]]$aweights)}
if(cores>1){clusterExport(clust, list("model", "X", "m"), envir=environment())}
mymo <- suppressMessages(suppressWarnings(pdredge(mymodel, m.lim=c((length(c(id,add.stratum))+2)*(model=="cox"),mymax+length(c(id,add.stratum))+1*(model=="cox")),rank=get(criterion),fixed=c(id,add.stratum), cluster=clust,...)))
if(print.warnings==TRUE & (criterion!="CV" & criterion!="GCV") & m==M){cat(paste("Model averaging was based on",dim(mymo)[1],"candidate models \n"))}
if(print.warnings==TRUE & (criterion=="CV" | criterion=="GCV") & m==M){cat(paste("Model selection was based on",dim(mymo)[1],"candidate models \n"))}
if(any(is.na(mymo[,"delta"]))==FALSE | model=="cox"){mym  <- suppressWarnings(get.models(mymo,subset=NA))} else {stop("Problems in model averaging function `dredge'. Possibly try `criterion='BIC+' '.")}
if(model!="cox"){av_m <- suppressWarnings(ma(model.avg(mym,subset=NA)))} else{
       if(model=="cox"){av_m <- try(suppressWarnings(ma.coxph(model.avg(mym,subset=NA))))}}
                        sel_m <- ms(coefficients(summary(mym[[1]])),names(av_m[[1]]),model)
  } else {
  if(criterion=="BIC+" & model!="cox"){mymodel   <- bic.glm(as.formula(myformula), data = ddata, na.action=na.fail, wt=X[[m]]$aweights, glm.family=model,...)}
  if(criterion=="BIC+" & model=="cox"){mymodel   <- suppressWarnings(bic.surv(as.formula(myformula), data = ddata, na.action=na.fail, wt=X[[m]]$aweights,...))}
  bicma <- ma.bma(mymodel)
  av_m  <- bicma[1:3]
  sel_m <- bicma[4:5]
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
if(print.time==TRUE & inference=="+boot"){cat(paste("The analysis will run for about another", inbetweentime*(B+1), "minute(s) \n"))}
if(!(model!="cox" & is.null(id)==FALSE)){mynames <- names(av_m[[1]])}else{
mynames<- names(av_m[[1]])#rownames(summary(mymodel)[[10]])
}
if((length(keepnames)-length(oc)+1+add.parameters)!=length(mynames)){mynames<-c("(Intercept)",mynames)}

if(inference=="+boot"){
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
  if(missing.data=="imputed" & imp.method=="amelia"){mydata.out <-  amelia(mydata, m = M, p2s = as.numeric(print.warnings==TRUE),frontend = FALSE, idvars = am$idvars,
                                  ts = am$ts, cs = am$cs, polytime = am$polytime, splinetime = am$splinetime, intercs = am$intercs,
                                  lags = am$lags, leads = am$leads, startvals = am$startvals, tolerance = am$tolerance,
                                  logs = am$logs, sqrts = am$sqrts, lgstc = am$lgstc, noms = am$noms, ords = am$ords,
                                  incheck = TRUE, collect = FALSE, arglist = am$arglist, empri = am$empri,
                                  priors = am$priors, autopri = am$autopri, emburn = c(0,0), bounds = am$bounds,
                                  max.resample = am$max.resample, overimp = am$overimp, boot.type = am$boot.type)} 
  if(missing.data=="imputed" & imp.method=="mice"){mydata.out <-  mice(mydata, m = M, method=mp[[1]], predictorMatrix=mp[[2]],
                                  visitSequence=mp[[3]], post= mp[[4]], printFlag=print.warnings, form=mp[[5]])
                                  mydata.out <- complete(mydata.out, action="long")
                                  mydata.out <- split(mydata.out[,-c(1:2)],mydata.out$.imp)
                                  mydata.out$imputations <- mydata.out
                                  } 
  if(missing.data!="imputed"){mydata.out <- NULL
                              mydata.out$imputations[[1]] <- mydata} 
if(any(!(outcome %in% 1))==TRUE | is.null(var.remove)==FALSE){for(m in 1:M){mydata.out$imputations[[m]] <- try(subset(mydata.out$imputations[[m]],select=c(names(mydata.out$imputations[[m]])[outcome],names(mydata.out$imputations[[m]])[-c(outcome,var.remove)])))}}
if(criterion!="BIC+"){
mycoefflist_b   <- matrix(NA,ncol=M,nrow=(length(keepnames)-length(outcome)+1+add.parameters))
mycoefflist_b.s <- matrix(NA,ncol=M,nrow=(length(keepnames)-length(outcome)+1+add.parameters))} else {
if(model!="cox"){mycoefflist_b      <- matrix(NA,ncol=M,nrow=(length(keepnames)-length(outcome)+1+add.parameters))
                          mycoefflist_b.s    <- matrix(NA,ncol=M,nrow=(length(keepnames)-length(outcome)+1+add.parameters)) }
if(model=="cox"){mycoefflist_b      <- matrix(NA,ncol=M,nrow=(length(keepnames)-length(outcome)  +add.parameters))
                          mycoefflist_b.s    <- matrix(NA,ncol=M,nrow=(length(keepnames)-length(outcome)  +add.parameters)) }
}
for(m in 1:M)try({  # try() is used here to let the bootstrap run even when a single bootstrap sample causes the specified model to crash. Will later be removed.
bddata <- as.data.frame(mydata.out$imputations[[m]])
if(print.warnings==TRUE & missing.data=="imputed" & (criterion!="CV" & criterion!="GCV")){cat("Busy with model averaging in imputed dataset no.",m, "\n")} 
if(print.warnings==TRUE & missing.data=="imputed" & (criterion=="CV" | criterion=="GCV")){cat("Busy with model selection in imputed dataset no.",m, "\n")}
if(criterion!="BIC+"){
if(model!="cox" & is.null(id)){       mymodel.b <- try(glm(as.formula(myformula), na.action=na.fail, data = bddata, family=model))} 
if(model!="cox" & is.null(id)==FALSE){mymodel.b <- suppressWarnings(glmer(as.formula(myformula), na.action=na.fail, data = bddata, family=model))}  # suppress warnings since we use glmer for gaussian family
if(model=="cox"){                     mymodel.b <- try(coxph(as.formula(myformula), na.action=na.fail, data = bddata))}
if(cores>1){clusterExport(clust, c("model", "X", "m"))}
mym.b <- suppressMessages(suppressWarnings(get.models(pdredge(mymodel.b, m.lim=c((length(c(id,add.stratum))+2)*(model=="cox"),mymax+length(c(id,add.stratum))+1*(model=="cox")),rank=get(criterion),fixed=c(id,add.stratum),cluster=clust,...),subset=NA)))
if(model!="cox"){av_m.b <- try(suppressWarnings(ma(model.avg(mym.b,subset=NA))))} else {
       if(model=="cox"){av_m.b <- try(suppressWarnings(ma.coxph(model.avg(mym.b,subset=NA))))}}
                        sel_mb <- try(ms(coefficients(summary(mym.b[[1]])),mynames,model))
  } else {
  if(criterion=="BIC+" & model!="cox"){mymodel.b   <- try(bic.glm(as.formula(myformula), data = as.data.frame(mydata.out$imputations[[m]]), na.action=na.fail, glm.family=model,...)) }
  if(criterion=="BIC+" & model=="cox"){mymodel.b   <- try(bic.surv(as.formula(myformula), data = as.data.frame(mydata.out$imputations[[m]]), na.action=na.fail,...))}
  bicma.b <- try(ma.bma(mymodel.b))
  av_m.b  <- bicma.b[1:3]
  sel_mb <- bicma.b[4:5]
  }
if(length(mycoefflist[,m])==length(av_m[[1]])){mycoefflist[,m] <- av_m[[1]]}
if(length(mycoefflist_b[,m])==length(av_m.b[[1]])){mycoefflist_b[,m] <- av_m.b[[1]]}
if(length(mycoefflist_b.s[,m])==length(sel_mb[[1]])){mycoefflist_b.s[,m] <- sel_mb[[1]]}
})
c(apply(mycoefflist_b.s,1,mean),apply(mycoefflist_b,1,mean))
}
result.boot <- try(boot(X.org,mami.boot,B))
if(any(apply(apply(result.boot$t,1,is.na),2,any))==TRUE & print.warnings==TRUE & method!="MS.criterion"){cat("There occurred fitting problems in one or many bootstrap samples. The relevant samples have been removed. \n")}
if(criterion=="BIC+" & model=="cox"){ result.boot.s <- result.boot$t[apply(apply(result.boot$t[,1:(length(keepnames)-length(outcome)  +add.parameters)],1,is.na),2,any)==F,1:(length(keepnames)-length(outcome)  +add.parameters)] } else {  
                                               result.boot.s <- result.boot$t[apply(apply(result.boot$t[,1:(length(keepnames)-length(outcome)+1+add.parameters)],1,is.na),2,any)==F,1:(length(keepnames)-length(outcome)+1+add.parameters)]    }
ests_boot.s <- result.boot.s
if(criterion=="BIC+" & model=="cox"){ est_boot.s <- result.boot$t0[1:(length(keepnames)-length(outcome)  +add.parameters)]} else {
                                               est_boot.s <- result.boot$t0[1:(length(keepnames)-length(outcome)+1+add.parameters)]}
est_boot2.s <- apply(result.boot.s,2,mean)
se_est_boot.s <- apply(result.boot.s,2,sd)
nquantile <- function(myvector){quantile(myvector,probs=c(((1-CI)/2),1-((1-CI)/2)))}
boot.quantiles.s <- apply(result.boot.s,2,nquantile)
lower_boot.s <- boot.quantiles.s[1,]
upper_boot.s <- boot.quantiles.s[2,]
if(criterion=="BIC+" & model=="cox"){result.boot.ma <- result.boot$t[apply(apply(result.boot$t[,c((length(keepnames)-length(outcome)+1  +add.parameters):dim(result.boot$t)[2])],1,is.na),2,any)==F,c((length(keepnames)-length(outcome)  +1+add.parameters):dim(result.boot$t)[2])]} else {
                                              result.boot.ma <- result.boot$t[apply(apply(result.boot$t[,c((length(keepnames)-length(outcome)+1+1+add.parameters):dim(result.boot$t)[2])],1,is.na),2,any)==F,c((length(keepnames)-length(outcome)+1+1+add.parameters):dim(result.boot$t)[2])]}
ests_boot <- result.boot.ma
if(criterion=="BIC+" & model=="cox"){ est_boot <- result.boot$t0[c((length(keepnames)-length(outcome)+1  +add.parameters):dim(result.boot$t)[2])]} else {
                                               est_boot <- result.boot$t0[c((length(keepnames)-length(outcome)+1+1+add.parameters):dim(result.boot$t)[2])]}
est_boot2 <- apply(result.boot.ma,2,mean)
se_est_boot <- apply(result.boot.ma,2,sd)
boot.quantiles <- apply(result.boot.ma,2,nquantile)
lower_boot <- boot.quantiles[1,]
upper_boot <- boot.quantiles[2,]
}   

if(cores>1){stopCluster(clust)}

if(criterion!="BIC+"){
rownames(mycoefflist)   <- names(av_m[[1]]) 
rownames(myvarlist)     <- names(av_m[[2]]) 
rownames(mycoefflist.s) <- colnames(sel_m[[1]]) 
rownames(myvarlist.s)   <- colnames(sel_m[[2]])} else{ 
mynames <- make.names(names(av_m[[1]]),unique=TRUE) 
rownames(mycoefflist)   <- mynames
rownames(myvarlist)     <- mynames
rownames(mycoefflist.s) <- mynames
rownames(myvarlist.s)   <- mynames               
}

}

#######################################
# Inference after multiple imputation #
#######################################

if(method!="MS.criterion" & method!="LASSO"){
est <- apply(mycoefflist,1,mean)
var_within <- apply(myvarlist^2,1,mean)
coeffdiff<-(matrix(cbind(rep(est,M)),ncol=M,nrow=length(est))-mycoefflist)^2
var_between <-  apply(coeffdiff,1,sum)
variance <- var_within + ((M+1)/(M*(M-1+1e-100)))*var_between
se_est <- round(sqrt(variance),digits=5)
if(M>1){mydf <- (M-1)*(1+((M*var_within)/(((M+1)/(M-1+0.000001))*var_between)))^2}
if(M==1){mydf <- Inf}  
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
if(M>1){mydf.s <- (M-1)*(1+((M*var_within.s)/((M+1)*var_between.s)))^2}
if(M==1){mydf.s <- Inf}  
mytq.s <- qt(1-((1-CI)/2), df = mydf.s)
upper.s <- est.s + mytq.s*se_est.s
lower.s <- est.s - mytq.s*se_est.s

if(method=="LAE" | method=="MA.criterion"){
variable_importance <- sort(round(apply(myvarimplist,1,mean),digits=2),decreasing=TRUE)
if(is.null(id)==FALSE & model=="cox"){variable_importance<-variable_importance[((substr(names(variable_importance),1,8)=="cluster(")==FALSE) & ((substr(names(variable_importance),1,7)=="strata(")==FALSE)]}
}

################################################################################

# Final results
mycoefficients=NULL
if(method!="MS.criterion" & method!="LASSO"){
mycoefficients <- round(as.matrix(cbind(est,se_est,lower,upper)),digits=6)
colnames(mycoefficients)<- c("Estimate","Std.Error","Lower CI","Upper CI")}
mycoefficients.b = NULL
if(inference=="+boot" & (method!="MS.criterion" & method!="LASSO")){mycoefficients.b <- round(as.matrix(cbind(est,lower_boot,upper_boot,est_boot2,se_est_boot)),digits=6)}
if(inference=="+boot" & (method!="MS.criterion" & method!="LASSO")){
colnames(mycoefficients.b)<- c("Estimate","Lower CI","Upper CI","(Boot. mean)","(Boot. SE)")
rownames(mycoefficients.b)<-rownames(mycoefficients)
}
mycoefficients.s <- round(as.matrix(cbind(est.s,se_est.s,lower.s,upper.s)),digits=6)
mycoefficients.s[is.na(mycoefficients.s)]<-0
colnames(mycoefficients.s)<- c("Estimate","Std.Error","Lower CI","Upper CI")
mycoefficients.b.s = NULL
if(inference=="+boot" & method!="MMA"){mycoefficients.b.s <- round(as.matrix(cbind(est.s,lower_boot.s,upper_boot.s,est_boot2.s,se_est_boot.s)),digits=6)}
if(inference=="+boot"  & method!="MMA"){
colnames(mycoefficients.b.s)<- c("Estimate","Lower CI","Upper CI","(Boot. mean)","(Boot. SE)")
rownames(mycoefficients.b.s)<-rownames(mycoefficients.s)
}
if(inference=="+boot"){colnames(result.boot.s) <- rownames(mycoefficients.s)}
if(inference=="+boot"){colnames(result.boot.ma) <- rownames(mycoefficients)}
if(model=="cox" & criterion!="BIC+"){
if((method!="MS.criterion" & method!="LASSO")){mycoefficients <- mycoefficients[-1,]}
if((method!="MS.criterion" & method!="LASSO")){mycoefficients.b <- mycoefficients.b[-1,]}
if((method!="MS.criterion" & method!="LASSO")){mycoefficients.s <- mycoefficients.s[-1,]}
if((method!="MS.criterion" & method!="LASSO")){mycoefficients.b.s <- mycoefficients.b.s[-1,]}
}
if(report.exp==TRUE){
if(is.null(mycoefficients)==FALSE){mycoefficients <- round(exp(mycoefficients),digits=6)}
if(is.null(mycoefficients.b)==FALSE){mycoefficients.b <- round(exp(mycoefficients.b),digits=6)}
if(is.null(mycoefficients.s)==FALSE){mycoefficients.s <- round(exp(mycoefficients.s),digits=6)}
if(is.null(mycoefficients.b.s)==FALSE){mycoefficients.b.s <- round(exp(mycoefficients.b.s),digits=6)}
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
          boot.results = list(result.boot.s,result.boot.ma),
          setup = list(missing.data,model,outcome,id,method,criterion,kfold,cvr,
                       inference,B,CI,var.remove,add.factor,add.interaction, 
                       add.transformation,add.stratum,screen,exclude,report.exp)
)

class(res) <- "mami"
res



##########################

}





