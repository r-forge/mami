jma <- function(y,x,ma.method=c("JMA","MMA"),model.subset=c("nested","all"),factor.variables=NULL,pd=T,
                calc.var=c("none","boot"), bsa=100){ 

ma.method            <- match.arg(ma.method)
model.subset         <- match.arg(model.subset) 
calc.var             <- match.arg(calc.var) 
if(model.subset=="all" & ncol(x)>20){stop("Too many model combinations: \n Use either model.subset='nested' or reduce number of covariates")}
if(calc.var=="boot" & pd==TRUE){cat(paste("Standard errors for", ma.method, "are based on", bsa, "bootstrap samples. \n"))}
y <- data.frame(y)
x <- data.frame(x)
  
tf <- function(mymatrix){colnames(mymatrix)[sapply(mymatrix,is.factor)]}
covariates <- colnames(x)
myindependent <- colnames(y)
if(length(tf(x))>0){
factor.variables <- c(factor.variables,tf(x))
factor.variables <- intersect(factor.variables,factor.variables)
factor.variables <- intersect(colnames(x),factor.variables)
}
if(is.null(factor.variables)==FALSE){
if(pd==T){cat(paste("Note: The following variables are treated as factors and are recoded into dummies:", paste(factor.variables,collapse=" ")," \n\n"))}
x <- as.data.frame(x)
for(i in match(factor.variables,covariates)){
x[,covariates[i]] <- as.factor(x[,covariates[i]])
}
myformula <- as.formula(paste(myindependent,"~", paste(covariates, collapse="+")))
x <- model.matrix(myformula,data=as.data.frame(cbind(y,x)))[,-1]
x<- as.matrix(x)
}
keepnames <- colnames(x)
x <- as.matrix(cbind(1,x))
y <- as.matrix(y)
n <- nrow(x)
p <- ncol(x)

if (model.subset == "nested"){
        s <- matrix(1,nrow=p,ncol=p)
        s[upper.tri(s)] <- 0
        zero <- matrix(0,nrow=1,ncol=p)
        s <- rbind(zero,s)
        s <- s[s[,1]==1,]  
     } 
if (model.subset == "all"){
        s <- matrix(0,nrow=2^p,ncol=p)
        s0 <- matrix(c(1,rep(0,p-1)),1,p)
        s1 <- matrix(c(rep(0,p)),1,p)
        for (i in 2:2^p){
           s1 <- s0 + s1
           for (j in 1:p){
              if (s1[1,j] == 2){
                 s1[1,j+1] <- s1[1,j+1]+1
                 s1[1,j] <- 0
              }
           }           
           s[i,] <- s1
        }
       s <- s[s[,1]==1,]    
     }   
    

  m <- nrow(s)
  bbeta <- matrix(0,nrow=p,ncol=m)
  if (ma.method == "JMA") ee <- matrix(0,nrow=n,ncol=m)

  for (j in 1:m)try({
     ss <- matrix(1,nrow=n,ncol=1) %*% s[j,]
     indx1 <- which(ss[,]==1)
     xs <- as.matrix(x[indx1])
     xs <- matrix(xs,nrow=n,ncol=nrow(xs)/n)
     if (sum(ss)==0){
        xs <- x
        betas <- matrix(0,nrow=p,ncol=1)
        indx2 <- matrix(c(1:p),nrow=p,ncol=1)  
     }  
     if (sum(ss)>0){   
        betas <- solve(t(xs)%*%xs)%*%t(xs)%*%y
        indx2 <- as.matrix(which(s[j,]==1))  
     }
     beta0 <- matrix(0,nrow=p,ncol=1)
     beta0[indx2] <- betas     
     bbeta[,j] <- beta0    
     if (ma.method == "JMA"){
        ei <- y - xs %*% betas
        hi <- diag(xs %*% solve(t(xs)%*%xs) %*% t(xs))
        ee[,j] <- ei*(1/(1-hi))
     }
  },silent=TRUE)

  if (ma.method == "MMA"){
     ee <- y %*% matrix(1,nrow=1,ncol=m) - x %*% bbeta
     ehat <- y - x %*% bbeta[,m]
     sighat <- (t(ehat) %*% ehat)/(n-p)
  }
  
  a1 <- t(ee) %*% ee
  if (qr(a1)$rank<ncol(ee)) a1 <- a1 + diag(m)*1e-10
  if (ma.method == "MMA") a2 <- matrix(c(-c(sighat)*rowSums(s)),m,1)  
  if (ma.method == "JMA") a2 <- matrix(0,nrow=m,ncol=1)
  a3 <- t(rbind(matrix(1,nrow=1,ncol=m),diag(m),-diag(m)))
  a4 <- rbind(1,matrix(0,nrow=m,ncol=1),matrix(-1,nrow=m,ncol=1))

  w0 <- matrix(1,nrow=m,ncol=1)/m 
  QP <- try(quadprog::solve.QP(a1,a2,a3,a4,1),silent=TRUE)
  if(class(QP)=="try-error"){QP <- try(quadprog::solve.QP(make.positive.definite(a1),a2,a3,a4,1),silent=TRUE)
  if(pd==T){cat("Note: the residual matrix was not positive definite and has been adjusted. \n")}
  }
  w <- QP$solution
  w <- as.matrix(w)
  w <- w*(w>0)
  w <- w/sum(w0)
  betahat <- matrix(bbeta %*% w,ncol=1,dimnames=list(c("Intercept",keepnames),"beta"))
  ybar <- mean(y)
  yhat <- x %*% betahat
  ehat <- y-yhat
  wsummary <- matrix(cbind(s,w),ncol=ncol(cbind(s,w)),dimnames=list(paste("model",1:length(w)),c("Intercept",keepnames,"weight")))

# Variance
se <- lce <- uce <- matrix(rep(NA,length(betahat)),nrow=1)
if(calc.var=="boot"){

jma.se.boot <- function(mydata,indices){
mydata <- mydata[indices,]
y <- matrix(mydata[,1],ncol=1)
x <- mydata[,-1]
bbeta <- matrix(0,nrow=p,ncol=m)
  if (ma.method == "JMA") ee <- matrix(0,nrow=n,ncol=m)
  for (j in 1:m)try({
     ss <- matrix(1,nrow=n,ncol=1) %*% s[j,]
     indx1 <- which(ss[,]==1)
     xs <- as.matrix(x[indx1])
     xs <- matrix(xs,nrow=n,ncol=nrow(xs)/n)
     if (sum(ss)==0){
        xs <- x
        betas <- matrix(0,nrow=p,ncol=1)
        indx2 <- matrix(c(1:p),nrow=p,ncol=1)  
     }  
     if (sum(ss)>0){   
        betas <- solve(t(xs)%*%xs)%*%t(xs)%*%y
        indx2 <- as.matrix(which(s[j,]==1))  
     }
     beta0 <- matrix(0,nrow=p,ncol=1)
     beta0[indx2] <- betas     
     bbeta[,j] <- beta0    
     if (ma.method == "JMA"){
        ei <- y - xs %*% betas
        hi <- diag(xs %*% solve(t(xs)%*%xs) %*% t(xs))
        ee[,j] <- ei*(1/(1-hi))
     }
  },silent=TRUE)
  if (ma.method == "MMA"){
     ee <- y %*% matrix(1,nrow=1,ncol=m) - x %*% bbeta
     ehat <- y - x %*% bbeta[,m]
     sighat <- (t(ehat) %*% ehat)/(n-p)
  }
  a1 <- t(ee) %*% ee
  if (qr(a1)$rank<ncol(ee)) a1 <- a1 + diag(m)*1e-10
  if (ma.method == "MMA") a2 <- matrix(c(-c(sighat)*rowSums(s)),m,1)  
  if (ma.method == "JMA") a2 <- matrix(0,nrow=m,ncol=1)
  a3 <- t(rbind(matrix(1,nrow=1,ncol=m),diag(m),-diag(m)))
  a4 <- rbind(1,matrix(0,nrow=m,ncol=1),matrix(-1,nrow=m,ncol=1))
  w0 <- matrix(1,nrow=m,ncol=1)/m 
  QP <- try(quadprog::solve.QP(a1,a2,a3,a4,1),silent=TRUE)
  if(class(QP)=="try-error"){QP <- try(quadprog::solve.QP(make.positive.definite(a1),a2,a3,a4,1),silent=TRUE)
  if(pd==T){cat("Note: the residual matrix was not positive definite and has been adjusted. \n")}
  }
  w <- QP$solution
  w <- as.matrix(w)
  w <- w*(w>0)
  w <- w/sum(w0)
  c(bbeta %*% w)
}
result.boot <- try(boot(cbind(y,x),jma.se.boot,bsa),silent=TRUE)
mysd <- function(myvalues){sd(na.omit(myvalues))}
se <- matrix(apply(result.boot$t,2,mysd),nrow=1)
lower95 <- function(xx){quantile(na.omit(xx),probs=0.025)}
upper95 <- function(xx){quantile(na.omit(xx),probs=0.975)}
lce <- matrix(apply(result.boot$t,2,lower95),nrow=1)
uce <- matrix(apply(result.boot$t,2,upper95),nrow=1)
} 

rownames(lce)<- paste("95% CI (lower)")
rownames(uce)<- paste("95% CI (upper)")  
rownames(se) <- "se"
colnames(lce) <- c("Intercept",keepnames)
colnames(uce) <- c("Intercept",keepnames)
colnames(se) <- c("Intercept",keepnames)
  
  # results
  res <- list(betahat=betahat,se=se,lci=lce,uci=uce,weight=wsummary,yhat=yhat,ehat=ehat,y=y,x=x)
  class(res) <- "jma"
  res  
}


