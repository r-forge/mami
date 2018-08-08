plot.mami <- function(x,shade.areas=TRUE,plots.p.page=c("1","4","9","16"),color=c("lightgrey","darkgrey"),adj.bw=1,ask=T,...){

if(is.null(x$boot.results[[1]]) & is.null(x$boot.results[[2]])){stop("Plots are only available after bootstrap estimation")}
plots.p.page = match.arg(plots.p.page)
ine <- function(x) {return(length(x)!=0)}         
par(mfrow=c(sqrt(as.numeric(plots.p.page)),sqrt(as.numeric(plots.p.page))))
         
op <- par(ask=ask)

if(is.null(x$boot.results[[1]])==FALSE){
for(i in (1+(as.numeric(colnames(x$boot.results[[1]])[1]=="(Intercept)"))):(dim(x$boot.results[[1]])[2])){
if(shade.areas==TRUE){
myx <-  c(density(x$boot.results[[1]][,i],adjust=adj.bw)$x[density(x$boot.results[[1]][,i],adjust=adj.bw)$x<0])
myy <-  c(density(x$boot.results[[1]][,i],adjust=adj.bw)$y[density(x$boot.results[[1]][,i],adjust=adj.bw)$x<0])
if(ine(myy)==TRUE){
myx <-  c(myx[1],myx,0)
myy <-  c(0,myy,0)
}
myx2 <- c(density(x$boot.results[[1]][,i],adjust=adj.bw)$x[density(x$boot.results[[1]][,i],adjust=adj.bw)$x>0])
myy2 <- c(density(x$boot.results[[1]][,i],adjust=adj.bw)$y[density(x$boot.results[[1]][,i],adjust=adj.bw)$x>0])
if(ine(myy2)==TRUE){
myx2 <-  c(myx2,myx2[length(myx2)],max(0,myx2[1]))
myy2 <-  c(myy2,0,0)
}
b <- density(x$boot.results[[1]][,i],adjust=adj.bw)
xb <- 0
yb <- 0
if(ine(b$y[b$x<0])==TRUE){
xb <- diff(b$x[b$x<0])
yb <- rollmean(b$y[b$x<0],2)}
ab <- round(sum(xb*yb),digits=2)*100
}
plot(density(x$boot.results[[1]][,i],adjust=adj.bw,...),main="Bootstrap distribution \n (after model selection)",xlab=colnames(x$boot.results[[1]])[i], ylab="Density",lwd=2,...)
if(shade.areas==TRUE){
polygon(myx,myy,col=color[1],border=color[1],...)  
polygon(myx2,myy2,col=color[2],border=color[2],...)
lines(density(x$boot.results[[1]][,i],adjust=adj.bw,...),col="black",lwd=1.25,...)
legend("topright",bty="n",legend=c(paste(ab,"%"),paste((100-ab),"%")), col=color, lwd=7.5, cex=1.25)
}
}}

par(mfrow=c(sqrt(as.numeric(plots.p.page)),sqrt(as.numeric(plots.p.page))))
if(is.null(x$boot.results[[2]])==FALSE){
for(i in (1+(as.numeric(colnames(x$boot.results[[2]])[1]=="(Intercept)"))):(dim(x$boot.results[[2]])[2])){
if(shade.areas==TRUE){
myx <-  c(density(x$boot.results[[2]][,i],adjust=adj.bw)$x[density(x$boot.results[[2]][,i],adjust=adj.bw)$x<0])
myy <-  c(density(x$boot.results[[2]][,i],adjust=adj.bw)$y[density(x$boot.results[[2]][,i],adjust=adj.bw)$x<0])
if(ine(myy)==TRUE){
myx <-  c(myx[1],myx,0)
myy <-  c(0,myy,0)
}
myx2 <- c(density(x$boot.results[[2]][,i],adjust=adj.bw)$x[density(x$boot.results[[2]][,i],adjust=adj.bw)$x>0])
myy2 <- c(density(x$boot.results[[2]][,i],adjust=adj.bw)$y[density(x$boot.results[[2]][,i],adjust=adj.bw)$x>0])
if(ine(myy2)==TRUE){
myx2 <-  c(myx2,myx2[length(myx2)],max(0,myx2[1]))
myy2 <-  c(myy2,0,0)
}
b <- density(x$boot.results[[2]][,i],adjust=adj.bw)
xb <- 0
yb <- 0
if(ine(b$y[b$x<0])==TRUE){
xb <- diff(b$x[b$x<0])
yb <- rollmean(b$y[b$x<0],2)}
ab <- round(sum(xb*yb),digits=2)*100
}
plot(density(x$boot.results[[2]][,i],adjust=adj.bw,...),main="Bootstrap distribution \n (after model averaging)",xlab=colnames(x$boot.results[[2]])[i], ylab="Density",lwd=2,...)
if(shade.areas==TRUE){
polygon(myx,myy,col=color[1],border=color[1],...)  
polygon(myx2,myy2,col=color[2],border=color[2],...)
lines(density(x$boot.results[[2]][,i],adjust=adj.bw,...),col="black",lwd=1.25,...)
legend("topright",bty="n",legend=c(paste(ab,"%"),paste((100-ab),"%")), col=color, lwd=7.5, cex=1.25)
}
}}

par(op)

}