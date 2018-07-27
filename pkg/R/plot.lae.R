plot.lae <- function(x,xaxis=c("index","realnumbers"),legend.place="topright",display.importance=TRUE,...){

xaxis <- match.arg(xaxis)

if(xaxis=="index"){
myx<-c(seq(1,length(x[[5]])))
myxlab <- "Index of complexity parameter"
cvchoice<-myx[x[[4]]==1]
}
if(xaxis=="realnumbers"){
myx<-x[[5]]
myxlab <- "complexity parameter"
cvchoice<-myx[x[[4]]==1]
}
y1 <- x[[3]]
y2 <- x[[4]]
mymin<-min(x[[1]][-1,c(1,3,5)])
mymax<-max(x[[1]][-1,c(1,3,5)])
mydiff<-max(x[[1]][-1,c(1,3,5)])-min(x[[1]][-1,c(1,3,5)])
t1 <- mymin
t2<- mymin+0.2*mydiff
t3<- mymin+0.4*mydiff
t4<- mymin+0.6*mydiff
t5<- mymin+0.8*mydiff
t6<- mymax

op <- par(ask=TRUE)

plot(myx,y1,type="l",xlab=myxlab,ylab="Shrinkage Averaging Weights",lwd=1.5)
axis(side = 1 , at = cvchoice, labels = F , tick = T , tcl = 1 , lwd.ticks = 2 , col.ticks = "blue")
legend("topright",bty="n",col=c("black","blue"),legend=c("averaging weights","selected by CV"), lty=c(1,1), lwd=1.5, cex=1.25)


par(mai=c(1,1,0.1,1))
plot(0.01 , , ylab="Estimates" , xlab = "Variables", type = "n" , axes=F, xlim=c(1,length(x[[1]][-1,1])) , ylim=c(min(x[[1]][-1,c(1,3,5)]),max(x[[1]][-1,c(1,3,5)])))
lines(c(1:length(rownames(x[[1]])[-1])),x[[1]][-1,1],type="p",pch=17,col="black")
lines(c(1:length(rownames(x[[1]])[-1])),x[[1]][-1,3],type="p",pch=18,col="red")
lines(c(1:length(rownames(x[[1]])[-1])),x[[1]][-1,5],type="p",pch=19,col="blue")
lines(c(0,length(x[[1]][-1,1])),c(0,0),col="grey")
axis(1,las=1,at=c(1:length(rownames(x[[1]])[-1])),labels=rownames(x[[1]])[-1])
axis(2,at=round((c(0,t1,t2,t3,t4,t5,t6)),digits=2))
if(display.importance==TRUE){
legend(legend.place,bty="n",legend=c("LASSO averaging","LASSO","OLS","Variable Importance"),col=c("black","red","blue","chartreuse4"), lwd=1, cex=1)
axis(4,las=1,at=(c(t1,t2,t3,t4,t5,t6)),labels=c(0,0.2,0.4,0.6,0.8,1),col="chartreuse4",col.ticks="chartreuse4",lwd=2,lwd.ticks=2)
lines(c(1:length(rownames(x[[1]])[-1])),x[[2]][-1]*mydiff+mymin,col="chartreuse4",type="p",pch="_",cex=2)
}
par(op)

}