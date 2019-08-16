
args <- commandArgs(trailingOnly = TRUE)
conFile <- args[1]
corFile <- args[2]
newFile <- args[3]
outFile1 <- args[4]
outFile2 <- args[5]
#setwd("/srv/netscratch/dep_coupland/grp_schneeberger/projects/AMPRILdenovo/genefamily/pangenome")
pdf(outFile1)

#par(mfrow=c(2,1))
conn <- read.table(conFile,header = FALSE)
plot(conn[,1],conn[,2] , xlab="Number of genomes", ylab= "Core genome  | Pan-genome  ", cex=0.3, xlim = c(0,12) ,ylim = c(24000,31000))
xvalues<-conn[,1]
yvalues<-conn[,2]
points(tapply(conn[,2],conn[,1],median),col="red",pch=15, cex=0.5)

####least squares fit of the exponential decay 
model_1 <- nls(yvalues ~ A*exp(B*xvalues)+C,start=list(A=800,B=-0.3,C=24000))
model_1
summary(model_1)

##least squares fit of the power low
model_2 <- nls(yvalues ~ A*xvalues^(-B),start=list(A=24000,B=-0.3))
model_2
summary(model_2)

new.data <- data.frame(xvalues = seq(min(xvalues),12,len = 100))
lines(new.data$xvalues,predict(model_1,newdata = new.data),col="blue")
lines(new.data$xvalues,predict(model_2,newdata = new.data),col="red")
legend(8, 31000, legend=c("exponential fit", "power law fit"),
       col=c("blue", "red"),  lty=c(1,1), cex=0.8, bty="n")
med<-c()


conn <- read.table(corFile,header = FALSE)
points(conn[,1],conn[,2],cex=0.4)
xvalues<-conn[,1]
yvalues<-conn[,2]
points(tapply(conn[,2],conn[,1],median),col="red",pch=15, cex=0.5)

####least squares fit of the exponential decay 
model_1 <- nls(yvalues ~ A*exp(B*xvalues)+C,start=list(A=800,B=-0.3,C=24000))
model_1
summary(model_1)

##least squares fit of the power low
model_2 <- nls(yvalues ~ A*xvalues^(-B),start=list(A=24000,B=-0.3))
model_2
summary(model_2)

new.data <- data.frame(xvalues = seq(min(xvalues),12,len = 100))
lines(new.data$xvalues,predict(model_1,newdata = new.data),col="blue")
lines(new.data$xvalues,predict(model_2,newdata = new.data),col="red")
dev.off()


############### new genes ##############
pdf(outFile2)
conn2 <- read.table(newFile,header = FALSE)
plot(conn2[,1],conn2[,2] , xlab="Number of genomes", ylab= "New genes", col="black", cex=0.4, xlim = c(0,12) ,ylim = c(0,1000))
points(seq(2,max(conn2[,1]),1),tapply(conn2[,2],conn2[,1],median),col="red",pch=15, cex =0.4)
xvalues<-conn2[,1]
yvalues<-conn2[,2]
##least squares fit of the exponential decay 
model_1 <- nls(yvalues ~ A*exp(B*xvalues)+C,start=list(A=800,B=-0.3,C=800))
model_1
summary(model_1)

##least squares fit of the power low
model_2 <- nls(yvalues ~ A*xvalues^(-B),start=list(A=800,B=-0.3))
model_2
summary(model_2)
new.data <- data.frame(xvalues = seq(min(xvalues),12,len = 100))
lines(new.data$xvalues,predict(model_1,newdata = new.data),col="blue")
lines(new.data$xvalues,predict(model_2,newdata = new.data),col="red")
legend("right", legend=c("exponential fit", "power law fit"),
       col=c("blue", "red"), lty=1:1, cex=0.8)

#abline(lm(conn[,2]~conn[,1]), col="red")
dev.off()

