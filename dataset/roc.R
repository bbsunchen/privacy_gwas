#Copyright 2015, Chen Sun(bbsunchen@outlook.com)

setwd("C:\\Users\\ChenSun\\Dropbox\\privacy")
library("flux") # calculate auc

# reproduce Sankararaman et al's Result
lr1000 <- read.table("lr1000.roc", sep="\t", header=TRUE)
lr10000<- read.table("lr10000.roc", sep="\t", header=TRUE)
plot(log(lr1000$fp,10), lr1000$tp, type="l", cex.main=2, cex.lab=1.6, cex.axis=1.6, lwd=7, col="red", xlab="False-positive rate(log 10)", ylab="Power(True-positive rate)")
lines(log(lr10000$fp,10), lr10000$tp, lwd=7, col="black")
auc_lr1000 <- auc(log(lr1000$fp,10), lr1000$tp)
auc_lr10000 <- auc(log(lr10000$fp,10), lr10000$tp)
text(x=-1.5,y=0.6, "SNPs=10000", col="black", cex=1.5)
text(x=-0.5,y=0.35, "SNPs=1000",col="red", cex=1.5)




# LR test for different frequency level
plot(log(lr1000$fp,10), lr1000$tp, type="l", cex.main=1.6, cex.lab=1.6, cex.axis=1.6, lwd=7, col="white", xlab="False-positive rate(log 10)", ylab="Power(True-positive rate)")
lines(log(lr10000$fp,10), lr10000$tp, lwd=7,col="black")
auc_lr10000 <- auc(log(lr10000$fp,10), lr10000$tp)

lr_round001<- read.table("lr_round001.roc", sep="\t", header=TRUE)
lines(log(lr_round001$fp,10), lr_round001$tp, lwd=7,col="red")
auc_lr_round001 <- auc(log(lr_round001$fp,10), lr_round001$tp)

lr_round01<- read.table("lr_round01.roc", sep="\t", header=TRUE)
lines(log(lr_round01$fp,10), lr_round01$tp, lwd=7,col="orange")
auc_lr_round01 <- auc(log(lr_round01$fp,10), lr_round01$tp)

lr_gaussian<- read.table("lr_gaussian.roc", sep="\t", header=TRUE)
lines(log(lr_gaussian$fp,10), lr_gaussian$tp, lwd=7,col="blue")
auc_lr_gaussian <- auc(log(lr_gaussian$fp,10), lr_gaussian$tp)

text(x=-0.8,y=0.3, cex=1.2,paste("AUC of Pretreatment P0=", round(auc_lr10000,6)), col="black")
text(x=-0.8,y=0.26, cex=1.2,paste("AUC of Pretreatment P1=", round(auc_lr_round001,6)),col="red")
text(x=-0.8,y=0.22, cex=1.2,paste("AUC of Pretreatment P2=", round(auc_lr_round01,6)),col="orange")
text(x=-0.8,y=0.18, cex=1.2,paste("AUC of Pretreatment P3=", round(auc_lr_gaussian,6)),col="blue")

##################################################################
# T1 test for different frequency pretreatment
# this is for the coordiate range, please do not delete
plot(log(lr1000$fp,10), lr1000$tp, type="l", cex.main=1.6, cex.lab=1.6, cex.axis=1.6, lwd=7, col="white", xlab="False-positive rate(log 10)", ylab="Power(True-positive rate)")

t1_direct<- read.table("t1_direct.roc", sep="\t", header=TRUE)
lines(log(t1_direct$fp,10), t1_direct$tp, lwd=7,col="black")
auc_t1_direct <- auc(log(t1_direct$fp,10), t1_direct$tp)

t1_round001<- read.table("t1_round001.roc", sep="\t", header=TRUE)
lines(log(t1_round001$fp,10), t1_round001$tp, lwd=7,col="red")
auc_t1_round001 <- auc(log(t1_round001$fp,10), t1_round001$tp)

t1_round01<- read.table("t1_round01.roc", sep="\t", header=TRUE)
lines(log(t1_round01$fp,10), t1_round01$tp, lwd=7,col="orange")
auc_t1_round01 <- auc(log(t1_round01$fp,10), t1_round01$tp)

t1_gaussian<- read.table("t1_gaussian.roc", sep="\t", header=TRUE)
lines(log(t1_gaussian$fp,10), t1_gaussian$tp, lwd=7,col="green")
auc_t1_gaussian <- auc(log(t1_gaussian$fp,10), t1_gaussian$tp)

text(x=-0.8,y=0.3, cex=1.2,paste("AUC of Pretreatment P0=", round(auc_t1_direct,6)), col="black")
text(x=-0.8,y=0.26, cex=1.2,paste("AUC of Pretreatment P1=", round(auc_t1_round001,6)),col="red")
text(x=-0.8,y=0.22, cex=1.2,paste("AUC of Pretreatment P2=", round(auc_t1_round01,6)),col="orange")
text(x=-0.8,y=0.18, cex=1.2,paste("AUC of Pretreatment P3=", round(auc_t1_gaussian,6)),col="green")

##################################################################
# T2 test for different frequency pretreatment
# this is for the coordiate range, please do not delete
plot(log(lr1000$fp,10), lr1000$tp, type="l", cex.main=1.6, cex.lab=1.6, cex.axis=1.6, lwd=7, col="white", xlab="False-positive rate(log 10)", ylab="Power(True-positive rate)")

t2_direct<- read.table("t2_direct.roc", sep="\t", header=TRUE)
lines(log(t2_direct$fp,10), t2_direct$tp, lwd=7,col="black")
auc_t2_direct <- auc(log(t2_direct$fp,10), t2_direct$tp)

t2_round001<- read.table("t2_round001.roc", sep="\t", header=TRUE)
lines(log(t2_round001$fp,10), t2_round001$tp, lwd=7,col="red")
auc_t2_round001 <- auc(log(t2_round001$fp,10), t2_round001$tp)

t2_round01<- read.table("t2_round01.roc", sep="\t", header=TRUE)
lines(log(t2_round01$fp,10), t2_round01$tp, lwd=7,col="orange")
auc_t2_round01 <- auc(log(t2_round01$fp,10), t2_round01$tp)

t2_gaussian<- read.table("t2_gaussian.roc", sep="\t", header=TRUE)
lines(log(t2_gaussian$fp,10), t2_gaussian$tp, lwd=7,col="green")
auc_t2_gaussian <- auc(log(t2_gaussian$fp,10), t2_gaussian$tp)

text(x=-0.8,y=0.3, cex=1.2,paste("AUC of Pretreatment P0=", round(auc_t2_direct,6)), col="black")
text(x=-0.8,y=0.26, cex=1.2,paste("AUC of Pretreatment P1=", round(auc_t2_round001,6)),col="red")
text(x=-0.8,y=0.22, cex=1.2,paste("AUC of Pretreatment P2=", round(auc_t2_round01,6)),col="orange")
text(x=-0.8,y=0.18, cex=1.2,paste("AUC of Pretreatment P3=", round(auc_t2_gaussian,6)),col="green")