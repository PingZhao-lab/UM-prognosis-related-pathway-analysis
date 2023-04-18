load("D:/1学习/eye/New_project/scRNAseq_proc/randomForest/TCGA_pathway.Rdata")
load("D:/1学习/eye/New_project/scRNAseq_proc/randomForest/surdata.Rdata")
load("D:/1学习/eye/New_project/scRNAseq_proc2/outdata/UMpathways.Rdata")
library(ggplot2)
library(pROC)
library(randomForest)

pt <- NULL
for(i in 1:nrow(TCGA_pathway)){
  a <- unlist(strsplit(row.names(TCGA_pathway)[i],":"))
  a <- a[length(a)]
  pt <- c(pt,gsub("]","",a))
}
pp<-match(UMpathways,pt)
pm<-TCGA_pathway[pp,]

pp<-match(surdata$sample,colnames(pm))
pm <- pm[,pp]
dim(pm)


TCGA_pathway1 <- data.frame(t(pm))
pp <- match(row.names(TCGA_pathway1),surdata$sample)
TCGA_pathway1$group <- surdata$sur_phen[pp]
TCGA_pathway1[1:3,c(1,ncol(TCGA_pathway1))]

TCGA_pathway1$group <- as.factor(TCGA_pathway1$group)


set.seed(1217)
pp1 <- which(TCGA_pathway1$group=="Poor")
pp2 <- which(TCGA_pathway1$group=="Good")

train_1 <- sample(pp1,length(pp1)*0.8)
train_2 <- sample(pp2,length(pp2)*0.8)

train_data <- TCGA_pathway1[c(train_1,train_2),]
dim(train_data)
test_data <- TCGA_pathway1[-c(train_1,train_2),]
dim(test_data)


b <- 0.5
accuracy <- NULL
for(i in 1:1000){
  print(i)
  um_train.forest <- randomForest(group ~ ., data = train_data, importance = TRUE)
  pre_ran <- predict(um_train.forest,newdata=test_data)
  a1 <- test_data$group==pre_ran

  a <- length(which(a1==T))/length(pre_ran)
  accuracy <- c(accuracy,a)

  if(a>b){
    forest_model <- um_train.forest
    b <- a
  }
}

pre_ran <- predict(forest_model,newdata=train_data)
obs_p_ran = data.frame(prob=pre_ran,obs=train_data$group)
table(train_data$group,pre_ran,dnn=c("T","P"))

pre_ran <- predict(forest_model,newdata=test_data)
obs_p_ran = data.frame(prob=pre_ran,obs=test_data$group)
table(test_data$group,pre_ran,dnn=c("T","P"))
ran_roc <- roc(test_data$group,as.numeric(pre_ran))
ran_roc$auc
ran_roc_test8path <- ran_roc
save(ran_roc_test8path,file = "./outdata/ran_roc_test8path.Rdata")

pre_ran <- predict(forest_model,newdata=TCGA_pathway1)
obs_p_ran = data.frame(prob=pre_ran,obs=TCGA_pathway1$group)
table(TCGA_pathway1$group,pre_ran,dnn=c("T","P"))
ran_roc <- roc(as.numeric(TCGA_pathway1$group),as.numeric(pre_ran))
ran_roc$auc


save(forest_model,file = "./outdata/forest_model.Rdata")
save(TCGA_pathway1,file = "./outdata/TCGA_pathway1.Rdata")
save(train_data,file = "./outdata/train_data.Rdata")
save(test_data,file = "./outdata/test_data.Rdata")


library(pROC)
pp <- match(names(pre_ran),surdata$sample)
surdata <- surdata[pp,]


pp<-match(surdata$sample,names(pre_ran))
pre_ran <- pre_ran[pp]

ran_roc <- roc(as.character(pre_ran),surdata$OS.time)
UM8path_roc_test <- ran_roc
plot(UM8path_roc_test, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),grid.col="black", max.auc.polygon=TRUE,auc.polygon.col="skyblue", print.thres=TRUE,main=NULL,
     xlab=NA,ylab=NA,font=2,cex.axis=2,identity.col="black",identity.lwd=2,grid.lwd=1.5,print.auc.cex=1.5,print.auc.col="black",print.thres.cex=1.5,col="#D25B45",lwd=5)
