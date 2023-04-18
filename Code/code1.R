#TCGA data####
data<-read.table("data/TCGA-UVM.htseq_fpkm.tsv",sep = "\t",header = T)
ann<-read.table("data/gencode.v22.annotation.gene.probeMap",sep = "\t",header = T)
pp<-match(data$Ensembl_ID,ann$id)
table(is.na(pp))
genes<-ann$gene[pp]
table(duplicated(genes))
table(is.na(genes))
data$Ensembl_ID<-genes
data_mean=aggregate(.~Ensembl_ID,mean,data=data)
row.names(data_mean)<-data_mean$Ensembl_ID
data_mean<-data_mean[,-1]

#Filter
qc<-NULL
for(i in 1:nrow(data_mean)){
  if(length(which(as.numeric(data_mean[i,])==0))/80>=0.33){
    qc<-c( qc,i)
  }
}
data_mean<-data_mean[-qc,]
save(data_mean,file = "outdata/data_mean.Rdata")

#K-M Survival Analysis
library(survival)
library(survminer)
surdata<-read.table("data/TCGA-UVM.survival.tsv",sep = "\t",header = T)

sursa<-gsub("-",".",surdata$sample)
Pvalue<-NULL
for(i in 1:nrow(data_mean)){
  print(i)
  group<-rep("nn",80)
  
  x<-as.numeric(data_mean[i,])
  zh<-median(x)
  h_g<-colnames(data_mean)[which(x>=zh)]
  l_g<-colnames(data_mean)[which(x<zh)]
  
  pp<-match(h_g,sursa)
  group[pp]<-"high"
  pp<-match(l_g,sursa)
  group[pp]<-"low"
  
  surdata$group<-group
  sfit <- survdiff(Surv(OS.time, OS)~group, data=surdata)
  #ggsurvplot(survfit(Surv(OS.time, OS)~group, data=surdata),pval=T)
  Pvalue[i] <- 1 - pchisq(sfit$chisq, length(sfit$n) -1)
}
survival_result<-data.frame(Genes=row.names(data_mean),Pvalue,FDR=p.adjust(Pvalue,method = "BH"))
length(which(survival_result$FDR<=0.05))
save(survival_result,file = "outdata/survival_result.Rdata")

#GSE176345####
GSE176345_data<-read.table("D:/eye/RNAseq_DEG_GSE176345/data/GSE176345_No92_1VSPig1.All.txt",header = T,sep = "\t")


#Survival information genes####
scgenes<-survival_result$Genes[which(survival_result$FDR<=0.05)]
dfgenes<-GSE176345_data$gene_id[which(GSE176345_data$padj<=0.05)]#这个是1480
SIGs<-intersect(scgenes,dfgenes)
save(SIGs,file = "outdata/SIGs.Rdata")

#Sample clustering####
library(factoextra)
library(gridExtra)
library(tidyverse)
library(cluster)

df <- scale(t(data_mean[match(SIGs,row.names(data_mean)),]))
dim(df)

fviz_nbclust(df, kmeans, method = "wss") + geom_vline(xintercept = 2, linetype = 2,colour="steelblue")+theme(text = element_text(size=20))#使用
fviz_nbclust(df, kmeans, method = "silhouette",print.summary=F)+theme(text = element_text(size=20))#使用
# fviz_nbclust(df, kmeans, method = "gap_stat",print.summary=F)+theme(text = element_text(size=20))


set.seed(123)
km_result <- kmeans(df, 2, nstart = 25)
fviz_cluster(km_result, data = df,
             palette = c("#E7D1A9", "#B2C69A"),ellipse=T,geom="point",pointsize=2.5,
             ellipse.type = "norm",
             star.plot = TRUE,
             repel = TRUE,
             ggtheme = theme_minimal()
)+theme(text = element_text(size=30))+labs(title = "UM samples")

tq<-which(km_result[["cluster"]]==1)
c1<-names(km_result[["cluster"]])[tq]
sum(surdata$OS.time[match(c1,surdata$sample)])

tq<-which(km_result[["cluster"]]==2)
c2<-names(km_result[["cluster"]])[tq]
sum(surdata$OS.time[match(c2,surdata$sample)])

t.test(surdata$OS.time[match(c1,surdata$sample)],surdata$OS.time[match(c2,surdata$sample)]) #P值为0.005346，2颗星
# wilcox.test(surdata$OS.time[match(c2,surdata$sample)],surdata$OS.time[match(c1,surdata$sample)],alternative  = "greater")#0.004

sample_data<-data.frame(V=c(surdata$OS.time[match(c1,surdata$sample)],surdata$OS.time[match(c2,surdata$sample)]),
                        C=c(rep("Cluster 1",length(c1)),rep("Cluster 2",length(c2))))
sample_data$C<-factor(sample_data$C,levels = c("Cluster 1","Cluster 2"))

ggplot(sample_data, aes(x=C, y=V)) +
  geom_boxplot(alpha=1,fill=c("#E7D1A9","#B2C69A"),size=1.5)+
  xlab(label = NULL)+
  ylab(label = "Survival time")+
  theme_classic()+
  theme(text = element_text(color ="black",size=30),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),axis.line=element_line(size=1.5))


#Gene clustering####
df <- scale(data_mean[match(SIGs,row.names(data_mean)),])
dim(df)

fviz_nbclust(df, kmeans, method = "wss") + geom_vline(xintercept = 2, linetype = 2,colour="steelblue")+theme(text = element_text(size=20))
fviz_nbclust(df, kmeans, method = "silhouette",print.summary=F)+theme(text = element_text(size=20))
# fviz_nbclust(df, kmeans, method = "gap_stat",print.summary=F)+theme(text = element_text(size=20))

set.seed(123)
km_result <- kmeans(df, 2, nstart = 25)
fviz_cluster(km_result, data = df,
             palette = c("#9DC7DB", "#DEA49D"),ellipse=T,geom="point",pointsize=2.5,
             ellipse.type = "norm",
             star.plot = TRUE,
             repel = TRUE,
             ggtheme = theme_minimal()
)+theme(text = element_text(size=30))+labs(title = "Survival information genes")


table(km_result[["cluster"]])
Cluster1_genes<-names(km_result[["cluster"]][which(km_result[["cluster"]]==1)])
Cluster2_genes<-names(km_result[["cluster"]][which(km_result[["cluster"]]==2)])
write.csv(Cluster1_genes,file = "outdata/Cluster1_genes.csv")
write.csv(Cluster2_genes,file = "outdata/Cluster2_genes.csv")

save(Cluster1_genes,file = "outdata/Cluster1_genes.Rdata")
save(Cluster2_genes,file = "outdata/Cluster2_genes.Rdata")

#Functional enrichment####
Cluster1_function <- read.table("./data/Cluster1_function.txt",sep = "\t",header = T)
Cluster1_function <- Cluster1_function[which(Cluster1_function$PValue<0.05),]
Cluster1_function$Term <- gsub("~",": ",Cluster1_function$Term)
Cluster1_function$PValue <- -log10(Cluster1_function$PValue)
Cluster1_function <- Cluster1_function[order(Cluster1_function$PValue),]
Cluster1_function$Term<-factor(Cluster1_function$Term,levels = Cluster1_function$Term)

library(ggplot2)
#Cluster1_genes_function
ggplot(Cluster1_function,aes(PValue,Term,size=Count))+geom_point()+geom_point(colour = "#99C2D5")+
  labs(x="-log10(P values)",y=NULL)+scale_y_discrete(labels = function(y) str_wrap(y, width = 40) )+theme(text=element_text(size=25))#16


Cluster2_function <- read.table("./data/Cluster2_function.txt",sep = "\t",header = T)
Cluster2_function <- Cluster2_function[which(Cluster2_function$PValue<0.05),]
Cluster2_function$Term <- gsub("~",": ",Cluster2_function$Term)
Cluster2_function$PValue <- -log10(Cluster2_function$PValue)
Cluster2_function <- Cluster2_function[order(Cluster2_function$PValue),]
Cluster2_function$Term<-factor(Cluster2_function$Term,levels = Cluster2_function$Term)

library(ggplot2)
#Cluster2_genes_function
ggplot(Cluster2_function,aes(PValue,Term,size=Count))+geom_point()+geom_point(colour = "#D8A099")+
  labs(x="-log10(P values)",y=NULL)+scale_y_discrete(labels = function(y) str_wrap(y, width = 30) )+theme(text=element_text(size=30))# 9.93 8


