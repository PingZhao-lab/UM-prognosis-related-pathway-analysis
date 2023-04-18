load("D:/1学习/eye/Scissor数据整合/GSM4147097/data/data_mean_jj.Rdata")
surdata<-read.table("D:/1学习/eye/RNAseq/data/TCGA-UVM.survival.tsv",sep = "\t",header = T)
surdata<-surdata[order(surdata$OS.time),]
sursa<-gsub("-",".",surdata$sample)


pp<-match(surdata$sample,colnames(data_mean_jj))
colnames(data_mean_jj)[pp]
data_mean_jj<-data_mean_jj[,pp]#样本按照生存时间从低到高排列
dim(data_mean_jj)

# library(pheatmap)
# pheatmap(data_mean_jj,cluster_rows=T,cluster_cols=F,show_rownames=F,show_colnames=F,scale = "column",
#          legend=F,border_color="black")#热图不好

#样本聚类
#PDF尺寸 9.76 8.05
library(factoextra)
library(gridExtra)
library(tidyverse)
library(cluster)

df <- scale(t(data_mean_jj))
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


#TCGA基因聚类####
df <- scale(data_mean_jj)
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

#两类基因top1功能气泡图
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


