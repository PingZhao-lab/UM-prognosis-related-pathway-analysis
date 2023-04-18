#从12个通路中筛选8个通路的折线图
load("D:/1学习/eye/New_project/scRNAseq_proc2/outdata/data_5cv1.Rdata")

library(ggplot2)
library(splines)

p <- ggplot(data_5cv1, aes(otus, value)) +
  geom_smooth(se = FALSE,  method = 'glm', formula = y~ns(x, 6),size=2,color="#81B5D8") +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'),
        axis.text=element_text(size=30,face="bold"),axis.title=element_text(size=30,face="bold"),panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid")) +
  labs(title = '',x = "Pathways", y = "Cross-validation error")

p
p + geom_vline(xintercept = 8,size =1.5, colour="#D25B45", linetype="dashed")


#重要性
load("D:/1学习/eye/New_project/scRNAseq_proc2/outdata/Path_importance.Rdata")
library(randomForest)
Path_importance <- Path_importance[order(Path_importance$MeanDecreaseGini,decreasing = T),]
a <- Path_importance[1:8,]
a$path <- c("Pyrimidine metabolism [PATH:hsa00240]","Epithelial cell signaling in Helicobacter pylori infection [PATH:hsa05120]"
                  ,"Transfer RNA biogenesis [BR:hsa03016]"
                  ,"Glycosaminoglycan binding proteins [BR:hsa00536]"
                  ,"Amino sugar and nucleotide sugar metabolism [PATH:hsa00520]"
                  ,"Peptidases [BR:hsa01002]"
                  ,"Cellular senescence [PATH:hsa04218]"
                  ,"Progesterone-mediated oocyte maturation [PATH:hsa04914]")

a$path <- factor(a$path,levels=a$path)
a <- a[order(a$MeanDecreaseGini),]
colnames(a )
library(ggplot2)
library(plotly)
#library(gapminder)
p <- ggplot(a, aes(MeanDecreaseGini, path, size = 30, color="continent")) +
  geom_point() +
  theme_bw()+theme(axis.text=element_text(size=10,face="bold"),axis.title=element_text(size=20),
                   legend.position = 'none',panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))+
  labs(title = '', y = NULL)+xlab('MeanDecreaseGini')+scale_y_discrete(labels = function(y) str_wrap(y, width = 30) )

p


#画ROC
library(ggplot2)
library(pROC)
library(randomForest)

load("D:/1学习/eye/New_project/结果图/ROC/dif8path_roc_test.Rdata")
load("D:/1学习/eye/New_project/结果图/ROC/Sig3path_roc.Rdata")
load("D:/1学习/eye/New_project/结果图/ROC/UM8path_roc_test.Rdata")

res <- list(UM=UM8path_roc_test, DF=dif8path_roc_test,Sig3=Sig3path_roc)
ROC_plot <- ggroc(res, legacy.axes = TRUE, size = 2)+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype=4)+
  theme_bw()+ggtitle('')+theme(axis.text=element_text(size=25,face="bold"),axis.title=element_text(size=30,face="bold"),
                                  legend.position = 'none',panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
ROC_plot
save(ROC_plot,file = "data/ROC_plot.Rdata")



