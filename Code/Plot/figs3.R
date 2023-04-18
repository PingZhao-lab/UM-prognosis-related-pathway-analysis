#COX模型的三个检测P值
#8path
# Likelihood ratio test= 33.65  on 8 df,   p=5e-05
# Wald test            = 20.39  on 8 df,   p=0.009
# Score (logrank) test = 26.88  on 8 df,   p=7e-04

#df8path
# Likelihood ratio test= 19.27  on 8 df,   p=4e-04
# Wald test            = 14.88  on 8 df,   p=0.01
# Score (logrank) test = 17.18  on 8 df,   p=7e-04

data <- data.frame(V=c(5e-05,0.009,7e-04,4e-04,0.01,7e-04),M=c("Likelihood ratio test","Wald test","Logrank test","Likelihood ratio test","Wald test","Logrank test"),
                   C=c(rep("Signatures",3),rep("Differences",3)))
data$V <- -log10(data$V)
data$C<-factor(data$C,levels = c("Differences","Signatures"))

library(ggplot2)
ggplot(data,aes(x=C, y=V, fill=C))+
  geom_bar(alpha=1,stat='identity',position="dodge",colour="black",size=1)+
  scale_fill_manual("legend", values = c("Signatures" = "#E27069", "Differences" = "#28ADB4"))+
  labs(x = NULL, y = "-log10(P values)")+
  theme(legend.position="none",text = element_text(color ="black",size=35),
        axis.text.x = element_text(color="black",size=30),
        axis.text.y = element_text(color="black",size=30),
        axis.ticks.x=element_line(color="black",size=2,lineend = 1),
        axis.line.x=element_line(linetype=1,color="black",size=2.5),
        axis.ticks.y=element_line(color="black",size=2,lineend = 1),
        axis.line.y=element_line(linetype=1,color="black",size=2.5),
        plot.title = element_text(size = 45),panel.background = element_blank(),
        strip.text.x = element_text(size = 45)
  )+scale_y_continuous(limits = c(0, 4.5))+facet_grid(. ~ M)

data <- data[which(data$C=="Signatures"),]
ggplot(data,aes(x=M, y=V, fill=C))+
  geom_bar(alpha=1,stat='identity',position="dodge",colour="black",size=1)+
  scale_fill_manual("legend", values = c("Signatures" = "#E27069"))+
  labs(x = NULL, y = "-log10(P values)")+
  theme(legend.position="none",text = element_text(color ="black",size=35),
        axis.text.x = element_text(color="black",size=30),
        axis.text.y = element_text(color="black",size=30),
        axis.ticks.x=element_line(color="black",size=2,lineend = 1),
        axis.line.x=element_line(linetype=1,color="black",size=1.5),
        axis.ticks.y=element_line(color="black",size=2,lineend = 1),
        axis.line.y=element_line(linetype=1,color="black",size=1.5),
        plot.title = element_text(size = 45),panel.background = element_blank(),
        strip.text.x = element_text(size = 45)
  )+scale_y_continuous(limits = c(0, 4.5))


load("D:/1学习/eye/New_project/scRNAseq_proc/randomForest/TCGA_pathway.Rdata")
load("D:/1学习/eye/New_project/scRNAseq_proc/randomForest/surdata.Rdata")
# load("D:/1学习/eye/New_project/scRNAseq_proc2/outdata/UMpathways.Rdata")
#
#
# pt <- NULL
# for(i in 1:nrow(TCGA_pathway)){
#   a <- unlist(strsplit(row.names(TCGA_pathway)[i],":"))
#   a <- a[length(a)]
#   pt <- c(pt,gsub("]","",a))
# }
# pp<-match(UMpathways,pt)
#
# UmPathways_8 <- row.names(TCGA_pathway)[pp]
# save(UmPathways_8,file = "data/UmPathways_8.Rdata")

load("D:/1学习/eye/New_project/结果图/data/UmPathways_8.Rdata")

pp<-match(UmPathways_8,row.names(TCGA_pathway))
TCGA_pathway_8 <- TCGA_pathway[pp,]


library(survival)
library(survminer)
i=1

UmPathways_8[i]
Phenotype<-rep("nn",80)
pp<-which(row.names(TCGA_pathway_8)==UmPathways_8[i])
x<-as.numeric(TCGA_pathway_8[pp,])
zh<-median(x)
h_g<-colnames(TCGA_pathway_8)[which(x>=zh)]
l_g<-colnames(TCGA_pathway_8)[which(x<zh)]
pp<-match(h_g,surdata$sample)
Phenotype[pp]<-"high"
pp<-match(l_g,surdata$sample)
Phenotype[pp]<-"low"

surdata$Phenotype<-Phenotype
sfit <- survdiff(Surv(OS.time, OS)~Phenotype, data=surdata)
ggsurvplot(survfit(Surv(OS.time, OS)~Phenotype, data=surdata),pval=T,title=UmPathways_8[i],palette=c("#DE5D47","#83BCE1"),
           font.title=32,font.x=27,font.y=27,font.tickslab=23,font.legend=27,pval.size=12)

i <- i+1
i
