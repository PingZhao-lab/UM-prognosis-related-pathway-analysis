#单细胞降维图


#VN图####
load("./data/Poor.markers.Rdata")
load("./data/Path_survival_result.Rdata")

for(i in 1:nrow(Poor.markers)){
  Poor.markers$path[i] <- unlist(strsplit(Poor.markers$path[i]," \\["))[1]
}

for(i in 1:nrow(Path_survival_result)){
  Path_survival_result$Genes[i] <- unlist(strsplit(Path_survival_result$Genes[i]," \\["))[1]
}


library(VennDiagram)

#Poor.markers是单细胞
D<-venn.diagram(list(TCGA= Path_survival_result$Genes[which(Path_survival_result$Pvalue<0.05)],
                     GSE139829=Poor.markers$path[which(Poor.markers$p_val_adj<0.0001)]),filename=NULL,reverse=TRUE,
                fill=c("#67D5B5","#EE7785"),lty=2,cat.cex=0,cex=2.5,output=TRUE
)
grid.draw(D)
plot(D)


#TCGA热图#####
#load("D:/1学习/eye/New_project/scRNAseq_proc/outdata/UM_pathway_Features.Rdata")
load("D:/1学习/eye/New_project/scRNAseq_proc/randomForest/Candidate_feature_pathway.Rdata")
load("D:/1学习/eye/New_project/scRNAseq_proc/randomForest/TCGA_pathway.Rdata")
load("D:/1学习/eye/New_project/scRNAseq_proc/randomForest/surdata.Rdata")

pp<-match(surdata$sample,colnames(TCGA_pathway))
TCGA_pathway <- TCGA_pathway[,pp]
head(colnames(TCGA_pathway))
head(row.names(TCGA_pathway))
head(UM_pathway_Features)


pn <- NULL
for(i in 1:nrow(TCGA_pathway)){
  pn <- c(pn,unlist(strsplit(row.names(TCGA_pathway)[i]," \\["))[1])
}

pn1 <- NULL
for(i in 1:length(UM_pathway_Features)){
  a <- unlist(strsplit(UM_pathway_Features[i],"\\.\\."))[1]
  a <- gsub("\\."," ",a)
  pn1 <- c(pn1,a)
}

pp <- match(pn1,pn)
pp <- na.omit(pp)
pp <- c(pp,184)
TCGA_pathway1 <- TCGA_pathway[pp,]

library(pheatmap)
pheatmap(TCGA_pathway1,cluster_rows=T,cluster_cols=F,show_rownames=F,show_colnames=F,scale = "column",
         legend=F,border_color="black")#热图不好
