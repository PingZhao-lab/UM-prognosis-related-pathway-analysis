load("D:/1学习/eye/New_project/scRNAseq_proc/data/scRNA_um_pagoda2.Rdata")
load("D:/1学习/eye/New_project/scRNAseq_proc/data/infos1.Rdata")
#
library(dplyr)
library(Seurat)
library(patchwork)
um_pathway_seurat <- CreateSeuratObject(counts = scRNA_um_pagoda2, project = "scRNA", min.cells = 0, min.features = 0)
names(cell_phen) <- colnames(scRNA_um_pagoda2)
um_pathway_seurat@active.ident <- as.factor(cell_phen)

table(cell_phen)
Poor.markers <- FindMarkers(um_pathway_seurat, ident.1 = "Poor",only.pos=T,logfc.threshold=0)
Poor.markers$path <- row.names(Poor.markers)
Pagoda2_scRNA_df_path <- Poor.markers

pathid <- NULL
for(i in 1:nrow(Pagoda2_scRNA_df_path)){
  a <- unlist(strsplit(Pagoda2_scRNA_df_path$path[i],"]"))[1]
  a <- unlist(strsplit(a,"\\["))[2]
  pathid <- c(pathid,a)
}
Pagoda2_scRNA_df_path$pathid <- pathid
pp <- which(is.na(Pagoda2_scRNA_df_path$pathid)==T)
Pagoda2_scRNA_df_path$path[80]
Pagoda2_scRNA_df_path$pathid[80] <- "General function prediction only"
save(Pagoda2_scRNA_df_path,file = "outdata/Pagoda2_scRNA_df_path.Rdata")

#K-M Survival Analysis
library(survival)
library(survminer)
load("D:/1学习/eye/New_project/scRNAseq_proc/randomForest/TCGA_pathway.Rdata")
surdata<-read.table("data/TCGA-UVM.survival.tsv",sep = "\t",header = T)

sursa<-gsub("-",".",surdata$sample)
Pvalue<-NULL
for(i in 1:nrow(TCGA_pathway)){
  print(i)
  group<-rep("nn",80)
  
  x<-as.numeric(TCGA_pathway[i,])
  zh<-median(x)
  h_g<-colnames(TCGA_pathway)[which(x>=zh)]
  l_g<-colnames(TCGA_pathway)[which(x<zh)]
  
  pp<-match(h_g,sursa)
  group[pp]<-"high"
  pp<-match(l_g,sursa)
  group[pp]<-"low"
  
  surdata$group<-group
  sfit <- survdiff(Surv(OS.time, OS)~group, data=surdata)
  #ggsurvplot(survfit(Surv(OS.time, OS)~group, data=surdata),pval=T)
  Pvalue[i] <- 1 - pchisq(sfit$chisq, length(sfit$n) -1)
}
Path_survival_result<-data.frame(Genes=row.names(TCGA_pathway),Pvalue,FDR=p.adjust(Pvalue,method = "BH"))
save(Path_survival_result,file = "outdata/Path_survival_result.Rdata")

pathid <- NULL
for(i in 1:nrow(Path_survival_result)){
  a <- unlist(strsplit(Path_survival_result$Genes[i]," \\["))[2]
  a <- gsub("]","",a)
  pathid <- c(pathid,a)
}
pp <- which(is.na(pathid)==T)
pathid[pp] <- Path_survival_result$Genes[pp]
Path_survival_result$pathid <- pathid
save(Path_survival_result,file = "outdata/Path_survival_result.Rdata")

Int_path_pagoda2_sur<-intersect(Pagoda2_scRNA_df_path$pathid[which(Pagoda2_scRNA_df_path$p_val_adj<0.0001)],
                                     Path_survival_result$pathid[which(Path_survival_result$Pvalue<0.05)])

pp <- match(Int_path_pagoda2_sur,Path_survival_result$pathid)
Int_path_pagoda2_sur_df <- data.frame(Pathid=Int_path_pagoda2_sur,Path_survival_result$Genes[pp])
save(Int_path_pagoda2_sur_df,file = "outdata/Int_path_pagoda2_sur_df.Rdata")


#
load("D:/scRNAseq_proc2/data/scapGNN_pathway.Rdata")
cell_phen <- rep("Poor",ncol(scapGNN_pathway))
pp <- match(infos1[["Scissor_neg"]],colnames(scapGNN_pathway))
cell_phen[pp] <- "Good"
table(cell_phen)
head(cell_phen)

um_pathway_seurat <- CreateSeuratObject(counts = scapGNN_pathway, project = "scRNA", min.cells = 0, min.features = 0)
names(cell_phen) <- colnames(scapGNN_pathway)
um_pathway_seurat@active.ident <- as.factor(cell_phen)

table(cell_phen)
scapGNN_pathway_df <- FindMarkers(um_pathway_seurat, ident.1 = "Poor",only.pos=T,logfc.threshold=0)
scapGNN_pathway_df$path <- row.names(scapGNN_pathway_df)
save(scapGNN_pathway_df,file = "outdata/scapGNN_pathway_df.Rdata")

Int_pathways<-intersect(scapGNN_pathway_df$path[which(scapGNN_pathway_df$p_val_adj<0.0001)],
                                Int_path_pagoda2_sur_df$Path_survival_result.Genes.pp.)
save(Int_pathways,file = "outdata/Int_pathways.Rdata")


#COX####
pp<-match(Int_pathways,row.names(TCGA_pathway))
pm<-TCGA_pathway[pp,]

px<-match(surdata$sample,colnames(pm))
pm<-pm[,px]
head(colnames(pm))
head(surdata$sample)
pm<-t(pm)
head(row.names(pm))
surdata1<-cbind(surdata,pm)
colnames(surdata1)

library(survival)
library("survminer")

#Single variable cox for each pathway
pathways <- colnames(surdata1)[c(6:29)]
colnames(surdata1)[c(6:29)] <- paste("Path",c(1:24),sep = "")
pathwaydf <- data.frame(id=colnames(surdata1)[c(6:29)],pathways=pathways)

covariates <- colnames(surdata1)[c(6:29)]
univ_formulas <- sapply(covariates,function(x) as.formula(paste('Surv(OS.time, OS)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = surdata1)})

univ_results <- lapply(univ_models,
                       function(x){
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR_all <- paste0(HR, " (",
                                          HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR,HR.confint.lower,HR.confint.upper,HR_all, wald.test, p.value)
                         names(res)<-c("beta","HR","HR.confint.lower","HR.confint.upper", "HR (95% CI for HR)", "wald.test",
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
res <-as.data.frame(res)
pp <- match(row.names(res),pathwaydf$id)

res$path <- pathwaydf$pathways[pp]
res$p.value <- as.numeric(res$p.value)
res$HR <- as.numeric(res$HR)
res$HR.confint.lower <- as.numeric(res$HR.confint.lower)
res$HR.confint.upper <- as.numeric(res$HR.confint.upper)
res <-res[which(res$p.value<0.05),]

#Multivariate cox regression
path <- res$path

pp<-match(path,row.names(TCGA_pathway))
pm<-TCGA_pathway[pp,]

px<-match(surdata$sample,colnames(pm))
pm<-pm[,px]
head(colnames(pm))
head(surdata$sample)
pm<-t(pm)
head(row.names(pm))
surdata1<-cbind(surdata,pm)
colnames(surdata1)
surdata1 <- surdata1[,-c(1,3,5)]

fit.cox <- coxph(Surv(OS.time, OS) ~ . , data = surdata1)
a <- summary(fit.cox)

library(broom)
coxphtable <- tidy(fit.cox)

#Multivariate cox regression
path <- coxphtable$term[which(coxphtable$p.value < 0.05)]

for(i in 1:length(path)){
  path[i] <- gsub("`","",path[i])
}


pp<-match(path,row.names(TCGA_pathway))
pm<-TCGA_pathway[pp,]

px<-match(surdata$sample,colnames(pm))
pm<-pm[,px]
head(colnames(pm))
head(surdata$sample)
pm<-t(pm)
head(row.names(pm))
surdata1<-cbind(surdata,pm)
colnames(surdata1)
surdata1 <- surdata1[,-c(1,3,5)]

fit.cox <- coxph(Surv(OS.time, OS) ~ . , data = surdata1)
a <- summary(fit.cox)

library(broom)
coxphtable2 <- tidy(fit.cox) 

for(i in 1:nrow(coxphtable2)){
  coxphtable2$term[i] <- gsub("`","",coxphtable2$term[i])
}

Pathway12 <- coxphtable2$term
save(Pathway12,file = "outdata/Pathway12.Rdata")

#Five replicate 10-fold cross-validation####
pp <- match(Pathway12,row.names(TCGA_pathway))
TCGA_pathway1 <- TCGA_pathway[pp,]

TCGA_pathway1 <- data.frame(t(TCGA_pathway1))
pp <- match(row.names(TCGA_pathway1),surdata$sample)
TCGA_pathway1$group <- surdata$sur_phen[pp]
TCGA_pathway1[1:3,c(1,ncol(TCGA_pathway1))]

TCGA_pathway1$group <- as.factor(TCGA_pathway1$group)

library(ggplot2)
library(splines)

data_5cv <- replicate(5, rfcv(TCGA_pathway1[-ncol(TCGA_pathway1)], TCGA_pathway1$group, cv.fold = 10,step = 1.5), simplify = FALSE)

data_5cv1<- data.frame(sapply(data_5cv, '[[', 'error.cv'))
data_5cv1$otus <- rownames(data_5cv1)
data_5cv1 <- reshape2::melt(data_5cv1, id = 'otus')
data_5cv1$otus <- as.numeric(as.character(data_5cv1$otus))

p <- ggplot(data_5cv1, aes(otus, value)) +
  geom_smooth(se = FALSE,  method = 'glm', formula = y~ns(x, 6),size=2,color="#81B5D8") +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'),
        axis.text=element_text(size=28,face="bold"),axis.title=element_text(size=25,face="bold")) +
  labs(title = '',x = NULL, y = NULL)

p
p + geom_vline(xintercept = 8,size =1.5, colour="#D25B45", linetype="dashed")

UMpathways <-row.names(Path_importance)[1:8]
