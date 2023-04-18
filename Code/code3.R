library(Scissor)
library(preprocessCore)

load("D:/myproject/E_UM_int/data/data_TCGA.Rdata")
load("D:/myproject/E_UM_int/data/phenotype.Rdata")
load("D:/myproject/E_UM_int/data/scdata.Rdata")

infos1<- Scissor(as.matrix(data_TCGA), scdata, phenotype,alpha =0.05,
                 family = "cox", Save_file = 'Scissor_LUAD_survival.RData')
save(infos1,file = "infos1.Rdata") #Scissor_neg对应的是生存好的细胞


Scissor_select <- rep(0, ncol(scdata))
names(Scissor_select) <- colnames(scdata)
Scissor_select[infos1$Scissor_pos] <- 1
Scissor_select[infos1$Scissor_neg] <- 2
scdata <- AddMetaData(scdata, metadata = Scissor_select, col.name = "scissor")
DimPlot(scdata, reduction = 'tsne', group.by = 'scissor', cols = c('grey','indianred1','royalblue'), pt.size = 1.2, order = c(2,1))

cellid <- row.names(scdata@meta.data)[which(scdata@meta.data$scissor!=0)]
scdata_scissor <- subset(scdata,cells=cellid)
save(scdata_scissor,file = "scdata_scissor.Rdata") 

VlnPlot(scdata_scissor, features = c("GARS"),group.by="scissor",slot = "counts", log = TRUE)
