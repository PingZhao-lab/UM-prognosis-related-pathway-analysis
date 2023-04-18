library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

list.files("/data1/lab/Hanxudong_data/myProject/ZPEYE_UM/")
f <- c("GSM4147091",    "GSM4147092",    "GSM4147093",
"GSM4147094" ,   "GSM4147095" ,   "GSM4147096" ,   "GSM4147097",
"GSM4147098"  ,    "GSM4147099_1",  "GSM4147100",
"GSM4147101")

GSM4147091 <- Read10X(data.dir = "/data1/lab/Hanxudong_data/myProject/ZPEYE_UM/GSM4147091/")
GSM4147092 <- Read10X(data.dir = "/data1/lab/Hanxudong_data/myProject/ZPEYE_UM/GSM4147092/")
GSM4147093 <- Read10X(data.dir = "/data1/lab/Hanxudong_data/myProject/ZPEYE_UM/GSM4147093/")
GSM4147094 <- Read10X(data.dir = "/data1/lab/Hanxudong_data/myProject/ZPEYE_UM/GSM4147094/")
GSM4147095 <- Read10X(data.dir = "/data1/lab/Hanxudong_data/myProject/ZPEYE_UM/GSM4147095/")
GSM4147096 <- Read10X(data.dir = "/data1/lab/Hanxudong_data/myProject/ZPEYE_UM/GSM4147096/")
GSM4147097 <- Read10X(data.dir = "/data1/lab/Hanxudong_data/myProject/ZPEYE_UM/GSM4147097/")
GSM4147098 <- Read10X(data.dir = "/data1/lab/Hanxudong_data/myProject/ZPEYE_UM/GSM4147098/")
GSM4147099 <- Read10X(data.dir = "/data1/lab/Hanxudong_data/myProject/ZPEYE_UM/GSM4147099_1/")
GSM4147100 <- Read10X(data.dir = "/data1/lab/Hanxudong_data/myProject/ZPEYE_UM/GSM4147100/")
GSM4147101 <- Read10X(data.dir = "/data1/lab/Hanxudong_data/myProject/ZPEYE_UM/GSM4147101/")

GSM4147091 <- CreateSeuratObject(counts = GSM4147091, project = "GSM4147091", min.cells = 3, min.features = 100)
GSM4147092 <- CreateSeuratObject(counts = GSM4147092, project = "GSM4147092", min.cells = 3, min.features = 100)
GSM4147093 <- CreateSeuratObject(counts = GSM4147093, project = "GSM4147093", min.cells = 3, min.features = 100)
GSM4147094 <- CreateSeuratObject(counts = GSM4147094, project = "GSM4147094", min.cells = 3, min.features = 100)
GSM4147095 <- CreateSeuratObject(counts = GSM4147095, project = "GSM4147095", min.cells = 3, min.features = 100)
GSM4147096 <- CreateSeuratObject(counts = GSM4147096, project = "GSM4147096", min.cells = 3, min.features = 100)
GSM4147097 <- CreateSeuratObject(counts = GSM4147097, project = "GSM4147097", min.cells = 3, min.features = 100)
GSM4147098 <- CreateSeuratObject(counts = GSM4147098, project = "GSM4147098", min.cells = 3, min.features = 100)
GSM4147099 <- CreateSeuratObject(counts = GSM4147099, project = "GSM4147099", min.cells = 3, min.features = 100)
GSM4147100 <- CreateSeuratObject(counts = GSM4147100, project = "GSM4147100", min.cells = 3, min.features = 100)
GSM4147101 <- CreateSeuratObject(counts = GSM4147101, project = "GSM4147101", min.cells = 3, min.features = 100)


Eye_combine <- merge(GSM4147091, y = c(GSM4147092,GSM4147093,GSM4147094,GSM4147095,GSM4147096,GSM4147097,GSM4147098,GSM4147099,GSM4147100,GSM4147101),
                  add.cell.ids = c("GSM4147091",    "GSM4147092",    "GSM4147093",
                                   "GSM4147094" ,   "GSM4147095" ,   "GSM4147096" ,   "GSM4147097",
                                   "GSM4147098"  ,    "GSM4147099",  "GSM4147100",
                                   "GSM4147101"), project = "Combine")

a <- Eye_combine@active.ident
a <- rep("um",length(a))
names(a) <- names(Eye_combine@active.ident)
Eye_combine@active.ident <- as.factor(a)
Eye_combine[["percent.mt"]] <- PercentageFeatureSet(Eye_combine, pattern = "^MT-")
VlnPlot(Eye_combine, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Eye_combine <- subset(Eye_combine, subset = nFeature_RNA > 100 & nFeature_RNA < 8000 & percent.mt < 10 & nCount_RNA > 400)

Eye_combine <- NormalizeData(Eye_combine)
Eye_combine <- FindVariableFeatures(Eye_combine, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(Eye_combine)
Eye_combine <- ScaleData(Eye_combine, features = all.genes)

Eye_combine <- RunPCA(Eye_combine, features = VariableFeatures(object = Eye_combine))

Eye_combine <- FindNeighbors(Eye_combine, dims = 1:20)
Eye_combine <- FindClusters(Eye_combine, resolution = 3)

Eye_combine <- RunTSNE(Eye_combine, dims = 1:20)
DimPlot(Eye_combine, reduction = "tsne")

save(Eye_combine,file = "Eye_combine.Rdata")

png("tsne.png" ,width = 1600, height = 1200)
DimPlot(Eye_combine, reduction = "tsne",label=T)
dev.off()


png("Tumor_xtq.png" ,width = 3000, height = 780)
VlnPlot(Eye_combine, features = c("MLANA", "MITF", "DCT"))
dev.off()

png("Tumor_PRAME.png" ,width = 1600, height = 780)
VlnPlot(Eye_combine, features = "PRAME")
dev.off()

DefaultAssay(Eye_combine) <- "RNA"
png("Tumor_marker_sd.png" ,width = 1600, height = 1600)
FeaturePlot(Eye_combine, features = c("MLANA", "MITF", "DCT", "PRAME", "GEP"))
dev.off()

png("Tumor_marker_GARS.png" ,width = 600, height = 600)
FeaturePlot(Eye_combine, features = "GARS")
dev.off()

load("Eye_combine.Rdata")
new.cluster.ids <- c(0:64)
Tid <- c(43,3,2,4,42,23,0,1,19,15,13,11,18,38,30,14,20,22,35,12,10,9,6,8,5,7,33,37,50,57,54)
pp <- match(Tid,new.cluster.ids)
new.cluster.ids[pp] <- "Tumor cells"
names(new.cluster.ids) <- levels(Eye_combine)
Eye_combine <- RenameIdents(Eye_combine, new.cluster.ids)
DimPlot(Eye_combine, reduction = "tsne",label=T)

save(Eye_combine,file = "Eye_combine_ann.Rdata")

Tumor_UM_cells <- subset(Eye_combine, idents = "Tumor cells")
DimPlot(Tumor_UM_cells, reduction = "tsne",label=T)
save(Tumor_UM_cells,file = "Tumor_UM_cells.Rdata")


