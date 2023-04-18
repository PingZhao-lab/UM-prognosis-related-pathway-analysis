#GSVA####
load_gs_data = function(gmt_file_path){
  # load pathway gene sets' gmt file
  # file format:
  #            first index: pathway's name/ID
  #            second index: pathway's url or others, it dosen't matter
  #            third to all: gene symbols in pathway
  #            sep: \t
  # return a list
  tmp = readLines(gmt_file_path)
  gsets = list()
  for(i in 1:length(tmp)){
    t = strsplit(tmp[i],'\t')[[1]]
    genes = t[3:length(t)]
    genes = genes[which(genes != "")]
    gsets[[t[1]]] = genes
  }
  return (gsets)
}
gs_list<-load_gs_data("./data/KEGG_human.gmt")

library(GSVA)
TCGA_pathway<-gsva(expr=as.matrix(data_TCGA),gset.idx.list=gs_list,
     min.sz=5)
save(TCGA_pathway,file = "outdata/TCGA_pathway.Rdata")

#Pagoda2####
cal_pagoda2 = function(counts,
                       gSets,
                       trim = 5,
                       n_cores){
  
  
  ### must be counts matrix !!!!!
  
  ### other parameters for knn.error.models
  # min.nonfailed = 5, min.count.threshold = 1,
  # max.model.plots = 50,
  # min.size.entries = 2000, min.fpm = 0, cor.method = "pearson",
  # verbose = 0, fpm.estimate.trim = 0.25, linear.fit = TRUE,
  # local.theta.fit = linear.fit, theta.fit.range = c(0.01, 100),
  # alpha.weight.power = 1/2
  
  ### other parameters for pagoda.varnorm
  # batch = NULL, prior = NULL,
  # fit.genes = NULL, minimize.underdispersion = FALSE,
  # n.cores = detectCores(), n.randomizations = 100, weight.k = 0.9,
  # verbose = 0, weight.df.power = 1, smooth.df = -1,
  # theta.range = c(0.01, 100), gene.length = NULL
  library(pagoda2)
  print(gc())
  nPcs = min(round(ncol(counts)/5),5)
  #counts = apply(counts,2,function(x) {storage.mode(x) = 'integer'; x})
  tryCatch({
    p2 = Pagoda2$new(counts, n.cores = n_cores,log.scale=T)
    print(gc())
    p2$adjustVariance(plot=F)
    print(gc())
    p2$calculatePcaReduction(nPcs = nPcs,use.odgenes=F,fastpath=F)
    print(gc())
    path_names = c()
    env = new.env(parent=globalenv())
    invisible(lapply(1:length(gSets),function(i) {
      genes = intersect(gSets[[i]],rownames(counts))
      name = paste0(names(gSets[i]),i)
      if(length(genes)>3){
        assign(name, genes, envir = env)
        path_names = c(path_names, name)
      }
    }))
    print(gc())
    p2$testPathwayOverdispersion(setenv = env, verbose = T,
                                 recalculate.pca = T,
                                 min.pathway.size = 1)
    print(gc())
    path_names = names(p2$misc$pwpca)
    score = matrix(NA,nrow=length(path_names),ncol=ncol(counts))
    rownames(score) = path_names
    colnames(score) = colnames(counts)
    for(i in 1:length(p2$misc$pwpca)){
      if(!is.null(p2$misc$pwpca[[i]]$xp$scores)){
        score[i,] = as.numeric(p2$misc$pwpca[[i]]$xp$scores)
      }
    }
    print(gc())
    return(score)
  },error = function(e){
    print(e)
    return("error")
  })
}


library(pagoda2)
load("D:/myproject/E_UM_int/outdata/scdata_scissor.Rdata")
library(scde)

scRNA_um_pagoda2<-cal_pagoda2(as.matrix(scdata_scissor@assays[["RNA"]]@counts),gs_list,n_cores=1)
save(scRNA_um_pagoda2,file = "outdata/scRNA_um_pagoda2.Rdata")


#scapGNN####
library(scapGNN)
load("D:/myproject/E_UM_int/outdata/scdata_scissor.Rdata")
library(Seurat)
data_s <- CreateSeuratObject(counts = scdata_scissor@assays[["RNA"]]@counts, project = "data")
data_s <- NormalizeData(data_s, normalization.method = "LogNormalize")

a <- colnames(data_s)
pp <- match(infos1[["Scissor_pos"]],a)
a1 <- a
a1[pp] <- "pos"

pp <- match(infos1[["Scissor_neg"]],a)
a1[pp] <- "neg"
table(a1)
names(a1) <- a

data_s@active.ident <- as.factor(a1)

data_s <- FindVariableFeatures(data_s, selection.method = "vst",nfeatures = 3000)

cluster2.markers <- FindMarkers(data_s, ident.1 = "pos", logfc.threshold = 0)
cluster2.markers$genes <- row.names(cluster2.markers)
row.names(cluster2.markers) <- 1:nrow(cluster2.markers)
genes <- cluster2.markers$genes[1:2000]
pp <- match(genes,row.names(data_s@assays[["RNA"]]@data))
exp1 <- data_s@assays[["RNA"]]@data[pp,]
save(exp1,file = "outdata/exp1.Rdata")

rm(scdata_scissor)
gc()

library(coop)
system.time(
  Prep_data1 <- Preprocessing(exp1,verbose=T)
)
save(Prep_data1,file = "outdata/Prep_data1.Rdata")

library(reticulate)
library(parallel)
system.time(
  ConNetGNN_data <- ConNetGNN(Prep_data1,python.path="C:/Users/dell/.conda/envs/GCN/python.exe")
)
save(ConNetGNN_data,file = "outdata/ConNetGNN_data.Rdata")


library(parallel)
library(compiler)
load("D:/hanxudong/BIB_xd/tp1/outdata/ConNetGNN_data.Rdata")
#scPathway<- cmpfun(scPathway)
system.time(
  scapGNN_pathway<-scPathway(ConNetGNN_data,gmt.path="./data/KEGG_human.gmt",nperm=100,
                       parallel.cores=30,pathway.min = 10)
)
save(scapGNN_pathway,file="outdata/scapGNN_pathway.Rdata")
