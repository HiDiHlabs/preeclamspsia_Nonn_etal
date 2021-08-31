#conda activate trajectory_analysis_R 
#R CMD BATCH all_INT.R &
library(Seurat)


libs<-read.table("/data/analysis/preeclampsia_2019/analysis/count_matrix/all_libs_filtered.txt", header=T)
dio<-read.table("/data/analysis/preeclampsia_2019/analysis/count_matrix/h5_filtered_feature_matrix.txt")[1]
names(dio)<-"path"
libs$path<-dio$path


grep("*.*",libs$group, value =T)       
all <- libs[libs$group %in% c(grep("*.*",libs$group, value=T)), ]    
cond.data <- vector("list", length(all$group))
for(i in 1:length(all$group)){
  print(paste0("read ", all$path[i]))
  cond.data[[i]]<-Read10X_h5(paste0(all$path[i]), use.names = TRUE, unique.features = TRUE)
  print(paste0("done with ", all$path[i]))
}


setupMultiobject<-function(x,y,z) {
  cond <- CreateSeuratObject(counts = z, project = y, min.cells = 3)
  cond$stim <- x
  cond <- subset(cond, subset = nFeature_RNA > 1)  
  cond <- NormalizeData(cond, verbose = T)
  cond <- FindVariableFeatures(cond, selection.method = "vst", nfeatures = 2000)
  cond <- PercentageFeatureSet(cond, assay = "RNA", col.name ="percent.mt" , pattern = "^MT-")
  cond
}



cond <- vector("list", length(all$group))
for(i in 1:length(all$group)){
  print(paste0("read ", all$library[i]))
  cond[[i]]<-setupMultiobject(all$library[i], "All", cond.data[[i]])
    print(paste0("edone ", all$library[i]))

}
remove(cond.data)

immune.anchors <- FindIntegrationAnchors(object.list = cond, dims = 1:20)
remove(cond)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)
remove(immune.anchors)

DefaultAssay(immune.combined) <- "integrated"
immune.combined<-CellCycleScoring(immune.combined,   g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes, set.ident = T)
# Run the standard workflow for visualization and clustering
#immune.combined <- ScaleData(immune.combined, verbose = FALSE, vars.to.regress = c("S.Score", "G2M.Score") )
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunTSNE(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

saveRDS(immune.combined, "/data/analysis/preeclampsia_2019/analysis/count_matrix/reproduced_seurat_cornelius/placenta_reproduced.rds")