Sys.setenv(LANG="en_US.UTF-8")
Sys.setenv(LC_ALL="en_US.UTF-8")

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(patchwork))

#Read the Seurat object to perform reclustering:
#In previous case, two tissues were integrated onto each other to perform comparative cell-typing. 

data <- readRDS(file= "/data/analysis/preeclampsia_2019/analysis/images/updated_seurat_clustering/placenta_seuratobj_231220.rds")

data

#Old cell-type (don't refer to this): 
Idents(object= data) <- 'integrated_snn_res.1'

table(Idents(data))

#Old annotation/unused: 
current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                        21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33)

new.cluster.ids <- c("SCT_classical_0", "SCT_classical_1", "SCT_growth_2", "SCT_endodysfun_3", "Mac_classical_4", "EpC_classical_5", "VCT_6", "SCT_lowhormone_7", "SCT_dysfinflamm_8", "dNK_resting_9", "dNK_resting_10", "Mac_activated_11","Treg_12", "LEC_13", "VEC_14", "Mac_classical_15", "DSC_MSC_16", "ILC_17", "dNK_activated_18", "FB_19","Tropho/MSC_20", "Cluster_21", "VCTprol_22", "EB_23", "PlasmaCell_24", "EpC_classical_25", "EVT_26", "dNK_prol_27", "TSC_28", "Mac_classical_29", "Mac_polarised_30", "EpC_inflamm_31", "CD141_DC_32", "EpC_classical_33")



#data@ident <- plyr::mapvalues(x = data@integrated_snn_res.1, from = current.cluster.ids, to = new.cluster.ids)
names(x = new.cluster.ids) <- levels(x = data)
data <- RenameIdents(object = data, new.cluster.ids)
table(Idents(data))

#Change to RNA slot:: 
#We need to separate "Decidua" & "Villi" individually and then, re-integrate the samples within tissue. 
DefaultAssay(data) <- "RNA"
#data
Idents(object= data) <- 'tissue_time'
table(Idents(data))

#Set idents to "tissue": 
Idents(object= data) <- 'tissue'
table(Idents(data))

#At first, we need to subset all the "Villi" cells for the subsequent integration: 

subset= "Villi"
seurat_villi = subset(data, idents = subset) 
seurat_villi

#Recheck if the active assay is set to "RNA slot:"
DefaultAssay(seurat_villi) <- "RNA"


#Change ident to "group" & look at the cell-composition:
Idents(object= seurat_villi) <- 'group'
table(Idents(seurat_villi))

#Proceed with Log-normalization in Seurat:
seurat_villi.CCA <- SplitObject(seurat_villi, split.by = "group")

#Re-normalize & find highly variable features using VST (recommended): 
#Alternatively, perform SCTransform:  
seurat_villi.CCA <- lapply(X = seurat_villi.CCA, FUN = function(x) {
    x <- NormalizeData(x)
    #x <- SCTransform(x, verbose = FALSE, variable.features.n = 3000)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})

seurat_villi.CCA

DimPlot(seurat_villi,
        reduction = "pca",
        group.by= "disease",
        split.by = "disease")

DimPlot(seurat_villi,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")

#The following function will call SelectIntegrationFeatures to select the provided number of features to be used in anchor finding. 
villi.anchors <- FindIntegrationAnchors(object.list = seurat_villi.CCA, dims = 1:30, reduction= "cca", anchor.features = 3000, normalization.method= "LogNormalize")


#Integrate within villi: 
villi.combined <- IntegrateData(anchorset = villi.anchors, dims = 1:30)
DefaultAssay(villi.combined) <- "integrated" #set assay to integrated prior to merging. 
villi.combined

#Save the decidua.combined for future use:
#The above step can be also done faster by submitting as a bsub R-script. 
saveRDS(villi.combined, file = "./updated_seurat_clustering/Villi_seuratobj_final_15012021.rds")

#Rerun the above steps for integration within the Decidua: 

subset= "Decidua"
#Subset only Decidua cells:
seurat_decidua= subset(data, idents = subset) 
seurat_decidua

#Proceed with Log-normalization in Seurat:
seurat_decidua.CCA <- SplitObject(seurat_decidua, split.by = "group")


#Re-normalize & find highly variable features using VST (recommended): 
#Alternatively, perform sctransform.  
seurat_decidua.CCA <- lapply(X = seurat_decidua.CCA, FUN = function(x) {
    x <- NormalizeData(x)
    #x <- SCTransform(x, verbose = FALSE, variable.features.n = 3000)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})

seurat_decidua.CCA

#The following will call SelectIntegrationFeatures to select the provided number of features to be used in anchor finding. 
decidua.anchors <- FindIntegrationAnchors(object.list = seurat_decidua.CCA, dims = 1:30, reduction= "cca", anchor.features = 3000, normalization.method= "LogNormalize")


#Integrate into a combined Decidua dataset: 
decidua.combined <- IntegrateData(anchorset = decidua.anchors, dims = 1:30)
DefaultAssay(decidua.combined) <- "integrated" #Set assay to integrated. 
decidua.combined

#Save the decidua.combined for future use
saveRDS(decidua.combined, file = "./updated_seurat_clustering/Decidua_seuratobj_final2_15012021.rds")

#Add "VL" & "DC" prefix to separate the villi & decidua cells clearly post merging. 
placenta.combined <- merge(villi.combined, y = decidua.combined, add.cell.ids = c("VL", "DC"), project = "Placenta")
placenta.combined
placenta.combined@assays

#notice the cell names now have an added identifier as described above: 
head(colnames(placenta.combined))
head(placenta.combined@meta.data)
table(placenta.combined$orig.ident)

#Extract the highly variable features from individual tissue objects: 
villi_hvg_features = VariableFeatures(object = villi.combined)
decidua_hvg_features = VariableFeatures(object = decidua.combined)

length(decidua_hvg_features)
length(villi_hvg_features)

#Find intersection of HVG (highly variable features) from both datasets: 
common.features <- intersect(villi_hvg_features, decidua_hvg_features)
length(x = common.features)

#Union of HVG from both datasets:: essentially, 4311 genes are stored in the "integrated" slot after merging. 
union_hvg <- union(villi_hvg_features, decidua_hvg_features)
length(x = union_hvg)

#Run the standard Seurat workflow for visualization and clustering: 
#Step1: Scale the merged data. 
placenta.combined.new <- ScaleData(placenta.combined, verbose = FALSE, features= union_hvg)

#Features used by PCA must be present in the "scaled.data"
#Default: features = VariableFeatures(object = placenta.combined.new)

placenta.combined.new <- RunPCA(placenta.combined.new, verbose = TRUE, features= union_hvg)
print(placenta.combined.new[["pca"]], dims = 11:16, nfeatures = 10)
ElbowPlot(placenta.combined.new)

#Randomly permute a subset of the data (1% by default) and rerun PCA, constructing a ‘null distribution’ of feature scores, and repeat this procedure. We identify ‘significant’ PCs as those who have a strong enrichment of low p-value features.
placenta.combined <- JackStraw(placenta.combined, num.replicate = 100)
placenta.combined <- ScoreJackStraw(placenta.combined, dims = 1:40)

#The JackStrawPlot() function provides a visualization tool for comparing the distribution of p-values for each PC with a uniform distribution (dashed line). ‘Significant’ PCs will show a strong enrichment of features with low p-values (solid curve above the dashed line).
JackStrawPlot(placenta.combined, dims = 1:20)
JackStrawPlot(placenta.combined, dims = 21:40)

#Viz Top 15 PC genes by loadings to explore the major (marker) genes encoding each PC: 
#If a PC seems to be confounded by contradictory markers (for example, contrasting lineages that shouldn't mix up), please exclude it from downstream analysis. 
DefaultAssay(placenta.combined.new) <- "integrated"

pdf("./Merged_reclustering_final/VizPCATopGenes_Merged01_Final.pdf", w=15, h= 20, paper='special') 
VizDimLoadings(placenta.combined.new, dims = 1:15, reduction = "pca")
dev.off()

#Viz Top 16 to 30 PC genes by loadings.
#DefaultAssay(placenta.combined.new) <- "integrated"

pdf("./Merged_reclustering_final/VizPCATopGenes_Merged02_Final.pdf", w=15, h= 20, paper='special') 
VizDimLoadings(placenta.combined.new, dims = 16:30, reduction = "pca")
dev.off()

#Viz top 31:40 PC genes by loadings.
#DefaultAssay(placenta.combined.new) <- "integrated"

pdf("./Merged_reclustering_final/VizPCATopGenes_Merged03_Final.pdf", w=15, h= 20, paper='special') 
VizDimLoadings(placenta.combined.new, dims = 31:40, reduction = "pca")
dev.off()

#placenta.combined.new
#features = VariableFeatures(object = placenta.combined.new)
#placenta.combined.new@assays$integrated@var.features
#placenta.combined.new

#Save PC features as heatmaps: PC1-10. 
pdf("./Merged_reclustering_final/Heatmap01_PCATopGenes_Merged01_Final.pdf", w=20, h= 20, paper='special') 
DimHeatmap(placenta.combined.new, dims = 1:10, cells = 800, balanced = TRUE)
dev.off()

#PC: 11-20. 
pdf("./Merged_reclustering_final/Heatmap02_PCATopGenes_Merged02_Final.pdf", w=20, h= 20, paper='special') 
DimHeatmap(placenta.combined.new, dims = 11:20, cells = 800, balanced = TRUE)
dev.off()

#PC: 21-30. 
pdf("./Merged_reclustering_final/Heatmap03_PCATopGenes_Merged03_Final.pdf", w=20, h= 20, paper='special') 
DimHeatmap(placenta.combined.new, dims = 21:30, cells = 800, balanced = TRUE)
dev.off()

pdf("./Merged_reclustering_final/Heatmap04_PCATopGenes_Merged04_Final.pdf", w=20, h= 20, paper='special') 
DimHeatmap(placenta.combined.new, dims = 31:40, cells = 800, balanced = TRUE)
dev.off()

#Calculating shared nearest neighbor or SNN graph: 
placenta.combined.new <- FindNeighbors(placenta.combined.new, reduction = "pca", dims = 1:35)
placenta.combined.new <- RunUMAP(placenta.combined.new, reduction = "pca", dims = 1:35, umap.method= "umap-learn")
#placenta.combined.new

#Rerun clustering: using the integration within tissue approach 
placenta.combined.new <- FindClusters(placenta.combined.new, resolution = 1, algorithm = 2)
placenta.combined.new <- StashIdent(object = placenta.combined.new, save.name = "Merged_cell_type_PC")
Idents(object= placenta.combined.new) <- 'Merged_cell_type_PC'
table(Idents(placenta.combined.new)) #Initial clusters; the cluster with 64 cells is a possibly doublet and lack specific/robust markers. Excluded from downstream analysis. 



#Save & rewrite the Seurat object as RDS file: further downstream analysis are based on this 
saveRDS(placenta.combined.new, file = "./updated_seurat_clustering/Merged_placenta_seuratobj_final_17012021.rds")


pdf("./Merged_reclustering_final/UMAP_mergedData_reclustering_Version01_Label.pdf", w=12, h=9)
UMAPPlot(object = placenta.combined.new, cols= "polychrome", label = TRUE)#recomputed later using a fixed seed to ensure reproducibility. 
dev.off()



union_hvg_new= as.data.frame(union_hvg)
#Write the HVG genes used for clustering:
write.csv(union_hvg_new, file= "./Merged_reclustering_final/Merged_HVG_clustering_4311.csv")



#Explore how the first PC dimension is distributed across clusters: also, see how it effects the Controls vs PE cells. 

PC_cluster <- VlnPlot(placenta.combined.new, features = "PC_4", sort= TRUE, split.by= "disease", pt.size = 0)
pdf("./Merged_reclustering_final/PCA_exploratory_merged/PC4_merged_VlnPlot_170120.pdf", w=15, h= 10)
PC_cluster
dev.off()






