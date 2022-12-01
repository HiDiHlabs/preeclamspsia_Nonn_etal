Sys.setenv(LANG="en_US.UTF-8")
Sys.setenv(LC_ALL="en_US.UTF-8")
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(tibble))

#The H5AD is converted to RDS using Seurat Disk. 
data <- readRDS(file = "/data/analysis/preeclampsia_2019/placenta_atlas_2022/villi_scVI_harmonization/villi_seurat/SP014_SP082_SP136_placenta_annotations_250322_noNaNs_KS.rds")
data <- NormalizeData(object = data)

DefaultAssay(data) <- "RNA"

#Set idents to the lastest "cluster annotations" having the polished cell annotations. 
Idents(object= data) <- 'leiden_subclusters_refined02'
table(Idents(data))


nCount_RNA = colSums(x = data, slot = 'counts')
nFeature_RNA = colSums(x = GetAssayData(object = data, slot = "counts") > 0)
data <- AddMetaData(data, nFeature_RNA, col.name = "nFeature_RNA")
data <- AddMetaData(data, nCount_RNA, col.name= "nCount_RNA")

#Calculate percent.MT inside Seurat:
mt.genes <- rownames(data)[grep("^MT-",rownames(data))]
C <-GetAssayData(object = data, slot = "counts")
percent.mito <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
data <- AddMetaData(data, percent.mito, col.name = "percent.mt")


#Run Logistic Regression for all cell types:
#Correct for latent confounders such as "nCount_RNA", "nFeature_RNA", "percent.mt", "disease", "library" to prevent technical effects. 
markers_LR <- FindAllMarkers(data, logfc.threshold = 0.25,
  only.pos = TRUE, min.pct = 0.25,
  test.use= "LR",
  latent.vars= c("nCount_RNA", "nFeature_RNA", "percent.mt", "disease", "library"))

write.csv(markers_LR, file= "Placenta_markers_LogisticRegression_290322_covars.csv")
#This LR table is reported under Supplementary Table-8 (for cell markers)

