Sys.setenv(LANG="en_US.UTF-8")
Sys.setenv(LC_ALL="en_US.UTF-8")

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(tibble))

#Read the decidua RDS file. Converted from anndata H5AD via SeuratDisk. 
data <- readRDS(file = "/data/analysis/preeclampsia_2019/placenta_atlas_2022/decidua_scVI_harmonization/decidua_seurat/Raw_decidua_OD_140322.rds")
data <- NormalizeData(object = data)

DefaultAssay(data) <- "RNA"

#Set idents to the lastest "cluster annotations"
Idents(object= data) <- 'celltype_v5'
table(Idents(data))

#Calculate nCount_RNA & nFeature_RNA inside Seurat:
nCount_RNA = colSums(x = data, slot = 'counts')
nFeature_RNA = colSums(x = GetAssayData(object = data, slot = "counts") > 0)
data <- AddMetaData(data, nCount_RNA, col.name= "nCount_RNA")
data <- AddMetaData(data, nFeature_RNA, col.name = "nFeature_RNA")

#Calculate percent.MT inside Seurat:
mt.genes <- rownames(data)[grep("^MT-",rownames(data))]
C <-GetAssayData(object = data, slot = "counts")
percent.mito <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
data <- AddMetaData(data, percent.mito, col.name = "percent.mt")


#Run Logistic Regression for all cell types:
#DEG(s) are corrected for categorical variables such as ~disease (PE vs controls); library (10X V3 vs 10X V2) and other continuous covars. 
markers_LR <- FindAllMarkers(data, logfc.threshold = 0.25,
  only.pos = TRUE, min.pct = 0.25,
  test.use= "LR",
  latent.vars= c("nCount_RNA", "nFeature_RNA", "percent.mt", "disease", "library"))

write.csv(markers_LR, file= "Decidua_markers_LogisticRegression_280322_covars.csv")

