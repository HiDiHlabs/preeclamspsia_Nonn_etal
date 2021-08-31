Sys.setenv(LANG="en_US.UTF-8")
Sys.setenv(LC_ALL="en_US.UTF-8")

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(patchwork))

#From Cellranger to initial Seurat object: Dr.Cornelius Fischer (BIH/MDC Berlin)
#The seurat object from scratch (using raw_h5_matrix) from Cellranger was created using: create_seurat_object.R (needs sample information & full path to the h5 files)
seuratObj  <- readRDS("/data/analysis/preeclampsia_2019/analysis/images/Annotated_placenta_object.rds") #renamed

seuratObj

DefaultAssay(seuratObj) <- "RNA"
seuratObj #Initial object contains 101067 cells. 

#Read oldest cell type annotation:CF & ON/02.12.20 
#unused for the study
clusters <- read.csv("celltype_20201202.csv")
head(clusters)
colnames(clusters) <- c("barcode", "cell_type_new")
barcodes <- strsplit(as.character(clusters$barcode),  "-", fixed = TRUE)
clusters$barcode <- paste0(unlist(lapply(barcodes, `[[`, 1)), "-1_", unlist(lapply(barcodes, `[[`, 2)))

suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(zinbwave))

counts <- GetAssayData(object = seuratObj, slot = "counts")
metadata= seuratObj@meta.data
metadata$barcode <- colnames(counts)
metadata <- join(metadata, clusters, by = "barcode")

head(seuratObj@meta.data)[1:6] #check the top few lines 

seuratObj[["percent.mt"]] <- PercentageFeatureSet(seuratObj, pattern = "^MT-")
VlnPlot(seuratObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size=0)

Idents(object= seuratObj) <- "orig.ident"
VlnPlot(seuratObj, features = "percent.mt", pt.size=0)
plot1 <- FeatureScatter(seuratObj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seuratObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

#Remove probable outliers: second filtering step: working Seurat object 
placenta_subset <- subset(seuratObj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
placenta_subset

#Previous filtering (directly after Cellranger) is already described in R-script: create_seurat_object

saveRDS(placenta_subset, "/data/analysis/preeclampsia_2019/analysis/updated_seurat_clustering/placenta_seuratobj_231220.rds")
#Input for step2 (integration & basic clustering)

#To check the results for integration between decidua & villi and comparative cell typing, please see: Reclustering_Seurat_Step1.ipynb
#unsed for manuscript. 
