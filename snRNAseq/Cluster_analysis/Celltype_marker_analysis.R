#conda activate trajectory_analysis_R
#Run as bsub script: bsub -n 16 -R "span[ptile=2]" R CMD BATCH Celltype_marker_analysis.R 


Sys.setenv(LANG="en_US.UTF-8")
Sys.setenv(LC_ALL="en_US.UTF-8")

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(tibble))

library(future)
plan("multiprocess", workers = 15)
data <- readRDS(file = "/data/analysis/preeclampsia_2019/analysis/images/updated_seurat_clustering/Placenta_seuratobj_merged_final_2021.rds")

DefaultAssay(data) <- "RNA"

#Set idents to the lastest cluster annotations
Idents(object= data) <- 'cell_type_semifinal_v2'
table(Idents(data))

#Run Negative Binomial model:
markers_NB <- FindAllMarkers(data, logfc.threshold = 0.25,
  only.pos = TRUE, min.pct = 0.25,
  test.use= "negbinom",
  latent.vars= c("nCount_RNA", "nFeature_RNA", "percent.mt"))

write.csv(markers_NB, file= "Signatures_NegBinom_1304_covars_corrected.csv")

#Run Negative Binomial model: add the categorical covariate ~disease. 
markers_NB <- FindAllMarkers(data, logfc.threshold = 0.25,
  only.pos = TRUE, min.pct = 0.25,
  test.use= "negbinom",
  latent.vars= c("nCount_RNA", "nFeature_RNA", "percent.mt", "disease"))

write.csv(markers_NB, file= "Signatures_NegBinom_1304_disease_covars.csv")

#Filter the top 30 genes per cluster:
NB_filtered= markers_NB %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
write.csv(NB_filtered, file= "Signatures_NegBinom_1304_disease_covars_filtered_top30.csv")

#Verify the markers additionally with Logistic Regression:
markers_LR <- FindAllMarkers(data, logfc.threshold = 0.25,
  only.pos = TRUE, min.pct = 0.25,
  test.use= "LR",
  latent.vars= c("nCount_RNA", "nFeature_RNA", "percent.mt", "disease"))

write.csv(markers_LR, file= "Signatures_LogReg_1304_disease_covars.csv")

#Filter the top 30 genes per cluster: LR 
LR_filtered= markers_LR %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
write.csv(LR_filtered, file= "Signatures_LogReg_1304_disease_covars_filtered_top30.csv")

#Correcting for nCount_RNA, nFeature_RNA, percent.MT (and, ~disease) doesn't significantly effect top cluster markers. To verify, run the same function without the confounding/latent variables. 