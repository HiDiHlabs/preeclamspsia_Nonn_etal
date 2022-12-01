Sys.setenv(LANG="en_US.UTF-8")
Sys.setenv(LC_ALL="en_US.UTF-8")

#Load libraries:
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(tibble))

#Part-I: Decidua PE vs term control differentially expressed signatures:
#Read the decidua RDS file. Converted from anndata H5AD via SeuratDisk.
#Note that, after Sept'22, the file path on Eils-HPC changed to: /dh-projects/preeclampsia_2022/analysis/placenta_atlas_2022/decidua_scVI_harmonization/decidua_seurat/Raw_decidua_OD_140322.rds
#Upon the Zenodo release, one can also simply use the decidua RDS file for analysis.
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

condition_oi = "late_preterm" #Condition of interest (disease/PE)
condition_reference = "late_term" #reference condition (late term)

#Subset each cell-type/state & run Logistic Regression for eoPE vs term-controls:
#Subset dNK1:
receiver = "dNK1"
nkc_obj_receiver= subset(data, idents = receiver)
nkc_obj_receiver = SetIdent(nkc_obj_receiver, value = nkc_obj_receiver[["time"]])
table(Idents(nkc_obj_receiver))

#Compute a score of preterm/labor specific features to adjust the confounding factors from analysis:
#We've particularly used features explained by Pique-Regi etl.al 2019(doi:10.7554/eLife.52004)
#Selected features were discussed with Dr. Olivia Nonn (04.04.2022)
#Since the maternal decidua is female, correction for Y-genes isn't required. 
nkc_obj_receiver <- PercentageFeatureSet(object = nkc_obj_receiver, features = c("SLC30A2", "GKN1", "SERPINE2", "GALNT6", "MYOZ1", "ARHGDIB", "GDPD3", "DSG2",
"IGFBP2", "HLA-DRB5", "HLA-DPB1"), col.name = 'NK_preterm_features')

#Find dNK1 DEGs after correcting for confounding factors (encoded by latent variables as implemented inside Seurat):
nk1_list = FindMarkers(object = nkc_obj_receiver, ident.1 = condition_oi,
            ident.2 = condition_reference, min.pct = 0.10,logfc.threshold=0.25,
latent.vars= c("nCount_RNA", "nFeature_RNA", "XIST", "NK_preterm_features","MALAT1", "percent.mt"),
test.use= "LR") %>% rownames_to_column("gene")
write.csv(nk1_list, file= "./PE_decidua_markers_logreg/dNK1_PE_vs_lateC_preterm_corrected_040422.csv")

#Subset dNK2:
receiver = "dNK2"
nkc_obj_receiver= subset(data, idents = receiver)
nkc_obj_receiver = SetIdent(nkc_obj_receiver, value = nkc_obj_receiver[["time"]])
table(Idents(nkc_obj_receiver))

#Compute a score of preterm/labor specific features to adjust the confounding factors from analysis (same as dNK1):
nkc_obj_receiver <- PercentageFeatureSet(object = nkc_obj_receiver, features = c("SLC30A2", "GKN1", "SERPINE2", "GALNT6", "MYOZ1", "ARHGDIB", "GDPD3", "DSG2",
"IGFBP2", "HLA-DRB5", "HLA-DPB1"), col.name = 'NK_preterm_features')

#Find dNK2 DEGs after correcting for confounding factors (encoded by latent variables as implemented inside Seurat):
nk2_list = FindMarkers(object = nkc_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10,
                       logfc.threshold=0.25, latent.vars= c("nCount_RNA", "nFeature_RNA", "XIST", "NK_preterm_features", "MALAT1", "percent.mt"),
                       test.use= "LR") %>% rownames_to_column("gene")
write.csv(nk2_list, file= "./PE_decidua_markers_logreg/dNK2_PE_vs_lateC_preterm_corrected_040422.csv")

#For dNKprol, the sparsity of cells per group & composition imbalance doesn't allow for meaningful p-value calculation & hence, DEG analysis.

#Subset dMAC1:
receiver = "dMAC1"
seurat_obj_receiver= subset(data, idents = receiver)
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["time"]])
seurat_obj_receiver
table(Idents(seurat_obj_receiver))

mac_preterm_features= c("CLDN5", "CD3E", "LILRB5", "FILIP1L", "RORA", "MMP19", "CD163", "SOD2", "BIRC3","SDC4", "ARL4C", "KDM6B", "ANGPTL4", "B4GALT1", "JUND", "NFKB1", "PDE4B", "B4GALT1", "JUND", "NFKB1", "PDE4B",
"MAFF", "NAMPT", "CD55", "C1QA", "CFLAR", "LMNA", "SAMSN1", "TIPARP", "CCNL1", "FABP5", "LDHB", "LSP1", "SPINT1", "HSPA6", "MATK", "SMAGP", "COL1A2", "COL5A1", "MMP2", "SLC30A2", "GKN1", "SERPINE2", "GALNT6", "MYOZ1",
"ARHGDIB", "GDPD3", "DSG2")
#Calculate score using dMAC based preterm/labor features:
seurat_obj_receiver <- PercentageFeatureSet(object = seurat_obj_receiver, features = mac_preterm_features,
  col.name = 'mac_preterm_features')

#dMAC1 DEGs:
mac1_list = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10,
                        logfc.threshold=0.25, latent.vars= c("nCount_RNA", "nFeature_RNA", "XIST", "mac_preterm_features",
                                            "MALAT1", "percent.mt"), test.use= "LR") %>% rownames_to_column("gene")
write.csv(mac1_list, file= "./PE_decidua_markers_logreg/dMAC1_PE_vs_lateC_preterm_corrected_040422.csv")

#Subset dMAC2:
receiver = "dMAC2"
seurat_obj_receiver= subset(data, idents = receiver)
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["time"]])
seurat_obj_receiver
table(Idents(seurat_obj_receiver))

#Calculate score using dMAC based preterm/labor features:
mac_preterm_features= c("CLDN5", "CD3E", "LILRB5", "FILIP1L", "RORA", "MMP19", "CD163", "SOD2", "BIRC3","SDC4", "ARL4C", "KDM6B", "ANGPTL4", "B4GALT1", "JUND", "NFKB1", "PDE4B", "B4GALT1", "JUND", "NFKB1", "PDE4B",
"MAFF", "NAMPT", "CD55", "C1QA", "CFLAR", "LMNA", "SAMSN1", "TIPARP", "CCNL1", "FABP5", "LDHB", "LSP1", "SPINT1", "HSPA6", "MATK", "SMAGP", "COL1A2", "COL5A1", "MMP2", "SLC30A2", "GKN1", "SERPINE2", "GALNT6", "MYOZ1",
"ARHGDIB", "GDPD3", "DSG2")
seurat_obj_receiver <- PercentageFeatureSet(object = seurat_obj_receiver, features = mac_preterm_features,
  col.name = 'mac_preterm_features')

#dMAC2 DEGs:
mac2_list = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10,
                        logfc.threshold=0.25, latent.vars= c("nCount_RNA", "nFeature_RNA", "XIST", "mac_preterm_features",
                                            "MALAT1", "percent.mt"), test.use= "LR") %>% rownames_to_column("gene")
write.csv(mac2_list, file= "./PE_decidua_markers_logreg/dMAC2_PE_vs_lateC_preterm_corrected_040422.csv")

#Subset dMono1 (or, dMono_LYZ):
receiver = "dMono_LYZ" #renamed as dMono1 in the manuscript.
seurat_obj_receiver= subset(data, idents = receiver)
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["time"]])
seurat_obj_receiver
table(Idents(seurat_obj_receiver))

#Calculate score using dMono/MAC based preterm & labor features:
dMono_preterm_features= c("CLDN5", "CD3E", "LILRB5", "FILIP1L", "RORA", "MMP19", "CD163", "SOD2", "BIRC3","SDC4", "ARL4C", "KDM6B", "ANGPTL4", "B4GALT1",
"JUND", "NFKB1", "PDE4B", "B4GALT1", "JUND", "NFKB1", "PDE4B", "MAFF", "NAMPT", "CD55", "C1QA", "CFLAR", "LMNA", "SAMSN1", "TIPARP", "CCNL1", "FABP5",
"LDHB", "LSP1", "SPINT1", "HSPA6", "MATK", "SMAGP", "COL1A2", "COL5A1", "MMP2", "SLC30A2", "GKN1", "SERPINE2", "GALNT6", "MYOZ1", "ARHGDIB", "GDPD3",
"DSG2", "TPSAB1", "SPP1")
seurat_obj_receiver <- PercentageFeatureSet(object = seurat_obj_receiver, features = dMono_preterm_features,
  col.name = 'dMono_preterm_features')

#dMono1 DEGs:
mono_list = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10,
                        logfc.threshold=0.25, latent.vars= c("nCount_RNA", "nFeature_RNA", "XIST", "dMono_preterm_features",
                        "MALAT1", "percent.mt"), test.use= "LR") %>% rownames_to_column("gene")
write.csv(mono_list, file= "./PE_decidua_markers_logreg/dMono_LYZ_PE_vs_lateC_preterm_corrected_040422.csv")
#write.csv(mono_list, file= "./PE_decidua_markers_logreg/dMono1_PE_vs_lateC_preterm_corrected_040422.csv")

#Subset dTcells:
receiver = "dTcell"
seurat_obj_receiver= subset(data, idents = receiver)
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["time"]])
seurat_obj_receiver
table(Idents(seurat_obj_receiver))

#T-cells preterm score:
seurat_obj_receiver <- PercentageFeatureSet(object = seurat_obj_receiver, features = c("SLC30A2", "GKN1", "SERPINE2", "GALNT6", "MYOZ1", "ARHGDIB", "GDPD3", "DSG2", "IFI27"),
  col.name = 'Tcells_preterm_features')

#DEGs of dTcells:
tc_list = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10,
                      logfc.threshold=0.25, latent.vars= c("nCount_RNA", "nFeature_RNA", "XIST", "Tcells_preterm_features",
                      "MALAT1", "percent.mt"), test.use= "LR") %>% rownames_to_column("gene")
write.csv(tc_list, file= "./PE_decidua_markers_logreg/dTcells_PE_vs_lateC_preterm_corrected_040422.csv")

#Subset DSC1:
receiver = "DSC1"
seurat_obj_receiver= subset(data, idents = receiver)
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["time"]])
seurat_obj_receiver
table(Idents(seurat_obj_receiver))

#Stromal preterm/labor score:
stromal_preterm_features= c("PLAC8", "ZBTB16", "IFI27", "SERPINE1", "SOD2", "PTGDS", "SLC30A2", "GKN1", "SERPINE2", "GALNT6", "MYOZ1", "ARHGDIB",
"GDPD3", "DSG2", "OLFML3", "MT2A") #ESAM is a drop-out in our dataset. 
seurat_obj_receiver <- PercentageFeatureSet(object = seurat_obj_receiver, features = stromal_preterm_features,
  col.name = 'stromal_preterm_features')

#DSC1 DEG list:
dsc1_list = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10,
            logfc.threshold=0.25, latent.vars= c("nCount_RNA", "nFeature_RNA",  "XIST", "stromal_preterm_features",  "MALAT1", "percent.mt"),
            test.use= "LR") %>% rownames_to_column("gene")
write.csv(dsc1_list, file= "./PE_decidua_markers_logreg/DSC1_PE_vs_lateC_preterm_corrected_040422.csv")

#Subset DSC_2:
receiver = "DSC2"
seurat_obj_receiver= subset(data, idents = receiver)
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["time"]])
seurat_obj_receiver
table(Idents(seurat_obj_receiver))

#Stromal preterm/labor score (same as DSC1):
stromal_preterm_features= c("PLAC8", "ZBTB16", "IFI27", "SERPINE1", "SOD2", "PTGDS", "SLC30A2", "GKN1", "SERPINE2", "GALNT6", "MYOZ1", "ARHGDIB",
"GDPD3", "DSG2", "OLFML3", "MT2A") #ESAM is a drop-out in our dataset. 
seurat_obj_receiver <- PercentageFeatureSet(object = seurat_obj_receiver, features = stromal_preterm_features, col.name = 'stromal_preterm_features')

#DSC2 DEG list:
dsc2_list = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10,
            logfc.threshold=0.25, latent.vars= c("nCount_RNA", "nFeature_RNA",  "XIST", "stromal_preterm_features",  "MALAT1", "percent.mt"),
            test.use= "LR") %>% rownames_to_column("gene")

write.csv(dsc2_list, file= "./PE_decidua_markers_logreg/DSC2_PE_vs_lateC_preterm_corrected_040422.csv")

#Subset dFB1:
receiver = "dFB1"
seurat_obj_receiver= subset(data, idents = receiver)
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["time"]])
seurat_obj_receiver
table(Idents(seurat_obj_receiver))

#Compute dFB preterm/labor score:
fb_preterm_features= c("PLAC8", "ZBTB16", "IFI27", "SERPINE1", "SOD2", "PTGDS", "SLC30A2", "GKN1", "SERPINE2", "GALNT6", "MYOZ1", "ARHGDIB", "GDPD3", "DSG2", "OLFML3", "MT2A", "DCN")
seurat_obj_receiver <- PercentageFeatureSet(object = seurat_obj_receiver, features = fb_preterm_features,
  col.name = 'fb_preterm_features')

#dFB1 DEGs:
fb1_list = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10,
           logfc.threshold=0.25, latent.vars= c("nCount_RNA", "nFeature_RNA",  "XIST", "fb_preterm_features", "MALAT1", "percent.mt"),
           test.use= "LR") %>% rownames_to_column("gene")

write.csv(fb1_list, file= "./PE_decidua_markers_logreg/dFB1_PE_vs_lateC_preterm_corrected_040422.csv")

#Subset dFB2:
receiver = "dFB2"
seurat_obj_receiver= subset(data, idents = receiver)
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["time"]])
seurat_obj_receiver
table(Idents(seurat_obj_receiver))

#Compute dFB preterm/labor score:
fb_preterm_features= c("PLAC8", "ZBTB16", "IFI27", "SERPINE1", "SOD2", "PTGDS", "SLC30A2", "GKN1", "SERPINE2", "GALNT6", "MYOZ1", "ARHGDIB", "GDPD3", "DSG2", "OLFML3", "MT2A", "DCN")
seurat_obj_receiver <- PercentageFeatureSet(object = seurat_obj_receiver, features = fb_preterm_features, col.name = 'fb_preterm_features')

#dFB2 DEGs:
fb2_list = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10,
           logfc.threshold=0.25, latent.vars= c("nCount_RNA", "nFeature_RNA",  "XIST", "fb_preterm_features", "MALAT1", "percent.mt"),
           test.use= "LR") %>% rownames_to_column("gene")

write.csv(fb2_list, file= "./PE_decidua_markers_logreg/dFB2_PE_vs_lateC_preterm_corrected_040422.csv")

#Subset dVEC:
receiver = "dVEC"
seurat_obj_receiver= subset(data, idents = receiver)
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["time"]])
seurat_obj_receiver
table(Idents(seurat_obj_receiver))
#LED preterm genes in Pique-Regi don't apply to dVEC- so, we prefer not to correct for them. 
seurat_obj_receiver <- PercentageFeatureSet(object = seurat_obj_receiver, features = c("PLAC8", "ZBTB16", "IFI27", "SERPINE1", "SOD2", "PTGDS", "OLFML3", "GKN1", "GLNT6"), 
                                            col.name = 'stromal_preterm_features')

#DEGs dVEC:
vec_list = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10,
           logfc.threshold=0.25, latent.vars= c("nCount_RNA", "nFeature_RNA", "XIST",  "MALAT1", "percent.mt", "stromal_preterm_features"),
           test.use= "LR") %>% rownames_to_column("gene")

#Save another file by running FindMarkers() without using the "stromal_preterm_features" as a covariate. The rest variables should stay.
write.csv(vec_list, file= "./PE_decidua_markers_logreg/dVEC_PE_vs_lateC_preterm_corrected_040422.csv")

#Case-II: Correcting for both stromal/LED don't cause major changes. File can be updated during revision. Alternatively, just correcting for core decidual genes is okay. 
seurat_obj_receiver <- PercentageFeatureSet(object = seurat_obj_receiver, features = c("PLAC8", "ZBTB16", "IFI27", "SERPINE1", "SOD2", "PTGDS"), 
                                            col.name = 'stromal_preterm_features') #DCN is a FB feature, so it's unused. 
seurat_obj_receiver <- PercentageFeatureSet(object = seurat_obj_receiver, features = c("FGL2", "EDN1", "OLFML3", "TXNRD2", "ANKRD1", "SLC30A2", "GKN1", "SERPINE2", "GALNT6", "MYOZ1", 
                "ARHGDIB", "GDPD3", "DSG2"), col.name = 'led_preterm_features')

vec_updated_list = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10,
           logfc.threshold=0.25, latent.vars= c("nCount_RNA", "nFeature_RNA", "XIST",  "MALAT1", "percent.mt", "stromal_preterm_features", "led_preterm_features"),
           test.use= "LR") %>% rownames_to_column("gene")

#Subset dLEC (or, dLECp):
receiver = "dLEC"
seurat_obj_receiver= subset(data, idents = receiver)
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["time"]])
seurat_obj_receiver
table(Idents(seurat_obj_receiver))

#Extract the LED preterm & labor signatures from Pique-Regi:
led_features= c("FGL2", "EDN1", "OLFML3", "TXNRD2", "ANKRD1", "SLC30A2", "GKN1", "SERPINE2", "GALNT6", "MYOZ1", "ARHGDIB", "GDPD3", "DSG2")
seurat_obj_receiver <- PercentageFeatureSet(object = seurat_obj_receiver, features = led_features, col.name = 'led_preterm_features')

#DEGs dLEC:
lec_list = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10,
           logfc.threshold=0.25, latent.vars= c("nCount_RNA", "XIST",  "MALAT1", "percent.mt", "led_preterm_features", "nGenes"),
                            test.use= "LR") %>% rownames_to_column("gene")

write.csv(lec_list, file= "./PE_decidua_markers_logreg/dLEC_PE_vs_lateC_preterm_corrected_040422.csv") #update Florian. 

#Subset dSMC:
receiver = "dSMC"
seurat_obj_receiver= subset(data, idents = receiver)
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["time"]])
seurat_obj_receiver
table(Idents(seurat_obj_receiver))

#FB preterm/labor features would apply to dSMC as well. 
#fb_preterm_features= c("PLAC8", "ZBTB16", "IFI27", "SERPINE1", "SOD2", "PTGDS", "SLC30A2", "GKN1", "SERPINE2", "GALNT6", "MYOZ1", "ARHGDIB", "GDPD3", "DSG2", "OLFML3", "MT2A", "DCN")
smc_features= c("SLC30A2", "GKN1", "SERPINE2", "GALNT6", "MYOZ1", "ARHGDIB", "GDPD3", "DSG2", "OLFML3", "MT2A",  "DCN", "IGF1")
seurat_obj_receiver <- PercentageFeatureSet(object = seurat_obj_receiver, features = smc_features, col.name = 'smc_preterm_features')

smc_list = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10,
           logfc.threshold=0.25, latent.vars= c("nCount_RNA", "nFeature_RNA",  "XIST", "smc_preterm_features", "MALAT1", "percent.mt"),
           test.use= "LR") %>% rownames_to_column("gene")
#Need to double check this file in publication folder or, update Florian about it. 
write.csv(smc_list, file= "./PE_decidua_markers_logreg/dSMC_PE_vs_lateC_preterm_corrected_040422.csv")

#Part-I: Placenta PE vs term control differentially expressed signatures:
#Read the placenta/villi RDS file. Converted from anndata H5AD via SeuratDisk.
#Note that, after Sept'22, the file path on Eils-HPC changed to: /dh-projects/preeclampsia_2022/analysis/placenta_atlas_2022/villi_scVI_harmonization/villi_seurat/SP014_SP082_SP136_placenta_updated_300322.rds
#Upon the Zenodo release, one can also simply use the villi RDS file for analysis.

data <- readRDS(file = "/data/analysis/preeclampsia_2019/placenta_atlas_2022/villi_scVI_harmonization/villi_seurat/SP014_SP082_SP136_placenta_updated_300322.rds")
data <- NormalizeData(object = data)
DefaultAssay(data) <- "RNA"

#Set idents to the lastest "cluster annotations"
Idents(object= data) <- 'leiden_subclusters_refined02'
table(Idents(data))

#Calculate nCount_RNA & nFeature_RNA:
nCount_RNA = colSums(x = data, slot = 'counts')
nFeature_RNA = colSums(x = GetAssayData(object = data, slot = "counts") > 0)
data <- AddMetaData(data, nFeature_RNA, col.name = "nFeature_RNA")
data <- AddMetaData(data, nCount_RNA, col.name= "nCount_RNA")

#Calculate percent.MT:
mt.genes <- rownames(data)[grep("^MT-",rownames(data))]
C <-GetAssayData(object = data, slot = "counts")
percent.mito <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
data <- AddMetaData(data, percent.mito, col.name = "percent.mt")

condition_oi = "late_preterm"
condition_reference = "late_term"

#Compile preterm genes for trophoblasts. Pique-Regi didn't focus on STB, so we restricted to what they showed for CTB.
#Additionally, SLIT2 & ROBO1 were included in the module-score based correction for as they are associated with risk for spontaneous preterm birth.
#Reference: doi:10.1371/journal.pgen.1008107 (2019).
tropho_preterm_features= c("HSPA1A", "CLIC3", "MYO7A", "SLIT2", "ROBO1","ATP1A4", "ZNF676", "CCND2", "SYDE1", "CDH5", "DACT3", "CST7", "CCL5",
          "KLRB1", "IL32", "CD44", "SLC27A2", "HSD11B1", "IL1RL1", "TNNI2", "RGCC")

#Subset vSCT_1 (named as vSTB1 in the manuscript):
receiver = "vSCT_1"
seurat_obj_receiver= subset(data, idents = receiver)
Idents(seurat_obj_receiver) <- "donor_id"
seurat_obj_receiver= subset(seurat_obj_receiver, idents= "Donor-557_2-villi", invert=TRUE) #remove the technical replicate 557_2.
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["time"]]) #Set the "condition" metadata as Ident.
seurat_obj_receiver
table(Idents(seurat_obj_receiver))

#Compute a preterm module score per cluster:
seurat_obj_receiver <- PercentageFeatureSet(object = seurat_obj_receiver, features = tropho_preterm_features,
  col.name = 'tropho_preterm_features')

#Since a villi can be either a male or a female, correcting for Y-genes is strongly recommended to prevent gender-specific bias in analysis.
stb1_list = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10,
            logfc.threshold=0.25, latent.vars= c("nCount_RNA", "nFeature_RNA", "XIST", "tropho_preterm_features", "MALAT1", "percent.mt", "pct_chrY"),
            test.use= "LR") %>% rownames_to_column("gene")

head(stb1_list)
write.csv(stb1_list, file= "./PE_decidua_markers_logreg/vSTB1_PE_vs_lateC_preterm_corrected_040422.csv")

#Subset vSCT_2 (named as vSTB2 in the manuscript):
receiver = "vSCT_2"
seurat_obj_receiver= subset(data, idents = receiver)
Idents(seurat_obj_receiver) <- "donor_id"
seurat_obj_receiver= subset(seurat_obj_receiver, idents= "Donor-557_2-villi", invert=TRUE) #remove the technical replicate 557_2.
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["time"]])
seurat_obj_receiver
table(Idents(seurat_obj_receiver))

#Compute a preterm module score per cluster:
seurat_obj_receiver <- PercentageFeatureSet(object = seurat_obj_receiver, features = tropho_preterm_features, col.name = 'tropho_preterm_features')

#Since a villi can be either a male or a female, correcting for Y-genes is strongly recommended to prevent gender-specific bias in analysis.
stb2_list = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10,
            logfc.threshold=0.25, latent.vars= c("nCount_RNA", "nFeature_RNA", "XIST", "tropho_preterm_features", "MALAT1", "percent.mt", "pct_chrY"),
            test.use= "LR") %>% rownames_to_column("gene")

head(stb2_list)
write.csv(stb2_list, file= "./PE_decidua_markers_logreg/vSTB2_PE_vs_lateC_preterm_corrected_040422.csv")

#Subset vSCT_juv (named as vSTBjuv in the manuscript):
receiver = "vSCTjuv"
seurat_obj_receiver= subset(data, idents = receiver)
Idents(seurat_obj_receiver) <- "donor_id"
seurat_obj_receiver= subset(seurat_obj_receiver, idents= "Donor-557_2-villi", invert=TRUE) #remove the technical replicate 557_2.
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["time"]])
seurat_obj_receiver
table(Idents(seurat_obj_receiver))

#Compute a preterm module score per cluster:
seurat_obj_receiver <- PercentageFeatureSet(object = seurat_obj_receiver, features = tropho_preterm_features, col.name = 'tropho_preterm_features')

#Since a villi can be either a male or a female, correcting for Y-genes is strongly recommended to prevent gender-specific bias in analysis.
juv_list = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10,
            logfc.threshold=0.25, latent.vars= c("nCount_RNA", "nFeature_RNA", "XIST", "tropho_preterm_features", "MALAT1", "percent.mt", "pct_chrY"),
            test.use= "LR", max.cells.per.ident=540) %>% rownames_to_column("gene")

head(juv_list)
write.csv(juv_list, file= "./PE_decidua_markers_logreg/vSTBjuv_PE_vs_lateC_preterm_corrected_040422.csv")


#Subset vVCT (named as vCTB in the manuscript):
receiver = "vVCT"
seurat_obj_receiver= subset(data, idents = receiver)
Idents(seurat_obj_receiver) <- "donor_id"
seurat_obj_receiver= subset(seurat_obj_receiver, idents= "Donor-557_2-villi", invert=TRUE) #remove the technical replicate 557_2.
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["time"]])
seurat_obj_receiver
table(Idents(seurat_obj_receiver))

#Compute a preterm module score per cluster:
seurat_obj_receiver <- PercentageFeatureSet(object = seurat_obj_receiver, features = tropho_preterm_features, col.name = 'tropho_preterm_features')

#Since a villi can be either a male or a female, correcting for Y-genes is strongly recommended to prevent gender-specific bias in analysis.
vct_list = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10,
            logfc.threshold=0.25, latent.vars= c("nCount_RNA", "nFeature_RNA", "XIST", "tropho_preterm_features", "MALAT1", "percent.mt", "pct_chrY"),
            test.use= "LR") %>% rownames_to_column("gene")

head(vct_list)
write.csv(vct_list, file= "./PE_decidua_markers_logreg/vCTB_PE_vs_lateC_preterm_corrected_040422.csv")

#vCTBp isn't suitable for DEG analysis given it's sparisty in both term controls & eoPE groups.

#In case of testing selected features in a group, such as in STB1, please do the following:
#gene_to_test: custom list of SASP (given as example)
gene_to_test= c('CD59', 'FLT1', 'HTRA4', 'PALM2-AKAP2', 'HBG2', 'SAT1', 'PAPPA2', 'AC093827.4', 'LEP', 'MIF', 'SH3PXD2A', 'SPECC1L', 'SLC2A11', 'LINC02484', 'AC026167.1', 'SASH1', 'SPATA13', 'LINC02109', 'DAB2', 'PHLDA2', 'SERPINB2', 'ATP6V0C', 'FNIP1', 'CGA', 'CARD18', 'EEF1G', 'SH3BP5', 'AL513164.1', 'POLR2J3', 'TJP2', 'MFSD4B', 'AC087857.1', 'FSTL3',
'LINC02163', 'RPL9', 'PRRG1', 'PEG3', 'MTRNR2L12', 'INSL4', 'ANKHD1', 'KDELR2', 'AFAP1', 'AL162718.1', 'HTRA1', 'FKBP1A', 'MICAL3', 'MYOF', 'PIK3CB', 'HBG1', 'LAMA5')

#Then, use, features= gene_to_test argument inside FindMarkers()
stb_list = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi,  ident.2 = condition_reference, min.pct = 0, logfc.threshold=0,
latent.vars= c("nCount_RNA", "nFeature_RNA",  "XIST", "tropho_preterm_features","MALAT1", "percent.mt", "pct_chrY"),
features= gene_to_test, test.use= "LR") %>% rownames_to_column("gene")
#Will return the avg_log2FC & adjusted p-values only for selected features (that can be later displayed as a heatmap, for instance)

#Subset vHBC:
receiver = "vHBC"
seurat_obj_receiver= subset(data, idents = receiver)
Idents(seurat_obj_receiver) <- "donor_id"
seurat_obj_receiver= subset(seurat_obj_receiver, idents= "Donor-557_2-villi", invert=TRUE) #remove the technical replicate 557_2.
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["time"]])
seurat_obj_receiver
table(Idents(seurat_obj_receiver))

#Even though vHBC is analogous to macrophages, Pique-Regi's MAC preterm/labour signatures were particularly explained for maternal decidua.
#Since the features weren't relevant for villi HBC, we didn't correct for it. Nonetheless, an user is free to try it out.
#Since a villi can be either a male or a female, correcting for Y-genes is strongly recommended to prevent gender-specific bias in analysis.
hbc_list = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10,
            logfc.threshold=0.25, latent.vars= c("nCount_RNA", "nFeature_RNA", "XIST", "MALAT1", "percent.mt", "pct_chrY"),
            test.use= "LR", max.cells.per.ident=548) %>% rownames_to_column("gene") #downsampled as reduction in #cells in PE is significant (FDR < 0.05)

#Results are not affected in any major way without downsampling. Save the results without too (since, higher #cells might indicate disease relevance).
hbc_list = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10,
            logfc.threshold=0.25, latent.vars= c("nCount_RNA", "nFeature_RNA", "XIST", "MALAT1", "percent.mt", "pct_chrY"),
            test.use= "LR") %>% rownames_to_column("gene")
head(hbc_list)
write.csv(hbc_list, file= "./PE_decidua_markers_logreg/vHBC_PE_vs_lateC_preterm_corrected_040422.csv")

#Subset vPAMM:
receiver = "vPAMM" 
seurat_obj_receiver= subset(data, idents = receiver)
Idents(seurat_obj_receiver) <- "donor_id"
seurat_obj_receiver= subset(seurat_obj_receiver, idents= "Donor-557_2-villi", invert=TRUE) #remove the technical replicate 557_2.
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["time"]])
seurat_obj_receiver
table(Idents(seurat_obj_receiver))

#Since a villi can be either a male or a female, correcting for Y-genes is strongly recommended to prevent gender-specific bias in analysis.
#Correction for dMAC_preterm genes (decidual) from Pique-Regi didn't have a major influence here. 
pamm_list = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10,
            logfc.threshold=0.25, latent.vars= c("nCount_RNA", "nFeature_RNA", "XIST", "MALAT1", "percent.mt", "pct_chrY"),
            test.use= "LR") %>% rownames_to_column("gene")

head(pamm_list)
write.csv(pamm_list, file= "./PE_decidua_markers_logreg/vPAMM_PE_vs_lateC_preterm_corrected_040422.csv")

#No preterm/labor signatures were described in Pique-Regi that could be used villi cell-types below:
#Subset vVEC:
receiver = "vVEC"
seurat_obj_receiver= subset(data, idents = receiver)
Idents(seurat_obj_receiver) <- "donor_id"
seurat_obj_receiver= subset(seurat_obj_receiver, idents= "Donor-557_2-villi", invert=TRUE) #remove the technical replicate 557_2.
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["time"]])
table(Idents(seurat_obj_receiver))

#Since a villi can be either a male or a female, correcting for Y-genes is strongly recommended to prevent gender-specific bias in analysis.
#Correcting for trophoblast based preterm genes is irrelevant here.
vec_list = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10,
            logfc.threshold=0.25, latent.vars= c("nCount_RNA", "nFeature_RNA", "XIST", "MALAT1", "percent.mt", "pct_chrY"),
            test.use= "LR") %>% rownames_to_column("gene")

head(vec_list)
write.csv(vec_list, file= "./PE_decidua_markers_logreg/vVEC_PE_vs_lateC_preterm_corrected_040422.csv")

#Subset vFB:
receiver = "vFB1" 
seurat_obj_receiver= subset(data, idents = receiver)
Idents(seurat_obj_receiver) <- "donor_id"
seurat_obj_receiver= subset(seurat_obj_receiver, idents= "Donor-557_2-villi", invert=TRUE)
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["time"]])
table(Idents(seurat_obj_receiver))

#Corrected as few genes make sense for villi:
fb_preterm_features= c("PLAC8", "ZBTB16", "IFI27", "SERPINE1", "SOD2", "PTGDS", "SLC30A2", "GKN1", "SERPINE2", "GALNT6", "MYOZ1", "ARHGDIB", "GDPD3", "DSG2", "OLFML3", "MT2A", "DCN")
seurat_obj_receiver <- PercentageFeatureSet(object = seurat_obj_receiver, features = fb_preterm_features, col.name = 'fb_preterm_features')

#Since a villi can be either a male or a female, correcting for Y-genes is strongly recommended to prevent gender-specific bias in analysis.
fb_list = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10,
            logfc.threshold=0.25, latent.vars= c("nCount_RNA", "nFeature_RNA", "fb_preterm_features", "MALAT1", "percent.mt", "pct_chrY"),
            test.use= "LR") %>% rownames_to_column("gene")
head(fb_list)
write.csv(fb_list, file= "./PE_decidua_markers_logreg/vFB_PE_vs_lateC_preterm_corrected_040422.csv")

#Subset vMC:
receiver = "vMC"
seurat_obj_receiver= subset(data, idents = receiver)
Idents(seurat_obj_receiver) <- "donor_id"
seurat_obj_receiver= subset(seurat_obj_receiver, idents= "Donor-557_2-villi", invert=TRUE) #remove the technical replicate 557_2.
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["time"]])
table(Idents(seurat_obj_receiver))

#Since a villi can be either a male or a female, correcting for Y-genes is strongly recommended to prevent gender-specific bias in analysis.
#Correcting for trophoblast based preterm genes is irrelevant here.
myocyte_list = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10,
            logfc.threshold=0.25, latent.vars= c("nCount_RNA", "nFeature_RNA", "XIST", "MALAT1", "percent.mt", "pct_chrY"),
            test.use= "LR", max.cells.per.ident=340) %>% rownames_to_column("gene")
head(myocyte_list)
write.csv(myocyte_list, file= "./PE_decidua_markers_logreg/vMC_PE_vs_lateC_preterm_corrected_040422.csv")

#Subset vTcell:
receiver = "vTcell"
seurat_obj_receiver= subset(data, idents = receiver)
Idents(seurat_obj_receiver) <- "donor_id"
seurat_obj_receiver= subset(seurat_obj_receiver, idents= "Donor-557_2-villi", invert=TRUE) #remove the technical replicate 557_2.
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["time"]])
table(Idents(seurat_obj_receiver))

#seurat_obj_receiver <- PercentageFeatureSet(object = seurat_obj_receiver, features = c("SLC30A2", "GKN1", "SERPINE2", "GALNT6", "MYOZ1", "ARHGDIB", "GDPD3", "DSG2", "IFI27"),
  #col.name = 'Tcells_preterm_features') #T-cells are maternal in origin. 
#Since a villi can be either a male or a female, correcting for Y-genes is strongly recommended to prevent gender-specific bias in analysis.
#Correcting for trophoblast based preterm genes is irrelevant here.
tcell_list = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10,
            logfc.threshold=0.25, latent.vars= c("nCount_RNA", "nFeature_RNA", "XIST", "MALAT1", "percent.mt", "pct_chrY"),
            test.use= "LR") %>% rownames_to_column("gene")
head(tcell_list)
write.csv(tcell_list, file= "./PE_decidua_markers_logreg/vTcell_PE_vs_lateC_preterm_corrected_040422.csv")
