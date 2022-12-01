Sys.setenv(LANG="en_US.UTF-8")
Sys.setenv(LC_ALL="en_US.UTF-8")

suppressPackageStartupMessages(library(tidyverse)) 
suppressPackageStartupMessages(library(tibble)) 
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(ggplot2)) 
suppressPackageStartupMessages(library(cowplot)) 
suppressPackageStartupMessages(library(ComplexHeatmap)) 
suppressPackageStartupMessages(library(Connectome)) 

#Main figure 4A: dDSTB ligand interactions with decidual endothelial receptors.
#Only dVEC & dSMC are shown in the manuscript given the biological feasibility of interactions. 
#Read the decidua RDS file. Converted from anndata H5AD via SeuratDisk.
#Note that, after Sept'22, the file path on Eils-HPC changed to: /dh-projects/preeclampsia_2022/analysis/placenta_atlas_2022/decidua_scVI_harmonization/decidua_seurat/Raw_decidua_OD_140322.rds
#Upon the Zenodo release, one can also simply use the decidua RDS file for analysis.
data <- readRDS(file = "/data/analysis/preeclampsia_2019/placenta_atlas_2022/decidua_scVI_harmonization/decidua_seurat/Raw_decidua_OD_140322.rds")
data <- NormalizeData(object = data)
DefaultAssay(data) <- "RNA" 
Idents(object= data) <- 'time' #Set idents to "time" (condition separating PE vs late controls)
data <- subset(data, idents= "late_preterm") #Only subset PE dataset. 
Idents(object= data) <- 'celltype_v5' #Set idents to current celltype annotations.  
data <- subset(data, idents= c("dVEC", "dSMC", "deported_SCT", "dLEC", "dLECp")) #subset clusters of interest. deported_SCT renamed as "dDSTB" in the manuscript. 
table(Idents(data)) #check the cell composition: dLEC=1830, dLECp=235, dSMC= 421, dVEC=1151, deported_SCT=980. 

#Here, after data normalization, we identify the ligand and receptor genes that have mapped to the dDSTB-endo object, 
#scale those features, and generate the connectome. 
connectome.genes <- union(Connectome::ncomms8866_human$Ligand.ApprovedSymbol,Connectome::ncomms8866_human$Receptor.ApprovedSymbol)
genes <- connectome.genes[connectome.genes %in% rownames(data)]
data <- ScaleData(data, features = genes)

#Create the R-L connectome, FANTOM5 database is used as default: 
#Compute DOR.  High DOR is an indicative of high specificity and sensitivity with a low rate of false positives and false negatives.
data.con <- CreateConnectome(data,species = 'human',min.cells.per.ident=75,p.values = T,calculate.DOR = TRUE) 
#To use custom DB like Cellchat, just read.csv() the file & do the following: 
#data.con <- CreateConnectome(data,species = 'human', LR.database= 'custom', custom.list= cellchat_db, calculate.DOR= T)

#Filter connectome to save results for ligands expressed by deported_SCT (or, dDSTB) & endothelial/SMC receptors:
data_filter <- FilterConnectome(data.con, sources.include= "deported_SCT", targets.include= c("dVEC", "dSMC")) 
#Note, the above data table unfiltered by further statistics like pct & DOR. 
#After filtering by source & targets of interest,we've 4844 interactions in total. Picks up pathologically relevant EBI3-IL6ST both in dVEC/dSMC
write.csv(data_filter, file= "dSTB_dVEC_dSMC_connectomics_FANTOM5_25052022.csv")

#For figure 4A, firstly, we took "data_filter.csv" table and then filtered the interaction list using pct.source (senders) >= 25% 
#and pct.target (receivers) >= 20% (empahsizes on cluster-specific communication). Next, we filtered by DOR.source > 3 and ligand expression > 1.5. 
#If a user is interested, they can apply z-score cut-offs too. Such as filter the data to only include edges with a ligand and receptor 
#z-score above 0.25 (focus on cell-type specific communication patterns) & where both candidates are expressed in at least 10% of the cells 
#in their respective cluster (hence ensures statistical confidence):
data.con2 <- FilterConnectome(data_filter,min.pct = 0.1,min.z = 0.25,remove.na = T) #20 edges would be retained post-filtering. 


#Main figure 4C: vSTB secreted ligands interaction with maternal vessels & dSMC. 
#Load the Seurat object (RDS); exclude "Donor-557_2-villi" technical replicate from analysis. 
data <- readRDS(file= "/data/analysis/preeclampsia_2019/placenta_atlas_2022/maternal_fetal_analysis/Decidua_villi_cellchat_310322.rds") #already subset for relevant clusters. 
Idents(data) <- "donor_id" #Set idents to current donor. 
data <- subset(data, idents= "Donor-557_2-villi", invert=TRUE) #Remove the TR 557_2:
table(Idents(data)) #confirm 557_2 is excluded. 
Idents(data) <- "celltype_final" #Set idents to cell type labels for performing RL interaction analysis:
#data <- subset(data, idents= c("vSCT_1", "vSCT_2", "vSCT_juv", "dVEC", "dSMC", "dLEC")) #RDS is already subsetted for desired clusters. 
table(Idents(data)) #look at the cell composition table. 

#Firstly, we identify ligands and receptors expressed in our dataset:
connectome.genes <- union(Connectome::ncomms8866_human$Ligand.ApprovedSymbol,Connectome::ncomms8866_human$Receptor.ApprovedSymbol)
genes <- connectome.genes[connectome.genes %in% rownames(data)]

#Split the object by condition (i.e., PE or late controls:
#In the metadata, "time" column describes if a cell is lateC or PE. 
data.list <- SplitObject(data,split.by = 'time')

#Normalize, scale, and create Connectome for both lateC & PE objects:
data.con.list <- list()
for (i in 1:length(data.list)){
  data.list[[i]] <- NormalizeData(data.list[[i]])
  data.list[[i]] <- ScaleData(data.list[[i]],features = rownames(data.list[[i]]))
  data.con.list[[i]] <- CreateConnectome(data.list[[i]],species = 'human')
  #data.con.list[[i]] <- CreateConnectome(data.list[[i]],species = 'human', LR.database= 'custom', custom.list= cellchat_df, calculate.DOR= T) #Run this line for custom DB such as Cellchat. 
}
names(data.con.list) <- names(data.list)

#Build the differential connectome: 
diff <- DifferentialConnectome(data.con.list[[1]], data.con.list[[2]]) 

#Find edges of statistical significance: 
celltypes <- as.character(unique(Idents(data))) 
celltypes.stim <- paste(celltypes, 'late_preterm', sep = '_')
celltypes.ctrl <- paste(celltypes, 'late_term', sep = '_')
data$celltype.condition <- paste(Idents(data), data$time, sep = "_")
data$celltype <- Idents(data)
Idents(data) <- "celltype.condition" #so, each cluster is further segregated by PE & late controls. 

#Identify which ligands and receptors, for which cell populations, have an adjusted p-value < 0.05 based on a LR test
#Do same preterm/labor gene correction as done for DEG. 
tropho_features= c("HSPA1A", "CLIC3", "MYO7A", "SLIT2", "ROBO1","ATP1A4", "ZNF676", "CCND2", "SYDE1", "CDH5", "DACT3", "CST7", "CCL5",
"KLRB1", "IL32", "CD44", "SLC27A2", "HSD11B1", "IL1RL1", "TNNI2", "RGCC") 
led_features= c("FGL2", "EDN1", "OLFML3", "TXNRD2", "ANKRD1", "SLC30A2", "GKN1", "SERPINE2", "GALNT6", "MYOZ1", "ARHGDIB", "GDPD3", "DSG2")

#Compute module scores using the preterm scores & add them to metadata:
data <- PercentageFeatureSet(object = data, features = tropho_features, col.name = 'tropho_preterm_features')
data <- PercentageFeatureSet(object = data, features = led_features, col.name = 'led_preterm_features')

diff.p <- data.frame() #Store the DEG results to a data frame. 
for (i in 1:length(celltypes)){
  temp <- FindMarkers(data,
                      ident.1 = celltypes.stim[i],
                      ident.2 = celltypes.ctrl[i],
                      verbose = FALSE,
                      features = genes,
                      min.pct = 0.1, logfc.threshold=0.1, 
                      latent.vars= c("nCount_RNA", "nFeature_RNA",  "XIST", "tropho_preterm_features", "led_preterm_features", 
                      "MALAT1", "percent.mt", "pct_chrY"), test.use= "LR")
  temp2 <- subset(temp, p_val_adj < 0.05)
  if (nrow(temp2)>0){
  temp3 <- data.frame(genes = rownames(temp2),cells = celltypes[i])
  diff.p <- rbind(diff.p, temp3)
  }
}
diff.p$cell.gene <- paste(diff.p$cells,diff.p$genes,sep = '.') #compute a cell.gene data frame which lists which RL maps to what cluster. 

#Alternatively, one can just input the LR DEG results like below. For a reference file, look at: DEG_input_example.csv (first column: "genes", second: "cells")
#The input list is already filtered by DEG criteria as explained for PE vs term controls. 
deg_file= read.csv("/dh-projects/preeclampsia_2022/analysis/placenta_atlas_2022/maternal_fetal_analysis/vSCT_endo_DEG_input.csv")
diff.p$cell.gene <- paste(deg_file$cells,deg_file$genes,sep = '.') 

#Filter differential connectome to only include significantly perturbed edges
diff$source.ligand <- paste(diff$source,diff$ligand,sep = '.')
diff$target.receptor <- paste(diff$target,diff$receptor,sep = '.')

#Subset interactions where both receptor & ligand are differentially expressed in our dataset (a/c to Logistic Regression)
diff.sub <- subset(diff,source.ligand %in% diff.p$cell.gene & target.receptor %in% diff.p$cell.gene)

#When ligand is upregulated (since we emphasized only on secreted ligand pressure for eoPE):
diff.up <- subset(diff.sub,ligand.norm.lfc > 0)
write.csv(diff.up, file= "vSTB_ligands_upregulated.csv")  
#Aggregated results for both FANTOM5 & CellchatDB in Excel & dropped duplicate interaction-pairs. 
#We emphasized only "secreted ligands" in the manuscript as they are cross maternal-fetal interface & can bind to vessels and SMC receptors in the mother. 

#Extended Data figure 8B: Intra villi communication among mesoendothelial & immune compartments. 
#One can read the villi RDS file upon Zenodo release. 
data <- readRDS(file = "/data/analysis/preeclampsia_2019/placenta_atlas_2022/villi_scVI_harmonization/villi_seurat/SP014_SP082_SP136_placenta_updated_300322.rds")
Idents(data) <- "donor_id" #Set idents to current donor. 
data <- subset(data, idents= "Donor-557_2-villi", invert=TRUE) #Remove the technical replicate 557_2:
Idents(data) <- "celltype_final" #Set idents to cell type labels for performing RL interaction analysis:
data <- subset(data, idents= c("vVEC", "vVCT", "vTcell", "vMC", "vHBC", "PAMM")) 
table(Idents(data)) #look at the cell composition table. 

#Firstly, we identify ligands and receptors expressed in our dataset:
connectome.genes <- union(Connectome::ncomms8866_human$Ligand.ApprovedSymbol,Connectome::ncomms8866_human$Receptor.ApprovedSymbol)
genes <- connectome.genes[connectome.genes %in% rownames(data)]

#Split the object by condition (i.e., PE or late controls:
#In the metadata, "time" column describes if a cell is lateC or PE. 
data.list <- SplitObject(data,split.by = 'time')

#Normalize, scale, and create Connectome for both lateC & PE objects:
data.con.list <- list()
for (i in 1:length(data.list)){
  data.list[[i]] <- NormalizeData(data.list[[i]])
  data.list[[i]] <- ScaleData(data.list[[i]],features = rownames(data.list[[i]]))
  data.con.list[[i]] <- CreateConnectome(data.list[[i]],species = 'human',calculate.DOR= T)
  #data.con.list[[i]] <- CreateConnectome(data.list[[i]],species = 'human', LR.database= 'custom', custom.list= cellchat_df, calculate.DOR= T) #Run this line for custom DB such as Cellchat. 
}
#Steps below applies independent of DB used. 
names(data.con.list) <- names(data.list)

diff <- DifferentialConnectome(data.con.list[[1]], data.con.list[[2]]) #Build the differential connectome.
#Filter differential connectome to only include significantly perturbed edges
diff$source.ligand <- paste(diff$source,diff$ligand,sep = '.')
diff$target.receptor <- paste(diff$target,diff$receptor,sep = '.') 

#Find edges of statistical significance:
#celltypes <- as.character(unique(Idents(data)))
#celltypes.stim <- paste(celltypes, 'late_preterm', sep = '_')
#celltypes.ctrl <- paste(celltypes, 'late_term', sep = '_')
#data$celltype.condition <- paste(Idents(data), data$time, sep = "_")
#data$celltype <- Idents(data)
#Idents(data) <- "celltype.condition" #so, each cluster is further segregated by PE & late controls. 
#Insert a loop that iterate over each cluster, calculates DEGs & store them into a df (see above for example)

#Instead of recalculating the differentially expressed receptors & ligands, simply input the already filtered DEGs for selected villi cell-types.
#DEGs were calculated using Logistic Regression. For reference, check: PE_vs_lateC_marker_analysis_LogisticRegression.r 
deg_file= read.csv("/dh-projects/preeclampsia_2022/analysis/placenta_atlas_2022/villi_scVI_harmonization/villi_immune_RL/villi_immune_DEGs.csv")
diff.p$cell.gene <- paste(deg_file$cells,deg_file$genes,sep = '.') 

#Subset RL interactions for pairs differentially expressed in our dataset: 
diff.sub <- subset(diff,source.ligand %in% diff.p$cell.gene & target.receptor %in% diff.p$cell.gene)
write.csv(diff.sub, file= "Villi_RL_FANTOM5_connectomics.csv") 
#write.csv(diff.sub, file= "Villi_RL_CCDB_connectomics.csv") #for saving results for CellchatDB. 

#Read the CSV file consisting of aggregated Villi FANTOM5 & CCDB results (dropped duplicated interaction-pairs in Excel):
villi_df= read.csv("/dh-projects/preeclampsia_2022/analysis/placenta_atlas_2022/Decidua_villi_interaction_extended/villi_connectomics_plot_090522.csv") 
diff.up <- subset(villi_df,ligand.norm.lfc > 0, receptor.norm.lfc > 0) #When both receptor & ligand are upregulated: edge activation 
write.csv(diff.up, file= "Villi_RL_upregulated.csv") 
diff.down <- subset(villi_df,ligand.norm.lfc < 0, receptor.norm.lfc < 0) #When both are downregulated: edge deactivation. 
write.csv(diff.down, file= "Villi_RL_downregulated.csv") 

#Compile a union of unique source & targets for both upregulated & downregulated villi df(s): 
if (is.null(cols.use)){  nodes <- as.character(unique(union(diff.up$source, diff.up$target)))  
cols.use <- hue_pal()(length(nodes))  names(cols.use) <- nodes  cols.use <- data.frame(cols.use)  
cols.use$cell <- rownames(cols.use)}else{  cols.use <- data.frame(cols.use)  cols.use$cell <- rownames(cols.use)}

if (is.null(cols.use)){  nodes <- as.character(unique(union(diff.down$source, diff.down$target)))  
cols.use <- hue_pal()(length(nodes))  names(cols.use) <- nodes  cols.use <- data.frame(cols.use)  
cols.use$cell <- rownames(cols.use)}else{  cols.use <- data.frame(cols.use)  cols.use$cell <- rownames(cols.use)}

grid.col <- as.vector(c("#fd96A9", "#00b3b3", "#a799b7", "#63264a", "#fe6776", "#c0c999")) #set colors as a vector. 
names(grid.col) <- c("vVEC", "vVCT", "vTcell", "vMC", "vHBC", "PAMM") #provide cluster mappings. 

#Extended Data Figure 8B: 
pdf("Villi_RL_upregulated_090522.pdf", paper= "special", h=10, w=10)
CircosDiff(diff.up, cols.use = grid.col) 
pdf("Villi_RL_downregulated_090522.pdf", paper= "special", h=10, w=10)
CircosDiff(diff.down, cols.use = grid.col) 


#Extended Data figure 8A: Intra decidua communication among meso-endothelial & immune compartments. 
data <- readRDS(file = "/data/analysis/preeclampsia_2019/placenta_atlas_2022/decidua_scVI_harmonization/decidua_seurat/Raw_decidua_OD_140322.rds")
data <- NormalizeData(object = data)
DefaultAssay(data) <- "RNA"
Idents(object= data) <- 'celltype_v5' #Set idents to current celltype annotations.  
data <- subset(data, idents= c("dMAC1", "dTcell", "dLEC", "dFB2", "DSC1", "dLECp", "dMAC2", "DSC2", "dSMC", "dNK2", "dMono_LYZ", "dVEC", "dFB1", "dNK1")) #subset clusters of interest. dMono_LYZ renamed as "dMono1" in the manuscript. 
table(Idents(data)) #cell composition. 

#Firstly, we identify ligands and receptors expressed in our dataset:
connectome.genes <- union(Connectome::ncomms8866_human$Ligand.ApprovedSymbol,Connectome::ncomms8866_human$Receptor.ApprovedSymbol)
genes <- connectome.genes[connectome.genes %in% rownames(data)]  

#Split the object by condition (i.e., PE or late controls:
#In the metadata, "time" column describes if a cell is lateC or PE. 
data.list <- SplitObject(data,split.by = 'time')

#Normalize, scale, and create Connectome for both lateC & PE objects: 
data.con.list <- list() 
for (i in 1:length(data.list)){
  data.list[[i]] <- NormalizeData(data.list[[i]])  
  data.list[[i]] <- ScaleData(data.list[[i]],features = rownames(data.list[[i]]))
  data.con.list[[i]] <- CreateConnectome(data.list[[i]],species = 'human',calculate.DOR= T)
  #data.con.list[[i]] <- CreateConnectome(data.list[[i]],species = 'human', LR.database= 'custom', custom.list= cellchat_df, calculate.DOR= T) #Run this line for custom DB such as Cellchat. 
}
#Steps below applies independent of DB used. 
names(data.con.list) <- names(data.list)
diff <- DifferentialConnectome(data.con.list[[1]], data.con.list[[2]]) #Build the differential connectome. 
#Filter differential connectome to only include significantly perturbed edges
diff$source.ligand <- paste(diff$source,diff$ligand,sep = '.')
diff$target.receptor <- paste(diff$target,diff$receptor,sep = '.')

#Find edges of statistical significance:
#celltypes <- as.character(unique(Idents(data)))
#celltypes.stim <- paste(celltypes, 'late_preterm', sep = '_')
#celltypes.ctrl <- paste(celltypes, 'late_term', sep = '_')
#data$celltype.condition <- paste(Idents(data), data$time, sep = "_")
#data$celltype <- Idents(data)
#Idents(data) <- "celltype.condition" #so, each cluster is further segregated by PE & late controls. 
#Insert a loop that iterate over each cluster, calculates DEGs & store them into a df (see above for example)

#Instead of recalculating the differentially expressed receptors & ligands, simply input the already filtered DEGs for selected decidual cell-types. 
deg_file= read.csv("/dh-projects/preeclampsia_2022/analysis/placenta_atlas_2022/decidua_scVI_harmonization/decidua_seurat/decidua_immune_input.csv")
diff.p$cell.gene <- paste(deg_file$cells,deg_file$genes,sep = '.') 

#Subset RL interactions for pairs differentially expressed in our dataset: 
diff.sub <- subset(diff,source.ligand %in% diff.p$cell.gene & target.receptor %in% diff.p$cell.gene)
write.csv(diff.sub, file= "Decidua_RL_FANTOM5_connectomics.csv") 
#write.csv(diff.sub, file= "Decidua_RL_CCDB_connectomics.csv") #for saving results for CellchatDB. 

#Read the CSV file consisting of aggregated FANTOM5 & CCDB results (dropped duplicates in Excel):
decidua_df= read.csv("/dh-projects/preeclampsia_2022/analysis/placenta_atlas_2022/Decidua_villi_interaction_extended/Decidua_connectomics_plot.csv") 
diff.up <- subset(decidua_df,ligand.norm.lfc > 0, receptor.norm.lfc > 0) #When both receptor & ligand are upregulated: edge activation 
write.csv(diff.up, file= "Decidua_RL_upregulated.csv") 
diff.down <- subset(decidua_df,ligand.norm.lfc < 0, receptor.norm.lfc < 0) #When both are downregulated: edge deactivation. 
write.csv(diff.down, file= "Decidua_RL_downregulated.csv") 

#Compile a union of unique source & targets for both upregulated & downregulated decidua df(s): 
if (is.null(cols.use)){  nodes <- as.character(unique(union(diff.up$source, diff.up$target)))  
cols.use <- hue_pal()(length(nodes))  names(cols.use) <- nodes  cols.use <- data.frame(cols.use)  
cols.use$cell <- rownames(cols.use)}else{  cols.use <- data.frame(cols.use)  cols.use$cell <- rownames(cols.use)}

if (is.null(cols.use)){  nodes <- as.character(unique(union(diff.down$source, diff.down$target)))  
cols.use <- hue_pal()(length(nodes))  names(cols.use) <- nodes  cols.use <- data.frame(cols.use)  
cols.use$cell <- rownames(cols.use)}else{  cols.use <- data.frame(cols.use)  cols.use$cell <- rownames(cols.use)}

grid.col <- as.vector(c('#b3b300', '#c2d6d6', '#ff66b3', '#994d00', '#33ccff', '#99004d', '#4d4d00', '#cc3300', '#000080', '#80d4ff', '#004d00', '#9966ff', '#ffcc99', '#006699')) 
names(grid.col) <- c("dMAC1", "dTcell", "dLEC", "dFB2", "DSC1", "dLECp", "dMAC2", "DSC2", "dSMC", "dNK2", "dMono_LYZ", "dVEC", "dFB1", "dNK1")

#Save figures for Extended Data Figure 8A. 
pdf("Decidua_RL_upregulated_090522", paper= "special", h=18, w=18)
CircosDiff(diff.up, cols.use = grid.col)  
pdf("Decidua_RL_downregulated_090522", paper= "special", h=15, w=15)
CircosDiff(diff.down, cols.use = grid.col)  

#Extended Data figure 8C: RL interaction analysis showing invasive dEVT (extravillous trophoblast) ligands communicating with decidual VEC & SMC in PE. 
data <- readRDS(file = "/data/analysis/preeclampsia_2019/placenta_atlas_2022/decidua_scVI_harmonization/decidua_seurat/Raw_decidua_OD_140322.rds")
data <- NormalizeData(object = data)
DefaultAssay(data) <- "RNA" 
Idents(object= data) <- 'time' #Set idents to "time" (condition separating PE vs late controls)
data <- subset(data, idents= "late_preterm") #Only subset PE dataset.  
Idents(object= data) <- 'celltype_v5' #Set idents to current celltype annotations.  
data <- subset(data, idents= c("dEVT", "dVEC", "dSMC", "dLEC", "dLECp", "dEpC", "DSC1", "DSC2", "dFB1", "dFB2", "dMSC", 
"dNK1", "dNK2", "dNK_prol", "dMAC1", "dMAC2", "dTcell", "dMono_LYZ", "dMono2", "dGranulocyte")) #dMono_LYZ= dMono1; dGranulocyte= dGranul, dNK_prol= dNKp (in manuscript)
table(Idents(data)) #check the cell composition 

#Here, after data normalization, we identify the ligand and receptor genes that have mapped to the dDSTB-endo object, 
#scale those features, and generate the connectome. 
connectome.genes <- union(Connectome::ncomms8866_human$Ligand.ApprovedSymbol,Connectome::ncomms8866_human$Receptor.ApprovedSymbol)
genes <- connectome.genes[connectome.genes %in% rownames(data)]
data <- ScaleData(data, features = genes)

#Create the R-L connectome, FANTOM5 database is used as default: 
#Compute DOR.  High DOR is an indicative of high specificity and sensitivity with a low rate of false positives and false negatives.
data.con <- CreateConnectome(data,species = 'human',min.cells.per.ident=75,p.values = T,calculate.DOR = TRUE) 
#To use custom DB like Cellchat, just read.csv() the file & do the following: 
#data.con <- CreateConnectome(data,species = 'human', LR.database= 'custom', custom.list= cellchat_db, calculate.DOR= T)

#Filter connectome to save results for ligands from dEVT to any targets (dEVT receptors were filtered out later in Excel).
#Minimum z-score=0.25 to emphasize on cell-type specific interaction patterns. 
data_filter <- FilterConnectome(data.con, sources.include= "dEVT", min.z = 0.25,remove.na = T) 
write.csv(data_filter, file= "PE_dEVT_ligands_connectomics_ZS0.25") 

#For the interactions with dEVT, relatively robust criteria were further used for narrowing down the important interaction partners (from an initial list of > 10K pairs). 
#Particularly, DOR.source=5, edge strength (weight_norm; product of the receptor and ligand expression)=3 
#and minimum percentage of ligand expressing source of 50% (percent.source) was used to ensure cell specific communication. 
#Aggregated results of FANTOM5 & CCDB in Excel & excluded duplicated RL interaction pairs. 
