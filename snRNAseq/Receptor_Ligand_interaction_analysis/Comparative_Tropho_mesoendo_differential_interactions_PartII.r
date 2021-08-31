Sys.setenv(LANG="en_US.UTF-8")
Sys.setenv(LC_ALL="en_US.UTF-8")

suppressPackageStartupMessages(library(CellChat))

#Read the cellchat object consisting of RL data (Trimean method). Refer to notebook:  
data <- readRDS(file= "/data/analysis/preeclampsia_2019/analysis/images/cellphonedb_analysis/TrophoL_mesenchymal_analysis/TriMean_cellchat_comparison_LateC_vs_PE.rds")

data

data@meta$datasets = factor(data@meta$datasets, levels = c("LC", "PE")) # set factor level

plotGeneExpression(data, signaling = "CALCR", split.by = "datasets", color.use = c("grey", "darkblue"))

pdf("./TrophoL_mesenchymal_analysis/comparison_violin_plots/PE_tropho_CALCR.pdf", paper= "special")

plotGeneExpression(data, signaling = "CALCR", split.by = "datasets", color.use = c("grey", "darkblue"))

dev.off()

pdf("./TrophoL_mesenchymal_analysis/comparison_violin_plots/PE_tropho_GDF.pdf", paper= "special")

plotGeneExpression(data, signaling = "GDF", split.by = "datasets", color.use = c("grey", "darkblue"))

dev.off()

pdf("./TrophoL_mesenchymal_analysis/comparison_violin_plots/PE_tropho_VEGF.pdf", paper= "special")

plotGeneExpression(data, signaling = "VEGF", split.by = "datasets", color.use = c("grey", "darkblue"))

dev.off()

#define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "PE"

#define a char name used for storing the results of differential expression analysis
features.name = pos.dataset

#perform differential expression analysis
cellchat <- identifyOverExpressedGenes(data, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, 
                                       thresh.pc = 0.1, thresh.fc = 0.1, 
                                       thresh.p= 0.01)


#>Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)

#Extract the ligand-receptor pairs with upregulated ligands & receptors in PE (disease) vs LC (normal)
net.up <- subsetCommunication(cellchat, net = net, datasets = "PE", sources.use= c('vSCT_1', 'vSCT_2', 'vtropho_15'), 
                             targets.use= c('dSMC', 'dVEC', 'dLEC'), 
                              ligand.logFC = 0.1, 
                              receptor.logFC = 0.1)

head(net.up)

write.csv(net.up, file= "./TrophoL_mesenchymal_analysis/vSCT_source_RL_up.csv")

#Only ligand upregulated: denotes case like "Ligand activation" (PE vs late controls)
net.up <- subsetCommunication(cellchat, net = net, datasets = "PE", sources.use= c('vSCT_1', 'vSCT_2', 'vtropho_15'), 
                             targets.use= c('dSMC', 'dVEC', 'dLEC'), 
                              ligand.logFC = 0.1, 
                              receptor.logFC = NULL)

head(net.up)

write.csv(net.up, file= "./TrophoL_mesenchymal_analysis/vSCT_source_ligands_overexp.csv")

#dVEC:: source
net.up <- subsetCommunication(cellchat, net = net, datasets = "PE", targets.use= c('vSCT_1', 'vSCT_2', 'vtropho_15'), 
                             sources.use= 'dVEC', 
                              ligand.logFC = 0.1, 
                              receptor.logFC = 0.1)

head(net.up)

#write.csv(net.up, file= "./TrophoL_mesenchymal_analysis/dVCT_source_ligand_overexp.csv")

#define a positive dataset, i.e., the dataset with positive fold change against the other dataset "LC"
pos.dataset = "LC"

#define a char name used for storing the results of differential expression analysis
features.name = pos.dataset

#perform differential expression analysis:
#thresh.pc: Threshold of the percent of cells expressed in one cluster. 
#thresh.fc: Threshold of Log Fold Change.
#thresh.p: Threshold of p-values
cellchat <- identifyOverExpressedGenes(data, group.dataset = "datasets", pos.dataset = pos.dataset, 
                                       features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, 
                                       thresh.fc = 0.1, thresh.p = 0.05)

#>Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)


#Extract the ligand-receptor pairs with upregulated ligands & receptors in LC (late controls)
net.up <- subsetCommunication(cellchat, net = net, datasets = "LC", sources.use= c('vSCT_1', 'vSCT_2', 'vtropho_15'), 
                             targets.use= c('dSMC', 'dVEC', 'dLEC'), ligand.logFC= NULL,
                             receptor.logFC=0.2)

head(net.up)

write.csv(net.up, file= "./TrophoL_mesenchymal_analysis/dSCTsource_receptor_overexp_LC.csv")



#Read the truncated data: not used in the manuscript (Extra analysis)
data <- readRDS(file= "/data/analysis/preeclampsia_2019/analysis/images/cellphonedb_analysis/TrophoL_mesenchymal_analysis/Truncated_cellchat_comparison_LateC_vs_PE.rds")

data

data@meta$datasets = factor(data@meta$datasets, levels = c("LC", "PE")) # set factor level

#define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "PE"

#define a char name used for storing the results of differential expression analysis
features.name = pos.dataset

#perform differential expression analysis
cellchat <- identifyOverExpressedGenes(data, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, 
                                       thresh.pc = 0.1, thresh.fc = 0.1, 
                                       thresh.p= 0.01)

#>Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)

#Extract the ligand-receptor pairs with upregulated ligands & receptors in PE (disease)
net.up <- subsetCommunication(cellchat, net = net, datasets = "PE", sources.use= c('vSCT_1', 'vSCT_2', 'vtropho_15'), 
                             targets.use= c('dSMC', 'dVEC', 'dLEC'), 
                              ligand.logFC = 0.1, 
                              receptor.logFC = 0.1)

head(net.up)

write.csv(net.up, file= "./TrophoL_mesenchymal_analysis/truncated_comparison/Truncated_vSCTsource_RL_upregulated.csv")

pdf("./TrophoL_mesenchymal_analysis/truncated_comparison/PE_tropho_IL12.pdf", paper= "special")

plotGeneExpression(data, signaling = "IL12", split.by = "datasets", color.use = c("grey", "darkblue"))

dev.off()

pdf("./TrophoL_mesenchymal_analysis/truncated_comparison/PE_tropho_ANGPTL.pdf", paper= "special")

plotGeneExpression(data, signaling = "ANGPTL", split.by = "datasets", color.use = c("grey", "darkblue"))

dev.off()

pdf("./TrophoL_mesenchymal_analysis/truncated_comparison/PE_tropho_SEMA3.pdf", paper= "special")

plotGeneExpression(data, signaling = "SEMA3", split.by = "datasets", color.use = c("grey", "darkblue"))

dev.off()

#Truncated: 

#define a positive dataset, i.e., the dataset with positive fold change against the other dataset "LC"
pos.dataset = "LC"

#define a char name used for storing the results of differential expression analysis
features.name = pos.dataset

#perform differential expression analysis:
#thresh.pc: Threshold of the percent of cells expressed in one cluster. 
#thresh.fc: Threshold of Log Fold Change.
#thresh.p: Threshold of p-values
cellchat <- identifyOverExpressedGenes(data, group.dataset = "datasets", pos.dataset = pos.dataset, 
                                       features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, 
                                       thresh.fc = 0.1, thresh.p = 0.05)

#>Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)


#Extract the ligand-receptor pairs with upregulated ligands & receptors in LC (late controls)
net.up <- subsetCommunication(cellchat, net = net, datasets = "LC", sources.use= c('vSCT_1', 'vSCT_2', 'vtropho_15'), 
                             targets.use= c('dSMC', 'dVEC', 'dLEC'), ligand.logFC= NULL,
                             receptor.logFC=0.1)

head(net.up)

write.csv(net.up, file= "./TrophoL_mesenchymal_analysis/truncated_comparison/trunc_vSCTsource_receptor_overexp_LC.csv")

pdf("./TrophoL_mesenchymal_analysis/truncated_comparison/LC_tropho_VEGF.pdf", paper= "special")

plotGeneExpression(data, signaling = "VEGF", split.by = "datasets", color.use = c("grey", "darkblue"))

dev.off()

pdf("./TrophoL_mesenchymal_analysis/truncated_comparison/LC_tropho_OSM.pdf", paper= "special")

plotGeneExpression(data, signaling = "OSM", split.by = "datasets", color.use = c("grey", "darkblue"))

dev.off()

gg1 <- rankNet(data, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(data, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

pdf("./Trophoblast_mesoendo_ranked_signaling.pdf", paper= "special", h=2)

gg1 <- rankNet(data, mode = "comparison", stacked = T, do.stat = TRUE, color.use = c("grey", "darkblue"))
gg2 <- rankNet(data, mode = "comparison", stacked = F, do.stat = TRUE, color.use = c("grey", "darkblue"))
gg1 + gg2

dev.off() 

sessionInfo()


