Sys.setenv(LANG="en_US.UTF-8")
Sys.setenv(LC_ALL="en_US.UTF-8")

suppressPackageStartupMessages(library(CellChat))

#Read the villi immune cellchat object: saved as RDS file
#To see how the RDS is created, refer to: Comparative_villi_immune_RL_analysis_PartI.R 
data <- readRDS(file= "/data/analysis/preeclampsia_2019/analysis/images/cellphonedb_analysis/villi_immune_cellchat/Villi_immune_LateC_vs_PE.rds")

data



#define a positive dataset, i.e., the dataset (PE; disease) with positive fold change against the other dataset (LC; late gestation)
pos.dataset = "PE"

#define a char name used for storing the results of differential expression analysis
features.name = pos.dataset

#perform differential expression analysis
cellchat <- identifyOverExpressedGenes(data, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, 
                                       thresh.pc = 0, 
                                       thresh.fc = 0, 
                                       thresh.p = 0.05)


net <- netMappingDEG(cellchat, features.name = features.name)

#Extract the ligand-receptor pairs with upregulated ligands & up/down receptors in PE (disease): Ligand activation 
net.up <- subsetCommunication(cellchat, net = net, datasets = "PE", ligand.logFC = 0.1, 
                              receptor.logFC= NULL,
                             targets.use = c('vTcell', 'vHBC', 'vDC'), 
                        sources.use= c('vVEC', 'vMC', 'vFB', 'vVCT'))

head(net.up)

write.csv(net.up, file= "./villi_immune_cellchat/differential_interactions/vImmune_targets_ligands_UP.csv")

net <- netMappingDEG(cellchat, features.name = features.name)

#Extract the ligand-receptor pairs with upregulated receptors in PE (disease) with up/down ligands. 
net.up <- subsetCommunication(cellchat, net = net, datasets = "PE", ligand.logFC = NULL, 
                              receptor.logFC= 0.1,
                             targets.use = c('vTcell', 'vHBC', 'vDC'), 
                        sources.use= c('vVEC', 'vMC', 'vFB', 'vVCT'))

head(net.up)

write.csv(net.up, file= "./villi_immune_cellchat/differential_interactions/vImmune_targets_receptors_UP.csv")

#Use the joint cell labels from the merged CellChat object
#Map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)

#Extract the ligand-receptor pairs with upregulated ligands & receptors in PE (disease)
net.up <- subsetCommunication(cellchat, net = net, datasets = "PE",ligand.logFC = 0.1, receptor.logFC= 0.1)

head(net.up)

write.csv(net.up, file= "./villi_immune_cellchat/vImmune_differential_RL_UP.csv")





#Extract the ligand-receptor pairs with upregulated ligands & receptors in PE (disease)
#Use logFC 0.1 for each receptor & ligand. Use the "immune cells" as "targets"/receivers. 
net.up <- subsetCommunication(cellchat, net = net, datasets = "PE", ligand.logFC = 0.1, receptor.logFC= 0.1,
                             targets.use = c('vTcell', 'vHBC', 'vDC'), 
                        sources.use= c('vVEC', 'vMC', 'vFB', 'vVCT'))

head(net.up)
write.csv(net.up, file= "./villi_immune_cellchat/differential_interactions/vImmune_targets_differential_RL_UP.csv")



#Extract the ligand-receptor pairs with upregulated receptors in PE (disease); receptor.logFC= 0.1
net.up <- subsetCommunication(cellchat, net = net, datasets = "PE", ligand.logFC = NULL, receptor.logFC= 0.1,
                             sources.use = c('vTcell', 'vHBC', 'vDC'), 
                        targets.use= c('vVEC', 'vMC', 'vFB', 'vVCT'))

head(net.up)
write.csv(net.up, file= "./villi_immune_cellchat/differential_interactions/vImmune_source_receptors_up.csv")

#Extract the ligand-receptor pairs with upregulated ligands in PE (disease)
#Use the immune cells as "source" 
net.up <- subsetCommunication(cellchat, net = net, datasets = "PE", ligand.logFC = 0.1, receptor.logFC= NULL,
                             sources.use = c('vTcell', 'vHBC', 'vDC'), 
                        targets.use= c('vVEC', 'vMC', 'vFB', 'vVCT'))

head(net.up)

write.csv(net.up, file= "./villi_immune_cellchat/differential_interactions/vImmune_source_ligands_up.csv")

#Extract the ligand-receptor pairs with upregulated ligands & receptors in PE (disease)
net.up <- subsetCommunication(cellchat, net = net, datasets = "PE", ligand.logFC = 0.1, receptor.logFC= 0.1,
                             sources.use = c('vTcell', 'vHBC', 'vDC'), 
                        targets.use= c('vVEC', 'vMC', 'vFB', 'vVCT'))

head(net.up)

write.csv(net.up, file= "./villi_immune_cellchat/differential_interactions/vImmune_source_differential_RL_UP.csv")



gene.up <- extractGeneSubsetFromPair(net.up, cellchat)

gene.up

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c('vTcell', 'vHBC', 'vDC'), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in PE"))
#> Comparing communications on a merged object
gg1

pairLR.use.up = net.up[, "interaction_name", drop = F]

pdf("./villi_immune_cellchat/villi_Immune_targets_differential_RL_UP.pdf", paper= "special")

netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, 
                        targets.use = c('vTcell', 'vHBC', 'vDC'), 
                        sources.use= c('vVEC', 'vMC', 'vFB', 'vVCT'),
                        comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in PE"))
#> Comparing communications on a merged object
dev.off()

cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("LC", "PE")) # set factor level

pdf("./villi_immune_cellchat/vImmune_THBS.pdf", paper= "special")

plotGeneExpression(cellchat, signaling = "THBS", split.by = "datasets", color.use = c("grey", "darkblue"))

dev.off()

pdf("./villi_immune_cellchat/vImmune_GDF.pdf", paper= "special")

plotGeneExpression(cellchat, signaling = "GDF", split.by = "datasets", color.use = c("grey", "darkblue"))

dev.off()

pdf("./villi_immune_cellchat/vImmune_VISFATIN.pdf", paper= "special")

plotGeneExpression(cellchat, signaling = "VISFATIN", split.by = "datasets", color.use = c("grey", "darkblue"))

dev.off()

pdf("./villi_immune_cellchat/vImmune_BMP.pdf", paper= "special", h=10)

plotGeneExpression(cellchat, signaling = "BMP", split.by = "datasets", color.use = c("grey", "darkblue"))

dev.off()

pdf("./villi_immune_cellchat/vImmune_PDGF.pdf", paper= "special")

plotGeneExpression(cellchat, signaling = "PDGF", split.by = "datasets", color.use = c("grey", "darkblue"))

dev.off()

pdf("./villi_immune_cellchat/vImmune_Collagen.pdf", paper= "special", h=10)

plotGeneExpression(cellchat, signaling = "COLLAGEN", split.by = "datasets", color.use = c("grey", "darkblue"))

dev.off()

pdf("./villi_immune_cellchat/vImmune_Laminin.pdf", paper= "special", h=10)

plotGeneExpression(cellchat, signaling = "LAMININ", split.by = "datasets", color.use = c("grey", "darkblue"))

dev.off()

#Compare information flow: when immune cell types are the "targets" or "ligand receivers"

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, targets.use = c('vTcell', 'vHBC', 'vDC'), 
                        sources.use= c('vVEC', 'vMC', 'vFB', 'vVCT'))

gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, targets.use = c('vTcell', 'vHBC', 'vDC'), 
                        sources.use= c('vVEC', 'vMC', 'vFB', 'vVCT'))

gg1 + gg2

#Compare information flow: when immune cell types are the "source" or "ligand senders"
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, sources.use = c('vTcell', 'vHBC', 'vDC'), 
                        targets.use= c('vVEC', 'vMC', 'vFB', 'vVCT'))

gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, sources.use = c('vTcell', 'vHBC', 'vDC'), 
                        targets.use= c('vVEC', 'vMC', 'vFB', 'vVCT'))

gg1 + gg2

pdf("./villi_immune_cellchat/vImmune_SPP1.pdf", paper= "special")

plotGeneExpression(cellchat, signaling = "SPP1", split.by = "datasets", color.use = c("grey", "darkblue"))

dev.off()

pdf("./villi_immune_cellchat/vImmune_CD46.pdf", paper= "special")

plotGeneExpression(cellchat, signaling = "CD46", split.by = "datasets", color.use = c("grey", "darkblue"))

dev.off()

pdf("./villi_immune_cellchat/vImmune_IGF.pdf", paper= "special")

plotGeneExpression(cellchat, signaling = "IGF", split.by = "datasets", color.use = c("grey", "darkblue"))

dev.off()



#define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "LC"

#define a char name used for storing the results of differential expression analysis
features.name = pos.dataset

#perform differential expression analysis
cellchat <- identifyOverExpressedGenes(data, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, 
                                       thresh.p = 0.05)

#>Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)

#Extract the ligand-receptor pairs with upregulated ligands & receptors in LC (controls)
net.up <- subsetCommunication(cellchat, net = net, datasets = "LC",ligand.logFC = 0.1, receptor.logFC= 0.1)

head(net.up)

#Extract the ligand-receptor pairs with upregulated ligands & receptors in LC (controls)
#Immune source: edge activation 
net.up <- subsetCommunication(cellchat, net = net, datasets = "LC", ligand.logFC = 0.1, receptor.logFC= 0.1,
                             sources.use = c('vTcell', 'vHBC', 'vDC'), 
                        targets.use= c('vVEC', 'vMC', 'vFB', 'vVCT'))

head(net.up)

write.csv(net.up, file= "./villi_immune_cellchat/differential_interactions/LC_vImmune_source_differential_RL_UP.csv")

#Extract the ligand-receptor pairs with upregulated ligands & receptors in PE (disease): immune targets
#Edge activation 
net.up <- subsetCommunication(cellchat, net = net, datasets = "LC", ligand.logFC = 0.1, receptor.logFC= 0.1,
                             targets.use = c('vTcell', 'vHBC', 'vDC'), 
                        sources.use= c('vVEC', 'vMC', 'vFB', 'vVCT'))

head(net.up)

write.csv(net.up, file= "./villi_immune_cellchat/differential_interactions/LC_vImmune_targets_differential_RL_UP.csv")

#Extract the ligand-receptor pairs with upregulated ligands)
net.up <- subsetCommunication(cellchat, net = net, datasets = "LC",ligand.logFC = NULL, receptor.logFC= 0.1,
                             targets.use = c('vTcell', 'vHBC', 'vDC'), 
                        sources.use= c('vVEC', 'vMC', 'vFB', 'vVCT'))

head(net.up)

#write.csv(net.up, file= "./villi_immune_cellchat/differential_interactions/LC_vImmune_targets_differential_ligands_up.csv")

net.up <- subsetCommunication(cellchat, net = net, datasets = "LC",ligand.logFC = 0.1, receptor.logFC= NULL,
                             targets.use = c('vTcell', 'vHBC', 'vDC'), 
                        sources.use= c('vVEC', 'vMC', 'vFB', 'vVCT'))

head(net.up)

#Extract the ligand-receptor pairs with upregulated ligands & receptors in LC (controls)
net.up <- subsetCommunication(cellchat, net = net, datasets = "LC")
head(net.up)
write.csv(net.up, file= "./villi_immune_cellchat/differential_interactions/LC_differential_all.csv")



Sys.setenv(LANG="en_US.UTF-8")
Sys.setenv(LC_ALL="en_US.UTF-8")

suppressPackageStartupMessages(library(CellChat))

data <- readRDS(file= "/data/analysis/preeclampsia_2019/analysis/images/cellphonedb_analysis/villi_immune_cellchat/Villi_immune_LateC_vs_PE.rds")

data

#data@net$LC

object.list <- list(LC = data@netP$LC, PE = data@netP$PE)



pdf("./villi_immune_cellchat/vImmune_ranked_signaling.pdf", paper= "special")

gg1 <- rankNet(data, mode = "comparison", stacked = T, do.stat = TRUE, color.use = c("grey", "darkblue"))
gg2 <- rankNet(data, mode = "comparison", stacked = F, do.stat = TRUE, color.use = c("grey", "darkblue"))
gg1 + gg2

dev.off() 


