Sys.setenv(LANG="en_US.UTF-8")
Sys.setenv(LC_ALL="en_US.UTF-8")

suppressPackageStartupMessages(library(CellChat))

#Read the saved Cellchat object as RDS file:
#To see how the object was created, refer to: Comparative_decidua_immune_RL_analysis_PartI.R 
data <- readRDS(file= "/data/analysis/preeclampsia_2019/analysis/images/cellphonedb_analysis/decidua_immune_analysis/comparative_decidua_immune_analysis.rds")

data

#define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "PE"

#define a char name used for storing the results of differential expression analysis
features.name = pos.dataset

#perform differential expression analysis
cellchat <- identifyOverExpressedGenes(data, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)


#Use the joint cell labels from the merged CellChat object
#Map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)

#Extract the ligand-receptor pairs with upregulated ligands in PE (disease)
net.up <- subsetCommunication(cellchat, net = net, datasets = "PE",ligand.logFC = 0.2, receptor.logFC = NULL)

head(net.up)

#Immune targets: overexpressed_ligands
net.up <- subsetCommunication(cellchat, net = net, datasets = "PE",ligand.logFC = 0.2, receptor.logFC = NULL,
                              sources.use = c('dFB', 'DSC_1', 'dVEC', 'dLEC', 'dSMC'), 
                        targets.use = c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1', 'dTcell', 'dDC', 'dGranulocyte'))

write.csv(net.up, file= "./decidua_immune_analysis/dImmune_targets_overexpressed_ligands.csv")

write.csv(net.up, file= "./decidua_immune_analysis/dImmune_targets_differentialUP.csv")

#Extract the ligand-receptor pairs with upregulated ligands & receptors in PE (disease)
net.up <- subsetCommunication(cellchat, net = net, datasets = "PE",ligand.logFC = 0.2, receptor.logFC = 0.1)

head(net.up)

write.csv(net.up, file= "./decidua_immune_analysis/dImmune_targets_differential_RL_UP.csv")

#Immune targets: when both ligand & receptor exp is up. 
net.up <- subsetCommunication(cellchat, net = net, datasets = "PE",ligand.logFC = 0.2, receptor.logFC = 0.1,
                              sources.use = c('dFB', 'DSC_1', 'dVEC', 'dLEC', 'dSMC'), 
                        targets.use = c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1', 'dTcell', 'dDC', 'dGranulocyte'))

write.csv(net.up, file= "./decidua_immune_analysis/dImmune_targets_differential_RL_UP.csv")

#Immune sources: when both ligand & receptor exp is up. 
net.up <- subsetCommunication(cellchat, net = net, datasets = "PE",ligand.logFC = 0.2, receptor.logFC = 0.1,
                              targets.use = c('dFB', 'DSC_1', 'dVEC', 'dLEC', 'dSMC'), 
                        sources.use = c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1', 'dTcell', 'dDC', 'dGranulocyte'))

write.csv(net.up, file= "./decidua_immune_analysis/dImmune_sources_differential_RL_UP.csv")

#extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in LC, i.e.,downregulated in PE. 
net.down <- subsetCommunication(cellchat, net = net, datasets = "LC",ligand.logFC = -0.1, receptor.logFC = -0.1)

head(net.down)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

#Case-I: Immune targets
sources.use= c('dFB', 'DSC_1', 'dVEC', 'dLEC', 'dSMC')

targets.use= c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1',
                                                            'dTcell', 'dDC', 'dGranulocyte')



#Bubble plots to visualize LR interaction pairs: 
pairLR.use.up = net.up[, "interaction_name", drop = F]

pdf("./decidua_immune_analysis/dImmune_targets_differential_RL_UP.pdf", paper= "special", w=15)

netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c('dFB', 'DSC_1', 'dVEC', 'dLEC', 'dSMC'), 
                        targets.use = c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1', 'dTcell', 'dDC', 'dGranulocyte'), 
                        comparison = c(1, 2),  
                        angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in PE"))

dev.off()

pairLR.use.down = net.down[, "interaction_name", drop = F]

gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, 
                        sources.use = c('dFB', 'DSC_1', 'dVEC', 'dLEC', 'dSMC'), 
                        targets.use = c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1', 'dTcell', 'dDC', 'dGranulocyte'), 
                        comparison = c(1, 2),  angle.x = 90, remove.isolate = T,
                        title.name = paste0("Down-regulated signaling in PE"))

gg2

#cellchat@net$PE

object= cellchat@net



object.list <- list(LC = cellchat@net$LC, PE = cellchat@net$PE)

netVisual_chord_gene(object.list[[2]], net = net.up, lab.cex = 0.8, 
                     small.gap = 3.5, title.name = paste0("Up-regulated signaling in PE"),
                    sources.use = c('dFB', 'DSC_1', 'dVEC', 'dLEC', 'dSMC'), 
                        targets.use = c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1', 'dTcell', 'dDC', 'dGranulocyte'))

cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("LC", "PE")) #set factor level


pdf("./decidua_immune_analysis/Immune_SEMA3.pdf", paper= "special", w=7, h=9)

plotGeneExpression(cellchat, signaling = "SEMA3", split.by = "datasets", color.use = c("grey", "darkblue"))

dev.off() 

cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("LC", "PE")) # set factor level

plotGeneExpression(cellchat, signaling = "ANNEXIN", split.by = "datasets", color.use = c("grey", "darkblue"))

pdf("./decidua_immune_analysis/Immune_FN1.pdf", paper= "special")

plotGeneExpression(cellchat, signaling = "FN1", split.by = "datasets", color.use = c("grey", "darkblue"),
                  colors.ggplot = T)

dev.off()

pdf("./decidua_immune_analysis/Immune_Collagen.pdf", w=7, h=9, paper= "special")

plotGeneExpression(cellchat, signaling = "COLLAGEN", split.by = "datasets", color.use = c("grey", "darkblue"),
                  colors.ggplot = T)

dev.off()

pdf("./decidua_immune_analysis/Immune_ANNEXIN.pdf", paper= "special")

plotGeneExpression(cellchat, signaling = "ANNEXIN", split.by = "datasets", color.use = c("grey", "darkblue"),
                  colors.ggplot = T)

dev.off()

pdf("./decidua_immune_analysis/Immune_Laminin.pdf", w=7, h=9, paper= "special")

plotGeneExpression(cellchat, signaling = "LAMININ", split.by = "datasets", color.use = c("grey", "darkblue"),
                  colors.ggplot = T)

dev.off()

pdf("./decidua_immune_analysis/Immune_MHC1.pdf", paper= "special")

plotGeneExpression(cellchat, signaling = "MHC-I", split.by = "datasets", color.use = c("grey", "darkblue"),
                  colors.ggplot = T)

dev.off()

pdf("./decidua_immune_analysis/Immune_APP.pdf", paper= "special")

plotGeneExpression(cellchat, signaling = "APP", split.by = "datasets", color.use = c("grey", "darkblue"),
                  colors.ggplot = T)

dev.off()

pdf("./decidua_immune_analysis/Immune_NOTCH_PE.pdf", paper= "special")

plotGeneExpression(cellchat, signaling = "NOTCH", split.by = "datasets", color.use = c("grey", "darkblue"),
                  colors.ggplot = T)

dev.off()

pdf("./decidua_immune_analysis/Immune_THBS_PE.pdf", paper= "special")

plotGeneExpression(cellchat, signaling = "THBS", split.by = "datasets", color.use = c("grey", "darkblue"),
                  colors.ggplot = T)

dev.off()

pairLR.use.up = net.up[, "interaction_name", drop = F]

pdf("./decidua_immune_analysis/dImmune_source_differential_RL_UP.pdf", paper= "special", w=10)

netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, targets.use = c('dFB', 'DSC_1', 'dVEC', 'dLEC', 'dSMC'), 
                        sources.use = c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1', 'dTcell', 'dDC', 'dGranulocyte'), 
                        comparison = c(1, 2),  
                        angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in PE"))

dev.off()

#extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in LC, i.e.,downregulated in PE. 
#Use logFC= -0.1 
net.down <- subsetCommunication(cellchat, net = net, datasets = "LC",ligand.logFC = -0.1, receptor.logFC = -0.1)

head(net.down)

write.csv(net.down, file= "./decidua_immune_analysis/dImmune_differential_downregulated_PE.csv")



#Extract the ligand-receptor pairs with upregulated receptors in PE (disease)
#Case: only receptor is up (so, include cases of both Ligand activation & ligand starvation)
net.up <- subsetCommunication(cellchat, net = net, datasets = "PE",ligand.logFC = NULL, receptor.logFC = 0.2)

head(net.up)

#Immune targets: overexpressed_receptors
net.up <- subsetCommunication(cellchat, net = net, datasets = "PE",ligand.logFC = NULL, receptor.logFC = 0.2,
                              sources.use = c('dFB', 'DSC_1', 'dVEC', 'dLEC', 'dSMC'), 
                        targets.use = c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1', 'dTcell', 'dDC', 'dGranulocyte'))

write.csv(net.up, file= "./decidua_immune_analysis/dImmune_targets_overexpressed_receptors.csv")

#Immune sources: overexpressed_receptors
net.up <- subsetCommunication(cellchat, net = net, datasets = "PE",ligand.logFC = NULL, receptor.logFC = 0.2,
                              targets.use = c('dFB', 'DSC_1', 'dVEC', 'dLEC', 'dSMC'), 
                        sources.use = c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1', 'dTcell', 'dDC', 'dGranulocyte'))

write.csv(net.up, file= "./decidua_immune_analysis/dImmune_sources_overexpressed_receptors.csv")

#Extract the ligand-receptor pairs with upregulated ligands & receptors in PE (disease): Edge activation 
net.up <- subsetCommunication(cellchat, net = net, datasets = "PE",ligand.logFC = 0.2, receptor.logFC = 0.2)

head(net.up)

#Immune targets: overexpressed_both (Edge activation)
net.up <- subsetCommunication(cellchat, net = net, datasets = "PE",ligand.logFC = 0.2, receptor.logFC = 0.2,
                              sources.use = c('dFB', 'DSC_1', 'dVEC', 'dLEC', 'dSMC'), 
                        targets.use = c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1', 'dTcell', 'dDC', 'dGranulocyte'))

write.csv(net.up, file= "./decidua_immune_analysis/dImmune_targets_overexpressed_both.csv")

#Immune targets: overexpressed_both (Edge activation); reduce receptor.logFC to 0.1 
net.up <- subsetCommunication(cellchat, net = net, datasets = "PE",ligand.logFC = 0.2, receptor.logFC = 0.1,
                              sources.use = c('dFB', 'DSC_1', 'dVEC', 'dLEC', 'dSMC'), 
                        targets.use = c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1', 'dTcell', 'dDC', 'dGranulocyte'))

write.csv(net.up, file= "./decidua_immune_analysis/dImmune_targets_overexpressed_both_lowcutoff.csv")

#Immune sources: overexpressed_both (Edge activation)
net.up <- subsetCommunication(cellchat, net = net, datasets = "PE",ligand.logFC = 0.2, receptor.logFC = 0.2,
                              targets.use = c('dFB', 'DSC_1', 'dVEC', 'dLEC', 'dSMC'), 
                        sources.use = c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1', 'dTcell', 'dDC', 'dGranulocyte'))

write.csv(net.up, file= "./decidua_immune_analysis/dImmune_sources_overexpressed_both.csv")

#Immune sources: overexpressed_both (Edge activation); reduce receptor.logFC to 0.1 
net.up <- subsetCommunication(cellchat, net = net, datasets = "PE",ligand.logFC = 0.2, receptor.logFC = 0.1,
                              targets.use = c('dFB', 'DSC_1', 'dVEC', 'dLEC', 'dSMC'), 
                        sources.use = c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1', 'dTcell', 'dDC', 'dGranulocyte'))

write.csv(net.up, file= "./decidua_immune_analysis/dImmune_sources_overexpressed_both_lowercutoff.csv")



#extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in LC, i.e.,downregulated in PE. 
net.down <- subsetCommunication(cellchat, net = net, datasets = "LC",ligand.logFC = -0.1, receptor.logFC = -0.1)

head(net.down)

#Edge deactivation: both receptor/ligands are down. 
#Immune targets (receivers)
net.down <- subsetCommunication(cellchat, net = net, datasets = "LC",ligand.logFC = -0.1, receptor.logFC = -0.1,
                                sources.use = c('dFB', 'DSC_1', 'dVEC', 'dLEC', 'dSMC'), 
                        targets.use = c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1', 'dTcell', 'dDC', 'dGranulocyte'))

write.csv(net.up, file= "./decidua_immune_analysis/dImmune_targets_downregulated_both_lowcutoff.csv")

#Edge deactivation: both receptor/ligands are down. 
#Immune source (senders)
net.down <- subsetCommunication(cellchat, net = net, datasets = "LC",ligand.logFC = -0.1, receptor.logFC = -0.1,
                                targets.use = c('dFB', 'DSC_1', 'dVEC', 'dLEC', 'dSMC'), 
                        sources.use = c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1', 'dTcell', 'dDC', 'dGranulocyte'))

write.csv(net.down, file= "./decidua_immune_analysis/dImmune_sources_downregulated_both_lowcutoff.csv")



#define a positive dataset, i.e., the dataset with positive fold change against the other dataset
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



#Extract the ligand-receptor pairs with upregulated ligands & receptors in PE (disease): logFC = 0.1
#Immune cell-types can be both senders & receivers. 
net.up <- subsetCommunication(cellchat, net = net, datasets = "LC",ligand.logFC = 0.1, receptor.logFC = 0.1)

head(net.up)

write.csv(net.up, file= "./decidua_immune_analysis/dImmune_overexp_LC.csv")

net.up <- subsetCommunication(cellchat, net = net, datasets = "LC", ligand.logFC = 0.1, receptor.logFC = NULL, 
                                sources.use = c('dFB', 'DSC_1', 'dVEC', 'dLEC', 'dSMC'), 
                        targets.use = c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1', 'dTcell', 'dDC', 'dGranulocyte'))
head(net.up)

write.csv(net.up, file= "./decidua_immune_analysis/dImmune_targets_overexp_ligands_LC.csv")

net.up <- subsetCommunication(cellchat, net = net, datasets = "LC", ligand.logFC = 0.1, receptor.logFC = NULL, 
                                targets.use = c('dFB', 'DSC_1', 'dVEC', 'dLEC', 'dSMC'), 
                        sources.use = c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1', 'dTcell', 'dDC', 'dGranulocyte'))
head(net.up)

write.csv(net.up, file= "./decidua_immune_analysis/dImmune_sources_overexp_ligands_LC.csv")



#Overexpressed receptors:: LC 
#Source: immune. 
net.up <- subsetCommunication(cellchat, net = net, datasets = "LC", ligand.logFC = NULL, receptor.logFC = 0.1, 
                                targets.use = c('dFB', 'DSC_1', 'dVEC', 'dLEC', 'dSMC'), 
                        sources.use = c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1', 'dTcell', 'dDC', 'dGranulocyte'))
head(net.up)

write.csv(net.up, file= "./decidua_immune_analysis/dImmune_sources_overexp_receptors_LC.csv")

#Overexpressed receptors:: LC 
#Source: immune. 
net.up <- subsetCommunication(cellchat, net = net, datasets = "LC", ligand.logFC = NULL, receptor.logFC = 0.1, 
                                sources.use = c('dFB', 'DSC_1', 'dVEC', 'dLEC', 'dSMC'), 
                        targets.use = c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1', 'dTcell', 'dDC', 'dGranulocyte'))
head(net.up)

write.csv(net.up, file= "./decidua_immune_analysis/dImmune_targets_overexp_receptors_LC.csv")



pdf("./decidua_immune_analysis/Immune_PARs.pdf", paper= "special")

plotGeneExpression(cellchat, signaling = "PARs", split.by = "datasets", color.use = c("grey", "darkblue"),
                  colors.ggplot = T)

dev.off()

pdf("./decidua_immune_analysis/Immune_IGF.pdf", paper= "special")

plotGeneExpression(cellchat, signaling = "IGF", split.by = "datasets", color.use = c("grey", "darkblue"),
                  colors.ggplot = T)

dev.off()

pdf("./decidua_immune_analysis/Immune_SPP1.pdf", paper= "special")

plotGeneExpression(cellchat, signaling = "SPP1", split.by = "datasets", color.use = c("grey", "darkblue"),
                  colors.ggplot = T)

dev.off()

pdf("./decidua_immune_analysis/dImmune_ranked_signaling.pdf", paper= "special")

gg1 <- rankNet(data, mode = "comparison", stacked = T, do.stat = TRUE, color.use = c("grey", "darkblue"))
gg2 <- rankNet(data, mode = "comparison", stacked = F, do.stat = TRUE, color.use = c("grey", "darkblue"))
gg1 + gg2

dev.off() 


