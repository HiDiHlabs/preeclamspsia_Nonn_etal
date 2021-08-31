Sys.setenv(LANG="en_US.UTF-8")
Sys.setenv(LC_ALL="en_US.UTF-8")

suppressPackageStartupMessages(library(Seurat))

#Read the Seurat object: 
data <- readRDS(file = "/data/analysis/preeclampsia_2019/analysis/images/PE_markers_sexcorrected/updated_placenta_0704.rds")

data

table(Idents(data))

#Subset only decidua. 
Idents(data) <- "tissue"

table(Idents(data))


subset = "Decidua" 
seurat_dec = subset(data, idents = subset) 

#seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["group"]])
seurat_dec
#Change the "idents" back to cell type:
Idents(seurat_dec) <- "cell_type_semifinal_v2"
table(Idents(seurat_dec))

seurat_dec



library(plyr)

#Merge both dFB_1 & dFB_2 into one:
seurat_dec <- RenameIdents(object = seurat_dec,  'dFB_1' = 'dFB', 'dFB_2' = 'dFB')

table(Idents(seurat_dec))



#Subset necessary cell-types from decidua & selectively remove them with invert=TRUE. 
##Subset necessary cell-types from decidua: dLECP/dLEC_dysfunctional is removed since it's not affected by PE in terms of DEG. 
subset = c("dNK_prol", "dPlasmaCell", "DSC_2", "dLEC_dysfunctional", "dSCT", "dMonocyte", "dEVT", "dEpC", "dMSC") 

seurat_chat = subset(seurat_dec, idents = subset, invert=TRUE) 

#seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["group"]])
seurat_chat
table(Idents(seurat_chat))



#Subset by Early, Late controls & PE samples: "group" metadata. 
Idents(seurat_chat) <- "group"

table(Idents(seurat_chat))

suppressPackageStartupMessages(library(CellChat))

#Subset "late controls" only: Late_Decidua_C. 
subset = c("Late_Decidua_C")
seurat_late = subset(seurat_chat, idents = subset) 
table(Idents(seurat_late))


#Change the "idents" back to cell type:
Idents(seurat_late) <- "cell_type_semifinal_v2"
#Merge both dFB_1 & dFB_2 into one: mainly to avoid class imbalance. 
seurat_late <- RenameIdents(object = seurat_late,  'dFB_1' = 'dFB', 'dFB_2' = 'dFB')
table(Idents(seurat_late))


#Get normalized data from the Seurat object: 
lateC.input <- GetAssayData(seurat_late, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(seurat_late)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels

#Create a CellChat object using data matrix as input: (from Seurat)
lateC.cellchat <- createCellChat(object = lateC.input, meta = meta, group.by = "group")

lateC.cellchat #LateC Cell chat object.

groupSize <- as.numeric(table(lateC.cellchat@idents)) # number of cells in each cell group
groupSize

CellChatDB <- CellChatDB.human 
#CellChatDB.use <- CellChatDB # simply use the default CellChatDB
CellChatDB.use <- CellChatDB
lateC.cellchat@DB <- CellChatDB.use

#Preprocessing & inferring interactions: LateC Cell-Chat. 
lateC.cellchat <- subsetData(lateC.cellchat) #subset the expression data of signaling genes for saving computation cost

#Identify over-expressed signaling genes associated with each cell group: 
lateC.cellchat <- identifyOverExpressedGenes(lateC.cellchat)

#Identify over-expressed ligand-receptor interactions (pairs) within the used CellChatDB: 
lateC.cellchat <- identifyOverExpressedInteractions(lateC.cellchat)
lateC.cellchat <- projectData(lateC.cellchat, PPI.human)

#Compute the communication probability and infer cellular communication network:
lateC.cellchat <- computeCommunProb(lateC.cellchat, raw.use = TRUE, type= "triMean")

#Filter out the cell-cell communication if there are only few number of cells in certain cell groups
lateC.cellchat <- filterCommunication(lateC.cellchat, min.cells = 10)

#Save all data. 
df.net <- subsetCommunication(lateC.cellchat)

head(df.net)

write.csv(df.net, file= "./decidua_immune_analysis/LateC_dImmune_all_interactions_triMean.csv")

#Extract the inferred cellular communication network as a data frame:

#Source: Immune. 
df.net <- subsetCommunication(lateC.cellchat, sources.use= c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1',
                                                            'dTcell', 'dDC', 'dGranulocyte'))

head(df.net)

write.csv(df.net, file= "./decidua_immune_analysis/LateC_dImmuneLig_allTargets_triMean.csv")

#Targets: "immune"
df.net <- subsetCommunication(lateC.cellchat, targets.use= c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1',
                                                            'dTcell', 'dDC', 'dGranulocyte'))

head(df.net)

write.csv(df.net, file= "./decidua_immune_analysis/LateC_dImmuneTargets_allLig_triMean.csv")

#Infer the cell-cell communication at a signaling pathway level:
lateC.cellchat <- computeCommunProbPathway(lateC.cellchat)

#Calculate the aggregated cell-cell communication network:
#We can calculate the aggregated cell-cell communication network by counting 
#the number of links or summarizing the communication probability:
lateC.cellchat <- aggregateNet(lateC.cellchat)

#Access all the signaling pathways showing significant communications
pathways.show.all <- lateC.cellchat@netP$pathways

pathways.show.all

#Targets: "immune" but mesoenothelial sources
df.net <- subsetCommunication(lateC.cellchat, sources.use= c('dFB', 'DSC_1', 'dVEC', 'dLEC', 'dSMC'), 
                              targets.use= c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1',
                                                            'dTcell', 'dDC', 'dGranulocyte'))

head(df.net)

write.csv(df.net, file= "./decidua_immune_analysis/LateC_dImmuneTargets_mesoendoSource_triMean.csv")

#Sources: "immune" but mesoenothelial targets
df.net <- subsetCommunication(lateC.cellchat, targets.use= c('dFB', 'DSC_1', 'dVEC', 'dLEC', 'dSMC'), 
                              sources.use= c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1',
                                                            'dTcell', 'dDC', 'dGranulocyte'))

head(df.net)

write.csv(df.net, file= "./decidua_immune_analysis/LateC_dImmuneSource_mesoendoTargets_triMean.csv")



pathways.show <- lateC.cellchat@netP$pathways

#Violin-plots: PE. 
for (i in 1:length(pathways.show)) {
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- plotGeneExpression(lateC.cellchat, signaling = pathways.show[i])
  ggsave(filename=paste0(pathways.show[i], "_LR_contribution_Violin_lateC.pdf"), plot=gg, dpi = 500)
}

groupSize <- as.numeric(table(lateC.cellchat@idents))

netVisual_circle(lateC.cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Number of interactions", targets.use= c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1',
                                                            'dTcell', 'dDC', 'dGranulocyte'),
                vertex.label.cex = 0.6, )

#Define both source & target: 
groupSize <- as.numeric(table(lateC.cellchat@idents))

netVisual_circle(lateC.cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Number of interactions", targets.use= c('dFB', 'DSC_1', 'dVEC', 'dLEC', 'dSMC'), 
                 sources.use= c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1',
                                                            'dTcell', 'dDC', 'dGranulocyte'),
                vertex.label.cex = 0.6)

#Define both source & target: 
groupSize <- as.numeric(table(lateC.cellchat@idents))

netVisual_circle(lateC.cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Number of interactions", sources.use= c('dFB', 'DSC_1', 'dVEC', 'dLEC', 'dSMC'), 
                 targets.use= c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1',
                                                            'dTcell', 'dDC', 'dGranulocyte'),
                vertex.label.cex = 0.6)

pathways.show <- pathways.show.all

lateC.cellchat <- netAnalysis_computeCentrality(lateC.cellchat, slot.name = "netP")

#Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(lateC.cellchat, pattern = "outgoing", width = 6, height = 12, color.heatmap = "GnBu")
ht2 <- netAnalysis_signalingRole_heatmap(lateC.cellchat, pattern = "incoming", width = 6, height = 12, color.heatmap = "GnBu")
ht1 + ht2



seurat_dec

table(Idents(seurat_chat))

#Subset "PE" only: 
subset = c("Late_Decidua_PE")
seurat_pe = subset(seurat_chat, idents = subset) 
table(Idents(seurat_pe))

#Change the "idents" back to cell type:
Idents(seurat_pe) <- "cell_type_semifinal_v2"
seurat_pe <- RenameIdents(object = seurat_pe,  'dFB_1' = 'dFB', 'dFB_2' = 'dFB')
table(Idents(seurat_pe))

#Analyze data for PE dataset: from Seurat 
pe.input <- GetAssayData(seurat_pe, assay = "RNA", slot = "data") # normalized data matrix

labels <- Idents(seurat_pe)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels

#Create a CellChat object using data matrix as input: (from Seurat)
pe.cellchat <- createCellChat(object = pe.input, meta = meta, group.by = "group")

pe.cellchat #PE Cell chat object. 
levels(pe.cellchat@idents) # show factor levels of the cell labels

groupSize <- as.numeric(table(pe.cellchat@idents)) # number of cells in each cell group

groupSize #PE cells per group

#set the used database in the object
pe.cellchat@DB <- CellChatDB.use

pe.cellchat <- subsetData(pe.cellchat) #subset the expression data of signaling genes for saving computation cost

pe.cellchat <- identifyOverExpressedGenes(pe.cellchat)
pe.cellchat <- identifyOverExpressedInteractions(pe.cellchat)
pe.cellchat <- projectData(pe.cellchat, PPI.human)

#Compute the communication probability and infer cellular communication network:
pe.cellchat <- computeCommunProb(pe.cellchat, raw.use = TRUE, type= "triMean")

#Filter out the cell-cell communication if there are only few number of cells in certain cell groups
pe.cellchat <- filterCommunication(pe.cellchat, min.cells = 10)

#Extract the inferred cellular communication network as a data frame:
#We provide a function subsetCommunication to easily access the inferred cell-cell communications of interest. For example,
df.net.pe <- subsetCommunication(pe.cellchat)

head(df.net.pe)

write.csv(df.net.pe, file= "./decidua_immune_analysis/PE_analysis/PE_dImmune_all_triMean.csv")

#Sources: "immune" but mesoenothelial targets
df.net <- subsetCommunication(pe.cellchat, targets.use= c('dFB', 'DSC_1', 'dVEC', 'dLEC', 'dSMC'), 
                              sources.use= c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1',
                                                            'dTcell', 'dDC', 'dGranulocyte'))

head(df.net)

write.csv(df.net, file= "./decidua_immune_analysis/PE_analysis/PE_dImmuneSource_mesoendoTargets_triMean.csv")

#targets: "immune" but mesoenothelial sources. 
df.net <- subsetCommunication(pe.cellchat, sources.use= c('dFB', 'DSC_1', 'dVEC', 'dLEC', 'dSMC'), 
                              targets.use= c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1',
                                                            'dTcell', 'dDC', 'dGranulocyte'))

head(df.net)

write.csv(df.net, file= "./decidua_immune_analysis/PE_analysis/PE_dImmuneTargets_mesoendoSources_triMean.csv")



#Infer the cell-cell communication at a signaling pathway level:
pe.cellchat <- computeCommunProbPathway(pe.cellchat)

#Calculate the aggregated cell-cell communication network:
#We can calculate the aggregated cell-cell communication network by counting 
#the number of links or summarizing the communication probability:
pe.cellchat <- aggregateNet(pe.cellchat)

#Access all the signaling pathways showing significant communications
pathways.show.pe <- pe.cellchat@netP$pathways

pathways.show.pe

pathways.show <- pathways.show.pe

pe.cellchat <- netAnalysis_computeCentrality(pe.cellchat, slot.name = "netP")

#Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(pe.cellchat, pattern = "outgoing", width = 6, height = 12, color.heatmap = "GnBu")
ht2 <- netAnalysis_signalingRole_heatmap(pe.cellchat, pattern = "incoming", width = 6, height = 12, color.heatmap = "GnBu")
ht1 + ht2

#Define both source & target: 
groupSize <- as.numeric(table(pe.cellchat@idents))

netVisual_circle(pe.cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Number of interactions", sources.use= c('dFB', 'DSC_1', 'dVEC', 'dLEC', 'dSMC'), 
                 targets.use= c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1',
                                                            'dTcell', 'dDC', 'dGranulocyte'),
                vertex.label.cex = 0.6)

groupSize <- as.numeric(table(pe.cellchat@idents))

netVisual_circle(pe.cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Number of interactions", targets.use= c('dFB', 'DSC_1', 'dVEC', 'dLEC', 'dSMC'), 
                 sources.use= c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1',
                                                            'dTcell', 'dDC', 'dGranulocyte'),
                vertex.label.cex = 0.6)

library(patchwork)

#Merge LateC & PE lists: 
object.list <- list(LC = lateC.cellchat, PE = pe.cellchat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

cellchat #combined

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

#decidua immune: targets
netVisual_diffInteraction(cellchat, weight.scale = T, vertex.weight=10, vertex.label.cex= 0.6,
                         edge.width.max= 5, alpha.edge = 0.8, vertex.size.max = 8)

#decidua immune: targets
netVisual_diffInteraction(cellchat, weight.scale = T, vertex.weight=10, vertex.label.cex= 0.6,
                         edge.width.max= 5, alpha.edge = 0.8, vertex.size.max = 8,
                         targets.use= c('dFB', 'DSC_1', 'dVEC', 'dLEC', 'dSMC'), 
                         sources.use= c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1',
                                                            'dTcell', 'dDC', 'dGranulocyte'))

netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]),
                 targets.use= c('dFB', 'DSC_1', 'dVEC', 'dLEC', 'dSMC'), 
                         sources.use= c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1',
                                                            'dTcell', 'dDC', 'dGranulocyte'))
}

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]),
                 sources.use= c('dFB', 'DSC_1', 'dVEC', 'dLEC', 'dSMC'), 
                         targets.use= c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1',
                                                            'dTcell', 'dDC', 'dGranulocyte'))
}

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2



#Immune targets: 
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, sources.use= c('dFB', 'DSC_1', 'dVEC', 'dLEC', 'dSMC'), 
                         targets.use= c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1',
                                                            'dTcell', 'dDC', 'dGranulocyte'))

gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, sources.use= c('dFB', 'DSC_1', 'dVEC', 'dLEC', 'dSMC'), 
                         targets.use= c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1',
                                                            'dTcell', 'dDC', 'dGranulocyte'))

gg1 + gg2

#Immune source: 
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, targets.use= c('dFB', 'DSC_1', 'dVEC', 'dLEC', 'dSMC'), 
                         sources.use= c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1',
                                                            'dTcell', 'dDC', 'dGranulocyte'))

gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, targets.use= c('dFB', 'DSC_1', 'dVEC', 'dLEC', 'dSMC'), 
                         sources.use= c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1',
                                                            'dTcell', 'dDC', 'dGranulocyte'))

gg1 + gg2

#cellchat@netP

library(ComplexHeatmap)

i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 6, height = 10)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 6, height = 10)
draw(ht1 + ht2, ht_gap = unit(1, "cm"))

#Overall signaling (incoming plus outgoing): LC vs PE (comparison)
#Summarized at the level of pathways. 

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 10, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 10, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 6, height = 10, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 6, height = 10, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(1, "cm"))

#Identify the upgulated and down-regulated signaling ligand-receptor pairs. 
pdf("./decidua_immune_analysis/dImmune_targets.pdf", paper= "special", w=15, h=13)

netVisual_bubble(cellchat, sources.use= c('dFB', 'DSC_1', 'dVEC', 'dLEC', 'dSMC'), 
                         targets.use= c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1', 'dTcell', 'dDC', 'dGranulocyte'), 
                 comparison = c(1, 2), angle.x = 45)
#> Comparing communications on a merged object

dev.off()

pdf("./decidua_immune_analysis/dImmune_sources.pdf", paper= "special", w=15)

netVisual_bubble(cellchat, targets.use= c('dFB', 'DSC_1', 'dVEC', 'dLEC', 'dSMC'), 
                         sources.use= c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1', 'dTcell', 'dDC', 'dGranulocyte'), 
                 comparison = c(1, 2), angle.x = 45)
#> Comparing communications on a merged object

dev.off()

#VP for PE signaling pathways: 

pathways.show <- pe.cellchat@netP$pathways

#Violin-plots: PE. 
for (i in 1:length(pathways.show)) {
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- plotGeneExpression(pe.cellchat, signaling = pathways.show[i])
  ggsave(filename=paste0(pathways.show[i], "_LR_contribution_Violin_PE.pdf"), plot=gg, dpi = 500)
}

pathways.show <- c("ANNEXIN") 

sources.use= c('dFB', 'DSC_1', 'dVEC', 'dLEC', 'dSMC')
targets.use= c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1', 'dTcell', 'dDC', 'dGranulocyte')

weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets

#par(mfrow = c(1,2), xpd=TRUE)

for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", targets.use = targets.use, sources.use=sources.use, 
                      edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

pathways.show <- c("FN1") 

sources.use= c('dFB', 'DSC_1', 'dVEC', 'dLEC', 'dSMC')
targets.use= c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1', 'dTcell', 'dDC', 'dGranulocyte')

weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets

#par(mfrow = c(1,2), xpd=TRUE)

for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", targets.use = targets.use, sources.use=sources.use, 
                      edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

pathways.show <- c("COLLAGEN") 

sources.use= c('dFB', 'DSC_1', 'dVEC', 'dLEC', 'dSMC')
targets.use= c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1', 'dTcell', 'dDC', 'dGranulocyte')

weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets

#par(mfrow = c(1,2), xpd=TRUE)

for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", targets.use = targets.use, sources.use=sources.use, 
                      edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

pathways.show <- c("LAMININ") 

sources.use= c('dFB', 'DSC_1', 'dVEC', 'dLEC', 'dSMC')
targets.use= c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1', 'dTcell', 'dDC', 'dGranulocyte')

weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets

#par(mfrow = c(1,2), xpd=TRUE)

for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", targets.use = targets.use, sources.use=sources.use, 
                      edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

pathways.show <- c("MHC-I") 

sources.use= c('dFB', 'DSC_1', 'dVEC', 'dLEC', 'dSMC')
targets.use= c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1', 'dTcell', 'dDC', 'dGranulocyte')

weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets

#par(mfrow = c(1,2), xpd=TRUE)

for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", targets.use = targets.use, sources.use=sources.use, 
                      edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

pathways.show <- c("APP") 

sources.use= c('dFB', 'DSC_1', 'dVEC', 'dLEC', 'dSMC')
targets.use= c('dMAC_activated', 'dMAC_classical', 'dNK_1', 'dNK_1', 'dTcell', 'dDC', 'dGranulocyte')

weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets

#par(mfrow = c(1,2), xpd=TRUE)

for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", targets.use = targets.use, sources.use=sources.use, 
                      edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

netVisual_chord_gene(pe.cellchat, sources.use = sources.use, targets.use = targets.use, lab.cex = 0.5, 
                     signaling = c("FN1", "LAMININ", "COLLAGEN"),legend.pos.x = 8,
                    legend.pos.y = 5)

netVisual_chord_gene(lateC.cellchat, sources.use = sources.use, targets.use = targets.use, lab.cex = 0.5, 
                     signaling = c("FN1", "LAMININ", "COLLAGEN"),legend.pos.x = 8,
                    legend.pos.y = 5)

netVisual_chord_gene(lateC.cellchat, sources.use = sources.use, targets.use = targets.use, lab.cex = 0.5,
                     signaling = c("APP", "MHC-I", "ANNEXIN"),legend.pos.x = 8,
                    legend.pos.y = 5)

netVisual_chord_gene(pe.cellchat, sources.use = sources.use, targets.use = targets.use, lab.cex = 0.5,
                     signaling = c("APP", "MHC-I", "ANNEXIN"),legend.pos.x = 8,
                    legend.pos.y = 5)

netVisual_chord_gene(lateC.cellchat, sources.use = sources.use, targets.use = targets.use, lab.cex = 0.4,
                     legend.pos.x = 15, legend.pos.y = 5, annotationTrackHeight = c(0.02),
                     big.gap = 8)

netVisual_chord_gene(pe.cellchat, sources.use = sources.use, targets.use = targets.use, lab.cex = 0.4,
                     legend.pos.x = 15, legend.pos.y = 5, annotationTrackHeight = c(0.02),
                     big.gap = 8)

pathways.show <- c("APP") 

pairLR.APP <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.APP[1,] # show one ligand-receptor pair



#Save the merged CellChat object:
saveRDS(cellchat, file = "/data/analysis/preeclampsia_2019/analysis/images/cellphonedb_analysis/decidua_immune_analysis/comparative_decidua_immune_analysis.rds")
