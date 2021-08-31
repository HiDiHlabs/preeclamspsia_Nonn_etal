Sys.setenv(LANG="en_US.UTF-8")
Sys.setenv(LC_ALL="en_US.UTF-8")

suppressPackageStartupMessages(library(Seurat))

#Load the Seurat object: 
data <- readRDS(file = "/data/analysis/preeclampsia_2019/analysis/images/PE_markers_sexcorrected/updated_placenta_0704.rds")

data

table(Idents(data))

#Subset necessary cell-types from villi:

subset = c("vTcell", "vHBC", "vDC", "vVEC", "vMC", "vFB", "vVCT") 
seurat_chat = subset(data, idents = subset) 
#seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["group"]])
seurat_chat

table(Idents(seurat_chat))

#Subset by Early, Late controls & PE samples: "group" metadata. 
Idents(seurat_chat) <- "group"

table(Idents(seurat_chat))



#Subset "PE" only: villi. 
subset = c("Late_Villi_PE")

seurat_pe = subset(seurat_chat, idents = subset) 

#seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["group"]])

table(Idents(seurat_pe))

#Change the "idents" back to cell type:
Idents(seurat_pe) <- "cell_type_semifinal_v2"

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

#Subset "late controls" only: villi. 
subset = c("Late_Villi_C")

seurat_late = subset(seurat_chat, idents = subset) 
#seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["group"]])
table(Idents(seurat_late))

#Change the "idents" back to cell type:
Idents(seurat_late) <- "cell_type_semifinal_v2"
table(Idents(seurat_late))

#Get the assay/normalized data matrix
lateC.input <- GetAssayData(seurat_late, assay = "RNA", slot = "data") 
labels <- Idents(seurat_late)

#create a dataframe of the cell labels
meta <- data.frame(group = labels, row.names = names(labels)) 

suppressPackageStartupMessages(library(CellChat))

#Create a CellChat object using data matrix as input: (from Seurat)
lateC.cellchat <- createCellChat(object = lateC.input, meta = meta, group.by = "group")

lateC.cellchat #LateC Cell chat object.

groupSize <- as.numeric(table(lateC.cellchat@idents)) # number of cells in each cell group

groupSize

CellChatDB <- CellChatDB.human 

#CellChatDB.use <- CellChatDB # simply use the default CellChatDB
CellChatDB.use <- CellChatDB

lateC.cellchat@DB <- CellChatDB.use



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

df.net <- subsetCommunication(lateC.cellchat)

head(df.net)

write.csv(df.net, file= "./villi_immune_cellchat/villi_all_source_triMean.csv")

#Extract the inferred cellular communication network as a data frame:
#We provide a function subsetCommunication to easily access the inferred cell-cell communications of interest. For example,
df.net <- subsetCommunication(lateC.cellchat, sources.use= c('vTcell', 'vHBC', 'vDC'))

head(df.net)

write.csv(df.net, file= "./villi_immune_cellchat/villi_immune_source_triMean.csv")

#Extract the inferred cellular communication network as a data frame:
#We provide a function subsetCommunication to easily access the inferred cell-cell communications of interest. For example,
df.net <- subsetCommunication(lateC.cellchat, targets.use= c('vTcell', 'vHBC', 'vDC'))

head(df.net)

write.csv(df.net, file= "./villi_immune_cellchat/villi_immune_target_triMean.csv")

#Infer the cell-cell communication at a signaling pathway level:
lateC.cellchat <- computeCommunProbPathway(lateC.cellchat)

#Calculate the aggregated cell-cell communication network:
#We can calculate the aggregated cell-cell communication network by counting 
#the number of links or summarizing the communication probability:
lateC.cellchat <- aggregateNet(lateC.cellchat)

#Access all the signaling pathways showing significant communications
pathways.show.all <- lateC.cellchat@netP$pathways

pathways.show.all

pathways.show <- lateC.cellchat@netP$pathways

#Violin-plots: PE. 
for (i in 1:length(pathways.show)) {
  #Visualize communication network associated with both signaling pathway and individual L-R pairs
  #netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway with Violins 
  gg <- plotGeneExpression(lateC.cellchat, signaling = pathways.show[i])
  ggsave(filename=paste0(pathways.show[i], "_LR_contribution_Violin_lateC.pdf"), plot=gg, dpi = 500)
}

groupSize <- as.numeric(table(lateC.cellchat@idents))

#Visualize the interactions with circular layouts 
#Color code represents cell type. 
netVisual_circle(lateC.cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Number of interactions", sources.use= c('vTcell', 'vHBC', 'vDC'),
                vertex.label.cex = 0.6)

groupSize <- as.numeric(table(lateC.cellchat@idents))

netVisual_circle(lateC.cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Number of interactions", targets.use= c('vTcell', 'vHBC', 'vDC'),
                vertex.label.cex = 0.6)

pathways.show <- pathways.show.all

lateC.cellchat <- netAnalysis_computeCentrality(lateC.cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
#netAnalysis_signalingRole_network(lateC.cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)


#Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(lateC.cellchat, pattern = "outgoing", width = 6, height = 12, color.heatmap = "GnBu")
ht2 <- netAnalysis_signalingRole_heatmap(lateC.cellchat, pattern = "incoming", width = 6, height = 12, color.heatmap = "GnBu")
ht1 + ht2



#Subset "PE" only: villi. (disease group)
subset = c("Late_Villi_PE")

seurat_pe = subset(seurat_chat, idents = subset) 

#seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["group"]])

table(Idents(seurat_pe))

#Change the "idents" back to cell type:
Idents(seurat_pe) <- "cell_type_semifinal_v2"

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

write.csv(df.net.pe, file= "./villi_immune_cellchat/PE_villi_immune_all_triMean.csv")

df.net.pe <- subsetCommunication(pe.cellchat, targets.use= c('vTcell', 'vHBC', 'vDC'))

head(df.net.pe)

write.csv(df.net.pe, file= "./villi_immune_cellchat/PE_villi_immune_targets_triMean.csv")

df.net.pe <- subsetCommunication(pe.cellchat, sources.use= c('vTcell', 'vHBC', 'vDC'))

head(df.net.pe)

write.csv(df.net.pe, file= "./villi_immune_cellchat/PE_villi_immune_sources_triMean.csv")

#Infer the cell-cell communication at a signaling pathway level:
pe.cellchat <- computeCommunProbPathway(pe.cellchat)

#Calculate the aggregated cell-cell communication network:
#We can calculate the aggregated cell-cell communication network by counting 
#the number of links or summarizing the communication probability:
pe.cellchat <- aggregateNet(pe.cellchat)

#Access all the signaling pathways showing significant communications
pathways.show.pe <- pe.cellchat@netP$pathways

pathways.show.pe

groupSize <- as.numeric(table(lateC.cellchat@idents))

#PE: 
netVisual_circle(pe.cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Number of interactions", targets.use= c('vTcell', 'vHBC', 'vDC'),
                vertex.label.cex = 0.6)

groupSize <- as.numeric(table(lateC.cellchat@idents))

#PE: sources
netVisual_circle(pe.cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Number of interactions", sources.use= c('vTcell', 'vHBC', 'vDC'),
                vertex.label.cex = 0.6)

pathways.show <- pathways.show.pe

#Violin-plots: PE. 
for (i in 1:length(pathways.show)) {
  #Visualize communication network associated with both signaling pathway and individual L-R pairs
  #netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- plotGeneExpression(pe.cellchat, signaling = pathways.show[i])
  ggsave(filename=paste0(pathways.show[i], "_LR_contribution_Violin.pdf"), plot=gg, dpi = 500)
}

#Compute the network centrality scores: all pathways. 
pathways.show <- pathways.show.pe

pe.cellchat <- netAnalysis_computeCentrality(pe.cellchat, slot.name = "netP") 

#Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(pe.cellchat, pattern = "outgoing", width = 6, height = 12, color.heatmap = "GnBu")
ht2 <- netAnalysis_signalingRole_heatmap(pe.cellchat, pattern = "incoming", width = 6, height = 12, color.heatmap = "GnBu")
ht1 + ht2

library(patchwork)

#Merge LateC & PE lists: 
object.list <- list(LC = lateC.cellchat, PE = pe.cellchat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

cellchat #combined

#Compare the #interactions: 
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

#villi immune: targets
netVisual_diffInteraction(cellchat, weight.scale = T, vertex.weight=10, vertex.label.cex= 0.6,
                         edge.width.max= 5, alpha.edge = 0.8, vertex.size.max = 8)

#villi immune: source. 
netVisual_diffInteraction(cellchat, weight.scale = T, vertex.weight=10, vertex.label.cex= 0.6,
                         edge.width.max= 5, alpha.edge = 0.8, vertex.size.max = 8, 
                         sources.use= c('2', '4', '6'))

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]),
                 targets.use= c('vTcell', 'vHBC', 'vDC'))
}

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]),
                 sources.use= c('vTcell', 'vHBC', 'vDC'))
}

gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

library(ComplexHeatmap)

i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 6, height = 10)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 6, height = 10)
draw(ht1 + ht2, ht_gap = unit(1, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 6, height = 10, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 6, height = 10, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(1, "cm"))

#Overall signaling heatmap: aggregating outgoing and incoming signaling together. 
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 10, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 10, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

#Identify the upgulated and down-regulated signaling ligand-receptor pairs. 
pdf("./villi_immune_cellchat/villi_immune_source_LR.pdf", paper= "special")

netVisual_bubble(cellchat, sources.use= c('vTcell', 'vHBC', 'vDC'), comparison = c(1, 2), angle.x = 45)
#> Comparing communications on a merged object

dev.off()

#Identify the upgulated and down-regulated signaling ligand-receptor pairs. 
pdf("./villi_immune_cellchat/villi_immune_targets_LR.pdf", paper= "special", w=10, h=12)

netVisual_bubble(cellchat, targets.use= c('vTcell', 'vHBC', 'vDC'), comparison = c(1, 2), angle.x = 45)
#> Comparing communications on a merged object

dev.off()

#define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "PE"
#define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)

#extract the ligand-receptor pairs with upregulated ligands in PE
net.up <- subsetCommunication(cellchat, net = net, datasets = "PE",ligand.logFC = 0.2, receptor.logFC = NULL)

#extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in PE.
net.down <- subsetCommunication(cellchat, net = net, datasets = "LC",ligand.logFC = -0.1, receptor.logFC = -0.1)

#head(net)

#net_up <- subsetCommunication(cellchat, net = net, datasets = "PE", sources.use = c('vTcell', 'vHBC', 'vDC'))

#head(net_up)

#head(net.up)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

gene.up

gene.down 

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c('vTcell', 'vHBC', 'vDC'), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, targets.use = c('vTcell', 'vHBC', 'vDC'), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1

pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c('vTcell', 'vHBC', 'vDC'), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg2

pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, targets.use = c('vTcell', 'vHBC', 'vDC'), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg2

#extract the ligand-receptor pairs with upregulated targetss in PE
net.up.targets <- subsetCommunication(cellchat, net = net, datasets = "PE",ligand.logFC = NULL, receptor.logFC = 0.1, targets.use = c('vTcell', 'vHBC', 'vDC'))

head(net.up.targets)

write.csv(net.up.targets, file= "./villi_immune_cellchat/PE_immune_targets_upregulated_LRpairs02.csv")

#extract the ligand-receptor pairs with upregulated targetss in PE
net.up.targets <- subsetCommunication(cellchat, net = net, datasets = "PE",ligand.logFC = 0.2, receptor.logFC = NULL, targets.use = c('vTcell', 'vHBC', 'vDC'))

head(net.up.targets)

write.csv(net.up.targets, file= "./villi_immune_cellchat/PE_immune_targets_upregulated_LRpairs.csv")

#extract the ligand-receptor pairs with upregulated ligands in PE
net.up.source <- subsetCommunication(cellchat, net = net, datasets = "PE",ligand.logFC = 0.2, receptor.logFC = NULL, sources.use = c('vTcell', 'vHBC', 'vDC'))

head(net.up.source)

write.csv(net.up.source, file= "./villi_immune_cellchat/PE_immune_source_upregulated_LRpairs.csv")

#Save the merged CellChat object:
saveRDS(cellchat, file = "/data/analysis/preeclampsia_2019/analysis/images/cellphonedb_analysis/villi_immune_cellchat/Villi_immune_LateC_vs_PE.rds")



pathways.show <- c("THBS") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", sources.use = c('vTcell', 'vHBC', 'vDC'),
                      edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

pathways.show <- c("THBS") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", targets.use = c('vTcell', 'vHBC', 'vDC'),
                      edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

pathways.show <- c("APP") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", targets.use = c('vTcell', 'vHBC', 'vDC'),
                      edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

pathways.show <- c("APP") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
#vertex.receiver = seq(1,5) # a numeric vector. 

netVisual_aggregate(lateC.cellchat, signaling = pathways.show, layout = "chord", 
                    targets.use = c('vTcell', 'vHBC', 'vDC'))

pathways.show <- c("APP") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
#vertex.receiver = seq(1,5) # a numeric vector. 

netVisual_aggregate(pe.cellchat, signaling = pathways.show, layout = "chord", 
                    targets.use = c('vTcell', 'vHBC', 'vDC'))

pathways.show <- c("LAMININ") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", targets.use = c('vTcell', 'vHBC', 'vDC'),
                      edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}


