Sys.setenv(LANG="en_US.UTF-8")
Sys.setenv(LC_ALL="en_US.UTF-8")

suppressPackageStartupMessages(library(Seurat))

#Read the Seurat object: 
data <- readRDS(file = "/data/analysis/preeclampsia_2019/analysis/images/PE_markers_sexcorrected/updated_placenta_0704.rds")

data

#Subset necessary cell-types from decidua & villi: dLECP is removed since it's not affected by PE in terms of DEG. 
#Villi: SCT clusters & decidua: dSMC, dVEC & dLEC. 

subset = c("vSCT_1", "vSCT_2", "vtropho_15", "dSMC", "dVEC", "dLEC") 

seurat_chat = subset(data, idents = subset) 

#seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["group"]])

seurat_chat

table(Idents(seurat_chat))

#Subset by Early, Late controls & PE samples: "group" metadata. 
Idents(seurat_chat) <- "group"

table(Idents(seurat_chat))

#Subset diseased group i.e., "PE" only: for decidua & villi. 
subset = c("Late_Decidua_PE", "Late_Villi_PE")

seurat_pe = subset(seurat_chat, idents = subset) 

#seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["group"]])

table(Idents(seurat_pe))

#Change the "idents" back to cell type:
Idents(seurat_pe) <- "cell_type_semifinal_v2"

table(Idents(seurat_pe))

#Subset normal late gestation group i.e., "late controls" only: decidua & villi. 
subset = c("Late_Decidua_C", "Late_Villi_C")

seurat_late = subset(seurat_chat, idents = subset) 

#seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["group"]])

table(Idents(seurat_late))

#Change the "idents" back to cell type:
Idents(seurat_late) <- "cell_type_semifinal_v2"

table(Idents(seurat_late))

#Get normalized data matrix from the seurat objects: lateC dataset
lateC.input <- GetAssayData(seurat_late, assay = "RNA", slot = "data") 

labels <- Idents(seurat_late)

#create a dataframe of the cell labels
meta <- data.frame(group = labels, row.names = names(labels)) 

#Load CellChat 
suppressPackageStartupMessages(library(CellChat))

#Create a CellChat object using data matrix as input: (from Seurat)
lateC.cellchat <- createCellChat(object = lateC.input, meta = meta, group.by = "group")

lateC.cellchat #LateC Cell chat object.

groupSize <- as.numeric(table(lateC.cellchat@idents)) #number of cells in each cell group

groupSize

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
#showDatabaseCategory(CellChatDB)

#use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling

#set the used database in the object
lateC.cellchat@DB <- CellChatDB.use

lateC.cellchat <- subsetData(lateC.cellchat) #subset the expression data of signaling genes for saving computation cost

#Identify over-expressed signaling genes associated with each cell group: 
lateC.cellchat <- identifyOverExpressedGenes(lateC.cellchat)

#Identify over-expressed ligand-receptor interactions (pairs) within the used CellChatDB: 
lateC.cellchat <- identifyOverExpressedInteractions(lateC.cellchat)
lateC.cellchat <- projectData(lateC.cellchat, PPI.human)

head(lateC.cellchat@LR$LRsig)

#Compute the communication probability and infer cellular communication network:
lateC.cellchat <- computeCommunProb(lateC.cellchat, raw.use = TRUE, type= "triMean")

#Filter out the cell-cell communication if there are only few number of cells in certain cell groups
lateC.cellchat <- filterCommunication(lateC.cellchat, min.cells = 10)

#Extract the inferred cellular communication network as a data frame:
#We provide a function subsetCommunication to easily access the inferred cell-cell communications of interest. For example,
df.net <- subsetCommunication(lateC.cellchat, sources.use= c('vSCT_1', 'vSCT_2', 'vtropho_15'), 
                             targets.use= c('dSMC', 'dVEC', 'dLEC'))

head(df.net)

#Infer the cell-cell communication at a signaling pathway level:
lateC.cellchat <- computeCommunProbPathway(lateC.cellchat)

#Calculate the aggregated cell-cell communication network:
#We can calculate the aggregated cell-cell communication network by counting 
#the number of links or summarizing the communication probability:
lateC.cellchat <- aggregateNet(lateC.cellchat)

#Access all the signaling pathways showing significant communications
pathways.show.all <- lateC.cellchat@netP$pathways

pathways.show.all

groupSize <- as.numeric(table(lateC.cellchat@idents))

netVisual_circle(lateC.cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Number of interactions")

#Compute the network centrality scores: 
pathways.show <- c("GDF", "VEGF", "CALCR")

lateC.cellchat <- netAnalysis_computeCentrality(lateC.cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(lateC.cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)


#Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(lateC.cellchat, pattern = "outgoing", width = 6, height = 12, color.heatmap = "GnBu")
ht2 <- netAnalysis_signalingRole_heatmap(lateC.cellchat, pattern = "incoming", width = 6, height = 12, color.heatmap = "GnBu")
ht1 + ht2

#Compute the communication probability and infer cellular communication network:
#Set population.size = TRUE: abundant cell populations tend to send collectively stronger signals than the rare cell populations.
lateC.cellchat.uncorr <- computeCommunProb(lateC.cellchat, raw.use = TRUE, type= "triMean", population.size = TRUE)

#Filter out the cell-cell communication if there are only few number of cells in certain cell groups
lateC.cellchat.uncorr <- filterCommunication(lateC.cellchat.uncorr, min.cells = 10)

#Extract the inferred cellular communication network as a data frame:
#We provide a function subsetCommunication to easily access the inferred cell-cell communications of interest. For example,
df.net <- subsetCommunication(lateC.cellchat.uncorr, sources.use= c('vSCT_1', 'vSCT_2', 'vtropho_15'), 
                             targets.use= c('dSMC', 'dVEC', 'dLEC'))

head(df.net)

write.csv(df.net, file= "./TrophoL_mesenchymal_analysis/trophoLigands_mesoReceptors_communication_triMean_SizeUncorrected.csv")

#Infer the cell-cell communication at a signaling pathway level:
lateC.cellchat <- computeCommunProbPathway(lateC.cellchat)

#Calculate the aggregated cell-cell communication network:
#We can calculate the aggregated cell-cell communication network by counting 
#the number of links or summarizing the communication probability:
lateC.cellchat <- aggregateNet(lateC.cellchat)

#Access all the signaling pathways showing significant communications
pathways.show.all <- lateC.cellchat@netP$pathways

pathways.show.all

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
df.net.pe <- subsetCommunication(pe.cellchat, sources.use= c('vSCT_1', 'vSCT_2', 'vtropho_15'), 
                             targets.use= c('dSMC', 'dVEC', 'dLEC'))

head(df.net.pe)

write.csv(df.net.pe, file= "./TrophoL_mesenchymal_analysis/PE_trophoLigands_mesoReceptors_triMean.csv")

#Infer the cell-cell communication at a signaling pathway level:
pe.cellchat <- computeCommunProbPathway(pe.cellchat)

#Calculate the aggregated cell-cell communication network:
#We can calculate the aggregated cell-cell communication network by counting 
#the number of links or summarizing the communication probability:
pe.cellchat <- aggregateNet(pe.cellchat)

#Access all the signaling pathways showing significant communications
pathways.show.pe <- pe.cellchat@netP$pathways

pathways.show.pe

groupSize <- as.numeric(table(pe.cellchat@idents))

netVisual_circle(pe.cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Number of interactions")

#Compute the network centrality scores: 
pathways.show <- c("GDF", "VEGF", "CALCR", "ANGPT", "LEP")

pe.cellchat <- netAnalysis_computeCentrality(pe.cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(pe.cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)


# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(pe.cellchat, pattern = "outgoing", width = 6, height = 12, color.heatmap = "GnBu")
ht2 <- netAnalysis_signalingRole_heatmap(pe.cellchat, pattern = "incoming", width = 6, height = 12, color.heatmap = "GnBu")
ht1 + ht2

library(patchwork)

#Merge LateC & PE lists: 
object.list <- list(LC = lateC.cellchat, PE = pe.cellchat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

cellchat

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

netVisual_diffInteraction(cellchat, weight.scale = T, vertex.weight=10, vertex.label.cex= 0.6,
                         edge.width.max= 5, alpha.edge = 0.8, vertex.size.max = 8)

#Remove isolated nodes:
netVisual_diffInteraction(cellchat, weight.scale = T, vertex.weight=10, vertex.label.cex= 0.6,
                         edge.width.max= 5, alpha.edge = 0.8, vertex.size.max = 8,
                         remove.isolate = TRUE)

#Differential interaction weight: 
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",  
                          vertex.weight=10, vertex.label.cex= 0.6,
                         edge.width.max= 5, alpha.edge = 0.8, vertex.size.max = 8,
                         top=1)

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))

for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, 
                   edge.weight.max = weight.max[2],  
                   vertex.weight=10, vertex.label.cex= 0.6,
                         edge.width.max= 5, alpha.edge = 0.8, vertex.size.max = 8,
                   title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

df.net.pe <- subsetCommunication(cellchat, sources.use= c('vSCT_1', 'vSCT_2', 'vtropho_15'), 
                             targets.use= c('dSMC', 'dVEC', 'dLEC'))

head(df.net.pe)

#write.csv(df.net.pe, file= "./TrophoL_mesenchymal_analysis/PE_trophoLigands_mesoReceptors_triMean.csv")

cellchat

#Identify signaling groups based on their functional similarity:
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2

rankSimilarity(cellchat, type = "functional")

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

g= gg1 + gg2

ggsave("./TrophoL_mesenchymal_analysis/Comparative_pathway_triMean01.pdf", g, dpi = 500, limitsize=FALSE) 

library(ComplexHeatmap)

i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 6, height = 7)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 6, height = 7)
draw(ht1 + ht2, ht_gap = unit(1, "cm"))

ht3 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 6, height = 7, color.heatmap = "GnBu")
ht4 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 6, height = 7, color.heatmap = "GnBu")
#draw(ht1 + ht2, ht_gap = unit(1, "cm"))
draw(ht3 + ht4, ht_gap = unit(1, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 6, height = 7, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 6, height = 7, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

netVisual_bubble(cellchat, sources.use = 'vSCT_1', comparison = c(1, 2), angle.x = 45)
#> Comparing communications on a merged object

netVisual_bubble(cellchat, sources.use = 'vtropho_15', comparison = c(1, 2), angle.x = 45)
#> Comparing communications on a merged object

#netVisual_bubble(cellchat, sources.use = 'vSCT_2', comparison = c(1, 2), angle.x = 45)
#> Comparing communications on a merged object

pathways.show <- c("CALCR") 

weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets

#par(mfrow = c(1,2), xpd=TRUE)

for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}

pathways.show <- c("GDF") 

weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets

#par(mfrow = c(1,2), xpd=TRUE)

for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}

pathways.show <- c("VEGF") 

weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets

#par(mfrow = c(1,2), xpd=TRUE)

for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}

plotGeneExpression(pe.cellchat, signaling = "CALCR")

pathways.show <- c("CALCR") 

weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets

#par(mfrow = c(1,2), xpd=TRUE)

for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list)[i]),
                     sources.use= c('vSCT_1', 'vSCT_2', 'vtropho_15'),
                     targets.use= c('dSMC', 'dVEC', 'dLEC'))
}

plotGeneExpression(lateC.cellchat, signaling = "GDF")

plotGeneExpression(pe.cellchat, signaling = "GDF")

pathways.show <- c("GDF") 

weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets

#par(mfrow = c(1,2), xpd=TRUE)

for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list)[i]),
                     sources.use= c('vSCT_1', 'vSCT_2', 'vtropho_15'),
                     targets.use= c('dSMC', 'dVEC', 'dLEC'))
}

plotGeneExpression(lateC.cellchat, signaling = "VEGF")

pathways.show <- c("VEGF") 

weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets

#par(mfrow = c(1,2), xpd=TRUE)

for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list)[i]),
                     sources.use= c('vSCT_1', 'vSCT_2', 'vtropho_15'),
                     targets.use= c('dSMC', 'dVEC', 'dLEC'))
}

plotGeneExpression(lateC.cellchat, signaling = "LEP")

plotGeneExpression(pe.cellchat, signaling = "LEP")

pathways.show= c("LEP")
netVisual_aggregate(pe.cellchat, signaling = pathways.show, layout = "circle", 
                   sources.use= c('vSCT_1', 'vSCT_2', 'vtropho_15'),
                     targets.use= c('dSMC', 'dVEC', 'dLEC'))

# Circle plot
pathways.show= c("ANGPT")
netVisual_aggregate(pe.cellchat, signaling = pathways.show, layout = "circle")

#Save the merged CellChat object:
saveRDS(cellchat, file = "/data/analysis/preeclampsia_2019/analysis/images/cellphonedb_analysis/TrophoL_mesenchymal_analysis/TriMean_cellchat_comparison_LateC_vs_PE.rds")

cellchat@net

mat <- cellchat@net$weight
#par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

cellchat_subset= subsetCommunication(cellchat, sources.use= c('vSCT_1', 'vSCT_2', 'vtropho_15'), 
                             targets.use= c('dSMC', 'dVEC', 'dLEC'))

cellchat_subset

g3= netVisual_bubble(cellchat, sources.use = 'vtropho_15', comparison = c(1, 2), angle.x = 45, targets.use= c('dSMC', 'dVEC', 'dLEC'))
#> Comparing communications on a merged object

ggsave("./TrophoL_mesenchymal_analysis/vSCTjuvenile_comparison_triMean_robust.pdf", g3, dpi = 400, limitsize=FALSE) 

g2= netVisual_bubble(cellchat, sources.use = 'vSCT_2', comparison = 2, angle.x = 45, targets.use= c('dSMC', 'dVEC', 'dLEC'))
#> Comparing communications on a merged object

ggsave("./TrophoL_mesenchymal_analysis/vSCT2_comparison_triMean_robust.pdf", g2, dpi = 400, limitsize=FALSE, height=4, width=4) 

g1= netVisual_bubble(cellchat, sources.use = 'vSCT_1', comparison = c(1, 2), angle.x = 45, targets.use= c('dSMC', 'dVEC', 'dLEC'))
#> Comparing communications on a merged object

ggsave("./TrophoL_mesenchymal_analysis/vSCT1_comparison_triMean_robust.pdf", g1, dpi = 400, limitsize=FALSE) 

#define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "PE"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "PE",ligand.logFC = 0.1, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "LC",ligand.logFC = -0.1, receptor.logFC = -0.1)

library(NMF)
library(ggalluvial)

selectK(pe.cellchat, pattern = "outgoing")

nPatterns = 2
pe.cellchat <- identifyCommunicationPatterns(pe.cellchat, pattern = "outgoing", k = nPatterns)

# river plot
netAnalysis_river(pe.cellchat, pattern = "outgoing")
#> Please make sure you have load `library(ggalluvial)` when running this function


#dot plot
netAnalysis_dot(pe.cellchat, pattern = "outgoing")


selectK(pe.cellchat, pattern = "incoming")

nPatterns = 3
pe.cellchat <- identifyCommunicationPatterns(pe.cellchat, pattern = "incoming", k = nPatterns)


#River plot
netAnalysis_river(pe.cellchat, pattern = "incoming")
#> Please make sure you have load `library(ggalluvial)` when running this function


#dot plot
netAnalysis_dot(pe.cellchat, pattern = "incoming")

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(pe.cellchat)

#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(pe.cellchat, signaling = c("GDF", "CALCR", "LEP"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg2


