Sys.setenv(LANG="en_US.UTF-8")
Sys.setenv(LC_ALL="en_US.UTF-8")

suppressPackageStartupMessages(library(Seurat))

#Read the Seurat object: 
data <- readRDS(file = "/data/analysis/preeclampsia_2019/analysis/images/PE_markers_sexcorrected/updated_placenta_0704.rds")

data

table(Idents(data))

#Subset necessary cell-types from decidua & DSCT (deported SCTs)

subset = c("dSCT", "dVEC", "dLEC", "dSMC", "dLEC_dysfunctional") 

seurat_chat = subset(data, idents = subset) 

#seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["group"]])

seurat_chat

table(Idents(seurat_chat))

#Subset by Early, Late controls & PE samples: "group" metadata. 
Idents(seurat_chat) <- "group"

table(Idents(seurat_chat))

#Subset disease group i.e., "late PE" only: in decidua 
subset = "Late_Decidua_PE" 
seurat_late = subset(seurat_chat, idents = subset) 
#seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["group"]])
seurat_late

#Change the "idents" back to cell type:
Idents(seurat_late) <- "cell_type_semifinal_v2"
table(Idents(seurat_late))

pe.input <- GetAssayData(seurat_late, assay = "RNA", slot = "data") #Normalized data matrix. 
labels <- Idents(seurat_late)
meta <- data.frame(group = labels, row.names = names(labels)) #create a dataframe of the cell labels. 

#pe.input

head(meta)

suppressPackageStartupMessages(library(CellChat))

#Create a CellChat object using data matrix as input: pe.input (from Seurat)
pe.cellchat <- createCellChat(object = pe.input, meta = meta, group.by = "group")

pe.cellchat #PE Cell chat object. 
levels(pe.cellchat@idents) # show factor levels of the cell labels

groupSize <- as.numeric(table(pe.cellchat@idents)) # number of cells in each cell group

groupSize #PE cells per group

CellChatDB <- CellChatDB.human #use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

#Take the whole database together:
#use all CellChatDB for cell-cell communication analysis
#CellChatDB.use <- CellChatDB # simply use the default CellChatDB
CellChatDB.use <- CellChatDB

pe.cellchat@DB <- CellChatDB.use

#pe.cellchat@DB

cellchat <- subsetData(pe.cellchat) #subset the expression data of signaling genes 

#Identify over-expressed signaling genes associated with each cell group: 
cellchat <- identifyOverExpressedGenes(cellchat)

#Identify over-expressed ligand-receptor interactions (pairs) within the used CellChatDB: 
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

head(cellchat@LR$LRsig)

write.csv(df.net, file= "./dSCT_cellchat_analysis/dSCT_cell_communication_raw.csv")

#Compute the communication probability and infer cellular communication network:
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, type= "triMean")

#Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

#Extract the inferred cellular communication network as a data frame:
#We provide a function subsetCommunication to easily access the inferred cell-cell communications of interest. For example,
df.net <- subsetCommunication(cellchat)

head(df.net)

write.csv(df.net, file= "./dSCT_cellchat_analysis/dSCT_cell_communication_triMean.csv")

#Infer the cell-cell communication at a signaling pathway level:
cellchat <- computeCommunProbPathway(cellchat)

#Calculate the aggregated cell-cell communication network:
#We can calculate the aggregated cell-cell communication network by counting 
#the number of links or summarizing the communication probability:
cellchat <- aggregateNet(cellchat)

#Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways

pathways.show.all

groupSize <- as.numeric(table(cellchat@idents))

netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Number of interactions")

netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
#par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}


pathways.show <- c("VEGF") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,5) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

netAnalysis_contribution(cellchat, signaling = pathways.show)

cellchat

pathways.show <- c("VEGF") 
# Hierarchy plot
#Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,2) #a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show, vertex.receiver=vertex.receiver)


pathways.show <- c("GDF") 
# Hierarchy plot
#Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,3) #a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show, vertex.receiver=vertex.receiver)


pathways.show <- c("VISFATIN") 
# Hierarchy plot
#Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,3) #a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show, vertex.receiver=vertex.receiver)


pathways.show.all

pathways.show <- c("TGFb")

#Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

netAnalysis_contribution(cellchat, signaling = pathways.show)

#show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')

netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1:5), remove.isolate = FALSE)

#> Comparing communications on a single object

netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1:5), remove.isolate = FALSE)

netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1:5), remove.isolate = FALSE)

netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:5), remove.isolate = FALSE)

netVisual_bubble(cellchat, sources.use = 5, targets.use = c(1:5), remove.isolate = FALSE)

for (i in 1:length(pathways.show.all)) {
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway with Violins
  gg <- plotGeneExpression(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_LR_contribution_Violin.pdf"), plot=gg, dpi = 500)
}

#Chord diagrams: 
for (i in 1:length(pathways.show.all)) {
  #Visualize communication network associated with both signaling pathway and individual L-R pairs
  #netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  netVisual_aggregate(cellchat, signaling = pathways.show.all[i], layout = "chord")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_LR_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}

#Plot the signaling gene expression distribution using violin/dot plot:
plotGeneExpression(cellchat, signaling = "TGFb")

plotGeneExpression(cellchat, signaling = "TGFb", enriched.only = FALSE)

#Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)


#Visualize the dominant senders (sources) and receivers (targets) in a 2D space:
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("GDF", "FN1"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2

#Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", width = 6, height = 12, color.heatmap = "GnBu")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", width = 6, height = 12, color.heatmap = "GnBu")
ht1 + ht2



cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)

library(NMF)
library(ggalluvial)

selectK(cellchat, pattern = "outgoing")


nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)

#dot-plot: 
netAnalysis_dot(cellchat, pattern = "outgoing")

#river plot
netAnalysis_river(cellchat, pattern = "outgoing")
#> Please make sure you have load `library(ggalluvial)` when running this function

selectK(cellchat, pattern = "incoming")

nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)

#river plot: 
netAnalysis_river(cellchat, pattern = "incoming")
#> Please make sure you have load `library(ggalluvial)` when running this function. 

#dot plot: 
netAnalysis_dot(cellchat, pattern = "incoming")

saveRDS(cellchat, file = "/data/analysis/preeclampsia_2019/analysis/images/cellphonedb_analysis/dSCT_cellchat_analysis/dSCT_cellchat_PE.rds")

pathways.show <- c("TENASCIN") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
#vertex.receiver = seq(1,5) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

netAnalysis_contribution(cellchat, signaling = pathways.show)


Sys.setenv(LANG="en_US.UTF-8")
Sys.setenv(LC_ALL="en_US.UTF-8")

suppressPackageStartupMessages(library(CellChat))

data <- readRDS(file = "/data/analysis/preeclampsia_2019/analysis/images/cellphonedb_analysis/dSCT_cellchat_analysis/dSCT_cellchat_PE.rds")

data

ht1 <- netAnalysis_signalingRole_heatmap(data, pattern = "outgoing", width = 6, height = 12, color.heatmap = "GnBu")
ht2 <- netAnalysis_signalingRole_heatmap(data, pattern = "incoming", width = 6, height = 12, color.heatmap = "GnBu")
ht1 + ht2

pdf("./dSCT_cellchat_analysis/dSCT_signaling_incoming_outgoing.pdf", paper= "special")
ht1 <- netAnalysis_signalingRole_heatmap(data, pattern = "outgoing", width = 6, height = 12, color.heatmap = "GnBu")
ht2 <- netAnalysis_signalingRole_heatmap(data, pattern = "incoming", width = 6, height = 12, color.heatmap = "GnBu")
ht1 + ht2
dev.off() 

#Used in the Extended manuscript figure: 
pdf("./dSCT_cellchat_analysis/dSCT_signaling_overall.pdf", paper= "special")
netAnalysis_signalingRole_heatmap(data, pattern = "all", width = 6, height = 12, color.heatmap = "OrRd")
dev.off() 

levels(data@idents)
groupSize <- as.numeric(table(data@idents)) # number of cells in each cell group

cols.use= c("#F87060", "#BFFF80", "#8EA604", "#FB4B4E", "#0B1D51")

#Access all the signaling pathways showing significant communications
pathways.show.all <- data@netP$pathways
pathways.show.all

for (i in 1:length(pathways.show.all)) {
  gg <- plotGeneExpression(data, signaling = pathways.show.all[i], color.use=cols.use)
  ggsave(filename=paste0(pathways.show.all[i], "_dSCT_Violin.pdf"), plot=gg, dpi = 500)
}



plotGeneExpression(data, signaling= "VEGF", color.use=cols.use)

sessionInfo()


