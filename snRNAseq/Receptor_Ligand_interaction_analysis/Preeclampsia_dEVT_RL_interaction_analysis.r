Sys.setenv(LANG="en_US.UTF-8")
Sys.setenv(LC_ALL="en_US.UTF-8")

suppressPackageStartupMessages(library(Seurat))

#Read the Seurat object: 
data <- readRDS(file = "/data/analysis/preeclampsia_2019/analysis/images/PE_markers_sexcorrected/updated_placenta_0704.rds")
data
table(Idents(data))

#Subset necessary cell-types from decidua: dEVT. 

subset = c("dEVT", "dEpC", "DSC_1", "DSC_2", "dSMC", "dFB_1", "dFB_2", "dVEC", "dLEC", "dLEC_dysfunctional",
          "dMonocyte", "dMAC_activated", "dMAC_classical", "dNK_1", "dNK_2", "dNK_prol", "dDC", "dGranulocyte",
          "dTcell", "dPlasmaCell") 

seurat_chat = subset(data, idents = subset) 
#seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["group"]])
seurat_chat
table(Idents(seurat_chat))

#Subset by Early, Late controls & PE samples: "group" metadata. 
Idents(seurat_chat) <- "group"

table(Idents(seurat_chat))

#Subset "late PE" only: decidua (diseased group only)
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

head(meta)

suppressPackageStartupMessages(library(CellChat))

#Create a CellChat object using data matrix as input: pe.input (from Seurat)
pe.cellchat <- createCellChat(object = pe.input, meta = meta, group.by = "group")

pe.cellchat #PE Cell chat object. 
levels(pe.cellchat@idents) #show factor levels of the cell labels

groupSize <- as.numeric(table(pe.cellchat@idents)) # number of cells in each cell group

CellChatDB <- CellChatDB.human #use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

#Take the whole database together:
#use all CellChatDB for cell-cell communication analysis
#CellChatDB.use <- CellChatDB # simply use the default CellChatDB
CellChatDB.use <- CellChatDB

pe.cellchat@DB <- CellChatDB.use

cellchat <- subsetData(pe.cellchat) #subset the expression data of signaling genes for saving computation cost

#Identify over-expressed signaling genes associated with each cell group: 
cellchat <- identifyOverExpressedGenes(cellchat)

#Identify over-expressed ligand-receptor interactions (pairs) within the used CellChatDB: 
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

head(cellchat@LR$LRsig)

#Compute the communication probability and infer cellular communication network:
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, type= "triMean")

#Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

#Extract the inferred cellular communication network as a data frame:
#We provide a function subsetCommunication to easily access the inferred cell-cell communications of interest. For example,
df.net <- subsetCommunication(cellchat)

head(df.net)

write.csv(df.net, file= "./dEVT_decidua/dEVT_cell_communication_triMean.csv")

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
                 title.name = "Number of interactions", sources.use= "dEVT", remove.isolate=TRUE,
                vertex.label.cex = 0.6, edge.width.max = 6)

groupSize <- as.numeric(table(cellchat@idents))

#dEVT:: targets. 
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Number of interactions", targets.use= "dEVT", remove.isolate=TRUE,
                vertex.label.cex = 0.5, edge.width.max = 6)

netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Interaction weights/strength", sources.use= "dEVT", remove.isolate=TRUE,
                vertex.label.cex = 0.6, edge.width.max = 6)

netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Interaction weights/strength", targets.use= "dEVT", remove.isolate=TRUE,
                vertex.label.cex = 0.5, edge.width.max = 6)

#Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways

pathways.show.all

cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)

df.net <- subsetCommunication(cellchat, sources.use= "dEVT")

head(df.net)

write.csv(df.net, file= "./dEVT_decidua/dEVT_source_outgoing_triMean.csv")

pathways.show <- c("CALCR") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
#vertex.receiver = seq(1,5) # a numeric vector. 

netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord", 
                    sources.use= "dEVT")

library(shiny)
library(cowplot)
library(ggplot2)
library(scales)
library(tidyverse)

pathways.show <- c("CALCR", "ANNEXIN", "COLLAGEN", "FN1", "LAMININ", "APP", "CD46", "PTPRM", "CDH", "SEMA6")

for (i in 1:length(pathways.show)) {
  #Visualize communication network associated with both signaling pathway and individual L-R pairs
  #netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  #P= netVisual_aggregate(cellchat, signaling = pathways.show[i], layout = "chord", sources.use= "dEVT")
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show[i])
  ggsave(filename=paste0(pathways.show[i], "_LR_contribution.pdf"), plot=gg, width = 5, height = 3, units = 'in', dpi = 300)
}

pathways.show <- c("CALCR", "ANNEXIN", "COLLAGEN", "FN1", "LAMININ", "APP", "CD46", "PTPRM", "CDH", "SEMA6")

#Violin plots: 
for (i in 1:length(pathways.show)) {
  #Visualize communication network associated with both signaling pathway and individual L-R pairs
  #netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- plotGeneExpression(cellchat, signaling = pathways.show[i])
  ggsave(filename=paste0(pathways.show[i], "_LR_contribution_Violin.pdf"), plot=gg, dpi = 500)
}



pathways.show <- c("CALCR", "ANNEXIN", "COLLAGEN", "FN1", "LAMININ", "APP", "CD46", "PTPRM", "CDH", "SEMA6")

#Dot-plots: 
for (i in 1:length(pathways.show)) {
  #Visualize communication network associated with both signaling pathway and individual L-R pairs
  #netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- plotGeneExpression(cellchat, signaling = pathways.show[i], type= "dot", cols = c("lightgrey", "blue"), assay= "RNA")
  ggsave(filename=paste0(pathways.show[i], "_LR_contribution_Violin.pdf"), plot=gg, dpi = 500)
}

pathways.show <- c("COLLAGEN") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
#vertex.receiver = seq(1,5) # a numeric vector. 

netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord", 
                    sources.use= "dEVT")

df.net <- subsetCommunication(cellchat, targets.use= "dEVT")

head(df.net)

write.csv(df.net, file= "./dEVT_decidua/dEVT_targets_incoming_triMean.csv")

pathways.show <- c("EDN", "NRXN", "PARs", "IGF")

#Dot-plots: 
for (i in 1:length(pathways.show)) {
  #Visualize communication network associated with both signaling pathway and individual L-R pairs
  #netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- plotGeneExpression(cellchat, signaling = pathways.show[i])
  ggsave(filename=paste0(pathways.show[i], "_LR_incoming_Violin.pdf"), plot=gg, dpi = 500)
  P <- netAnalysis_contribution(cellchat, signaling = pathways.show[i])
  ggsave(filename=paste0(pathways.show[i], "_LR_contribution.pdf"), plot=P, width = 5, height = 3, units = 'in', dpi = 300)
}



pdf("./dEVT_decidua/dEVT_source_bubblePlot_v1.pdf", paper= "special")

netVisual_bubble(cellchat, sources.use= "dEVT", remove.isolate = FALSE)

dev.off()

pdf("./dEVT_decidua/dEVT_targets_bubblePlot_v1.pdf", paper= "special")

netVisual_bubble(cellchat, targets.use= "dEVT", remove.isolate = FALSE)

dev.off()

saveRDS(cellchat, file = "/data/analysis/preeclampsia_2019/analysis/images/cellphonedb_analysis/dEVT_decidua/dEVT_decidua_cellchat_PE.rds")



Sys.setenv(LANG="en_US.UTF-8")
Sys.setenv(LC_ALL="en_US.UTF-8")

suppressPackageStartupMessages(library(CellChat))

data <- readRDS(file= "/data/analysis/preeclampsia_2019/analysis/images/cellphonedb_analysis/dEVT_decidua/dEVT_decidua_cellchat_PE.rds")

data

data <- netAnalysis_computeCentrality(data, slot.name = "netP")

ht1 <- netAnalysis_signalingRole_heatmap(data, pattern = "outgoing", width = 6, height = 12, color.heatmap = "GnBu")
ht2 <- netAnalysis_signalingRole_heatmap(data, pattern = "incoming", width = 6, height = 12, color.heatmap = "GnBu")
ht1 + ht2

#Overall signaling heatmap used in the manuscript: 
pdf("./dEVT_decidua/dEVT_signaling_overall.pdf", paper= "special")
netAnalysis_signalingRole_heatmap(data, pattern = "all", width = 6, height = 12, color.heatmap = "OrRd")
dev.off() 

pdf("./dEVT_decidua/dEVT_signaling_incoming_outgoing.pdf", paper= "special")

ht1 <- netAnalysis_signalingRole_heatmap(data, pattern = "outgoing", width = 6, height = 12, color.heatmap = "GnBu")
ht2 <- netAnalysis_signalingRole_heatmap(data, pattern = "incoming", width = 6, height = 12, color.heatmap = "GnBu")
ht1 + ht2

dev.off() 

levels(data@idents)

groupSize <- as.numeric(table(data@idents)) # number of cells in each cell group

cols.use= c("#F87060", "#BFFF80", "#FFF44F", "#CA6680", "#713E5A", "#BF3100", "#D76A03", "#EC9F05", "#8EA604", "#998650", "#E01A4F", "#7A9CC6",
            "#BDE4A7", "#009900", "#7D8CC4", "#EAC5D8", "#56CBF9", "#7C0B2B", "#0B1D51", "#646881")

pathways.show <- c("CALCR", "ANNEXIN", "COLLAGEN", "FN1", "LAMININ", "APP", "CD46", "PTPRM", "CDH", "SEMA6",
                  "EDN", "NRXN", "PARs", "IGF")

for (i in 1:length(pathways.show)) {
  gg <- plotGeneExpression(data, signaling = pathways.show[i], color.use=cols.use)
  ggsave(filename=paste0(pathways.show[i], "_LR_incoming_Violin.pdf"), plot=gg, dpi = 500)
}

plotGeneExpression(data, signaling= "CALCR", color.use=cols.use)


