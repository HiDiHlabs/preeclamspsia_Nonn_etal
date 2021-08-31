Sys.setenv(LANG="en_US.UTF-8")
Sys.setenv(LC_ALL="en_US.UTF-8")

suppressPackageStartupMessages(library(Seurat))

#Read the Seurat object: 
data <- readRDS(file = "/data/analysis/preeclampsia_2019/analysis/images/PE_markers_sexcorrected/updated_placenta_0704.rds")

data

#Subset necessary cell-types from decidua dMSC. (plus EVT from the trophoblasts)

subset = c("dVEC", "dLEC", "dSMC", "dMSC", "dFB_1", "dFB_2", "DSC_1", "DSC_2", "vEVT") 

seurat_chat = subset(data, idents = subset) 
seurat_chat
table(Idents(seurat_chat))

#Subset by Early, Late controls & PE samples: "group" metadata. 
Idents(seurat_chat) <- "group"

table(Idents(seurat_chat))

#Subset Early_Decidua_C for dMSC matrisome remodelling. 
#We won't do for late controls & PE as EVT is severely depleted there. 
subset = c("Early_Decidua_C", "Early_Villi_C") 

seurat_early = subset(seurat_chat, idents = subset) 
table(Idents(seurat_early))

#Change the "idents" back to cell type:
Idents(seurat_early) <- "cell_type_semifinal_v2"
table(Idents(seurat_early))

#Merge dFB1 & dFB2 to dFB to avoid class imbalance 
seurat_early <- RenameIdents(object = seurat_early,  'dFB_1' = 'dFB', 'dFB_2' = 'dFB')
table(Idents(seurat_early))

#old_order= ['dFB', 'dSMC', 'dVEC', 'dLEC']
order= ['dVEC', 'dSMC', 'dLEC', 'dFB', 'DSC_1', 'DSC_2', 'dMSC']

#old_cols_list= ['#bf3100', '#0b1d51', '#f87060', '#bfff80']
cols_list= ['#f87060', '#0b1d51', '#bfff80', '#bf3100', '#56cbf9', '#7c0b2b', '#c0c999']

cols.evt= c("#BF3100", "#C0C999", "#F87060", "#bFFF80", "#56CBF9", "#7C0B2B", "#0B1d51", "#7D8CC4")

seurat_early.input <- GetAssayData(seurat_early, assay = "RNA", slot = "data") #Extract normalized data matrix

labels <- Idents(seurat_early)
meta <- data.frame(group = labels, row.names = names(labels)) #create a dataframe of the cell labels

suppressPackageStartupMessages(library(CellChat))

#Create a CellChat object using data matrix as input: (from Seurat)
early.cellchat <- createCellChat(object = seurat_early.input, meta = meta, group.by = "group")

early.cellchat

CellChatDB <- CellChatDB.human #use CellChatDB.mouse if running on mouse data
CellChatDB.use <- CellChatDB

#set the used database in the object
early.cellchat@DB <- CellChatDB.use

groupSize <- as.numeric(table(early.cellchat@idents)) # number of cells in each cell group
groupSize

early.cellchat <- subsetData(early.cellchat) #subset the expression data of signaling genes for saving computation cost

#Identify over-expressed signaling genes associated with each cell group: 
early.cellchat <- identifyOverExpressedGenes(early.cellchat)

#Identify over-expressed ligand-receptor interactions (pairs) within the used CellChatDB: 
early.cellchat <- identifyOverExpressedInteractions(early.cellchat)
early.cellchat <- projectData(early.cellchat, PPI.human)

#Compute the communication probability and infer cellular communication network:
early.cellchat <- computeCommunProb(early.cellchat, raw.use = TRUE, type= "triMean")

#Filter out the cell-cell communication if there are only few number of cells in certain cell groups
early.cellchat <- filterCommunication(early.cellchat, min.cells = 10)

#Extract the inferred cellular communication network as a data frame:
#We provide a function subsetCommunication to easily access the inferred cell-cell communications of interest. For example,
df.net <- subsetCommunication(early.cellchat)

head(df.net)

write.csv(df.net, file= "./EVT_matrisome_analysis/EVT_early_matrisome_all_trimean.csv") 

#vEVT: source is the ligands from the extravillous trophoblast that bind to the decidual mesoendothelial receptors to possibly remodel them. 
df.net2 <- subsetCommunication(early.cellchat, sources.use= "vEVT")

head(df.net2)

write.csv(df.net2, file= "./EVT_matrisome_analysis/EVTsource_early_matrisome_trimean.csv") 

df.net3 <- subsetCommunication(early.cellchat, targets.use= "vEVT")

head(df.net3)

write.csv(df.net3, file= "./EVT_matrisome_analysis/EVT_target_early_matrisome_trimean.csv") 

#Infer the cell-cell communication at a signaling pathway level:
early.cellchat <- computeCommunProbPathway(early.cellchat)

#Calculate the aggregated cell-cell communication network:
#We can calculate the aggregated cell-cell communication network by counting 
#the number of links or summarizing the communication probability:
early.cellchat <- aggregateNet(early.cellchat)

#Access all the signaling pathways showing significant communications
pathways.show.all <- early.cellchat@netP$pathways
pathways.show.all

plotGeneExpression(early.cellchat, signaling= "FGF", color.use= cols.evt)

plotGeneExpression(early.cellchat, signaling= "ANGPT", color.use= cols.evt)

pathways.show <- early.cellchat@netP$pathways

#Violin-plots: LateC. 
for (i in 1:length(pathways.show)) {
  gg <- plotGeneExpression(early.cellchat, signaling = pathways.show[i], color.use= cols.evt)
  ggsave(filename=paste0(pathways.show[i], "_LR_Violin_Early.pdf"), plot=gg, dpi = 600)
}

netVisual_bubble(early.cellchat, sources.use = 'vEVT', remove.isolate= TRUE)

pdf("./EVT_matrisome_analysis/vEVT_source_matrisome_bubble_v2.pdf", w=4, h=9, paper= "special")

netVisual_bubble(early.cellchat, sources.use = 'vEVT', remove.isolate= TRUE)

dev.off() 

netVisual_bubble(early.cellchat, targets.use = 'vEVT', remove.isolate= TRUE)

pdf("./EVT_matrisome_analysis/vEVT_targets_matrisome_bubble_v2.pdf", w=4, h=7, paper= "special")
netVisual_bubble(early.cellchat, targets.use = 'vEVT', remove.isolate= TRUE)

dev.off() 

#Save the merged CellChat object:
saveRDS(early.cellchat, file = "/data/analysis/preeclampsia_2019/analysis/images/cellphonedb_analysis/EVT_matrisome_analysis/vEVT_matrix_early_cellchat.rds")

#Compute centrality (early):
early.cellchat <- netAnalysis_computeCentrality(early.cellchat, slot.name = "netP")

early.cellchat

pdf("./EVT_matrisome_analysis/vEVT_matrisome_signaling_heatmap01.pdf", paper= "special")

#Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(early.cellchat, pattern = "outgoing", width = 6, height = 12, color.heatmap = "GnBu")
ht2 <- netAnalysis_signalingRole_heatmap(early.cellchat, pattern = "incoming", width = 6, height = 12, color.heatmap = "GnBu")
ht1 + ht2

dev.off()

pdf("./EVT_matrisome_analysis/vEVT_matrisome_signaling_heatmap02.pdf", paper= "special")

#Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways in early pregnancy. 
ht3 <- netAnalysis_signalingRole_heatmap(early.cellchat, pattern = "all", width = 6, height = 12, color.heatmap = "OrRd")
ht3

dev.off()

Sys.setenv(LANG="en_US.UTF-8")
Sys.setenv(LC_ALL="en_US.UTF-8")
suppressPackageStartupMessages(library(CellChat))

data <- readRDS(file = "/data/analysis/preeclampsia_2019/analysis/images/cellphonedb_analysis/EVT_matrisome_analysis/vEVT_matrix_early_cellchat.rds")

data

#data@netP

netVisual_bubble(data, sources.use = 'vEVT', targets.use= c('dFB', 'dVEC', 'dLEC', 'DSC_1', 'DSC_2', 'dMSC'), 
                 remove.isolate= TRUE)

pdf("./EVT_matrisome_analysis/vEVT_source_matrisome_bubble_v3.pdf", w=3.5, h=5.5, paper= "special")

netVisual_bubble(data, sources.use = 'vEVT', targets.use= c('dFB', 'dVEC', 'dLEC', 'DSC_1', 'DSC_2', 'dMSC'), 
                 remove.isolate= TRUE, thresh = 0.001, color.heatmap = "Spectral", font.size=8)

dev.off() 

pdf("./EVT_matrisome_analysis/vEVT_source_matrisome_bubble_v3_maxquantile.pdf", w=3.5, h=5.5, paper= "special")

#Bubble plot: Ligands (vEVT) modulating targets/receptors (mesoendothelial): used for manuscript. 
netVisual_bubble(data, sources.use = 'vEVT', targets.use= c('dFB', 'dVEC', 'dLEC', 'DSC_1', 'DSC_2', 'dMSC'), 
                 remove.isolate= TRUE, thresh = 0.001, color.heatmap = "Spectral", max.quantile=0.8)

dev.off() 

pairLR.use <- extractEnrichedLR(data, signaling = c("VEGF", "VISFATIN", "PTN", "FN1", "COLLAGEN", "CALCR", "LAMININ", 
                                                      "ADGRE5", "APP", "CDH5", "EPHA", "DESMOSOME", "NECTIN"))

pdf("./EVT_matrisome_analysis/vEVTsource_matrisome_bubble_v3_enriched_maxquantile.pdf", w=3.5, h=5.5, paper= "special")


netVisual_bubble(data, sources.use = 'vEVT', targets.use= c('dFB', 'dVEC', 'dLEC', 'DSC_1', 'DSC_2', 'dMSC'), 
                 remove.isolate= TRUE, thresh = 0.001, color.heatmap = "Spectral", max.quantile=0.8,
                pairLR.use = pairLR.use)

dev.off() 

netVisual_bubble(data, sources.use = 'vEVT', targets.use= c('dFB', 'dVEC', 'dLEC', 'DSC_1', 'DSC_2', 'dMSC'), 
                 remove.isolate= TRUE, thresh = 0.001, color.heatmap = "Spectral", max.quantile=0.8,
                signaling = c("FN1", "COLLAGEN", "LAMININ", "DESMOSOME", 
                                                       "VEGF", "VISFATIN", "PTN",  "CALCR", 
                                                      "ADGRE5", "APP", "CDH5", "EPHA", "NECTIN"))





#Compute centrality (early):
data <- netAnalysis_computeCentrality(data, slot.name = "netP")

data

#Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways

pdf("./EVT_matrisome_analysis/vEVTsource_pathways_heatmap_v2.pdf", paper= "special")

ht3 <- netAnalysis_signalingRole_heatmap(data, pattern = "all", 
                                         signaling= c("VEGF", "VISFATIN", "PTN", "FN1", "COLLAGEN", "CALCR", "LAMININ", 
                                                      "ADGRE5", "APP", "CDH5", "EPHA", "DESMOSOME", "NECTIN"), 
                                         width = 6, height = 10, cluster.rows=TRUE, 
                                         color.heatmap = "OrRd")
ht3

dev.off() 

pdf("./EVT_matrisome_analysis/vEVTsource_bubble_Vascular01.pdf", w=3, h=5.5, paper= "special")


netVisual_bubble(data, sources.use = 'vEVT', targets.use= c('dFB', 'dVEC', 'dLEC'), 
                 remove.isolate= TRUE, thresh = 0.001, color.heatmap = "Spectral", max.quantile=0.8,
                pairLR.use = pairLR.use, font.size=8)

dev.off() 

data

pdf("./EVT_matrisome_analysis/vEVTsource_bubble_Stromal01.pdf", w=3, h=5.5, paper= "special")


netVisual_bubble(data, sources.use = 'vEVT', targets.use= c('dMSC', 'DSC_1', 'DSC_2'), 
                 remove.isolate= TRUE, thresh = 0.001, color.heatmap = "Spectral", max.quantile=0.8,
                pairLR.use = pairLR.use, font.size=8)

dev.off() 

#Identify signaling groups based on their functional similarity:
data <- computeNetSimilarity(data, type = "functional")
data <- netEmbedding(data, type = "functional")

#Manifold learning of the signaling networks for a single dataset
data <- netClustering(data, type = "functional")

#Visualization in 2D-space
netVisual_embedding(data, type = "functional", label.size = 3.5)

#Identify signaling groups based on structure similarity:
data <- computeNetSimilarity(data, type = "structural")
data <- netEmbedding(data, type = "structural")

#Manifold learning of the signaling networks for a single dataset
data <- netClustering(data, type = "structural")

#Visualization in 2D-space
netVisual_embedding(data, type = "structural", label.size = 3.5)

#data@net

#Don't refer or run this: 

netVisual_bubble <- function(object, sources.use = NULL, targets.use = NULL, signaling = NULL, pairLR.use = NULL, color.heatmap = c("Spectral","viridis"), n.colors = 10, direction = -1, thresh = 0.05,
                             comparison = NULL, group = NULL, remove.isolate = FALSE, max.dataset = NULL, min.dataset = NULL,
                             min.quantile = 0, max.quantile = 1, line.on = TRUE, line.size = 0.2, color.text.use = TRUE, color.text = NULL,
                             title.name = NULL, font.size = 10, font.size.title = 10, show.legend = TRUE,
                             grid.on = TRUE, color.grid = "grey90", angle.x = 90, vjust.x = NULL, hjust.x = NULL,
                             return.data = FALSE){
  color.heatmap <- match.arg(color.heatmap)
  if (is.list(object@net[[1]])) {
    message("Comparing communications on a merged object \n")
  } else {
    message("Comparing communications on a single object \n")
  }
  if (is.null(vjust.x) | is.null(hjust.x)) {
    angle=c(0, 45, 90)
    hjust=c(0, 1, 1)
    vjust=c(0, 1, 0.5)
    vjust.x = vjust[angle == angle.x]
    hjust.x = hjust[angle == angle.x]
  }
  if (length(color.heatmap) == 1) {
    color.use <- tryCatch({
      RColorBrewer::brewer.pal(n = n.colors, name = color.heatmap)
    }, error = function(e) {
      scales::viridis_pal(option = color.heatmap, direction = -1)(n.colors)
    })
  } else {
    color.use <- color.heatmap
  }
  if (direction == -1) {
    color.use <- rev(color.use)
  }

  if (is.null(comparison)) {
    cells.level <- levels(object@idents)
    if (is.numeric(sources.use)) {
      sources.use <- cells.level[sources.use]
    }
    if (is.numeric(targets.use)) {
      targets.use <- cells.level[targets.use]
    }
    df.net <- subsetCommunication(object, slot.name = "net",
                                  sources.use = sources.use, targets.use = targets.use,
                                  signaling = signaling,
                                  pairLR.use = pairLR.use,
                                  thresh = thresh)
    df.net$source.target <- paste(df.net$source, df.net$target, sep = " -> ")
    source.target <- paste(rep(sources.use, each = length(targets.use)), targets.use, sep = " -> ")
    source.target.isolate <- setdiff(source.target, unique(df.net$source.target))
    if (length(source.target.isolate) > 0) {
      df.net.isolate <- as.data.frame(matrix(NA, nrow = length(source.target.isolate), ncol = ncol(df.net)))
      colnames(df.net.isolate) <- colnames(df.net)
      df.net.isolate$source.target <- source.target.isolate
      df.net.isolate$interaction_name_2 <- df.net$interaction_name_2[1]
      df.net.isolate$pval <- 1
      a <- stringr::str_split(df.net.isolate$source.target, " -> ", simplify = T)
      df.net.isolate$source <- as.character(a[, 1])
      df.net.isolate$target <- as.character(a[, 2])
      df.net <- rbind(df.net, df.net.isolate)
    }

    df.net$pval[df.net$pval > 0.05] = 1
    df.net$pval[df.net$pval > 0.01 & df.net$pval <= 0.05] = 2
    df.net$pval[df.net$pval <= 0.01] = 3
    df.net$prob[df.net$prob == 0] <- NA
    df.net$prob.original <- df.net$prob
    df.net$prob <- -1/log(df.net$prob)

    idx1 <- which(is.infinite(df.net$prob) | df.net$prob < 0)
    if (sum(idx1) > 0) {
      values.assign <- seq(max(df.net$prob, na.rm = T)*1.1, max(df.net$prob, na.rm = T)*1.5, length.out = length(idx1))
      position <- sort(prob.original[idx1], index.return = TRUE)$ix
      df.net$prob[idx1] <- values.assign[match(1:length(idx1), position)]
    }
    # rownames(df.net) <- df.net$interaction_name_2

    df.net$source <- factor(df.net$source, levels = cells.level[cells.level %in% unique(df.net$source)])
    df.net$target <- factor(df.net$target, levels = cells.level[cells.level %in% unique(df.net$target)])
    group.names <- paste(rep(levels(df.net$source), each = length(levels(df.net$target))), levels(df.net$target), sep = " -> ")

    df.net$interaction_name_2 <- as.character(df.net$interaction_name_2)
    #df.net <- with(df.net, df.net[order(interaction_name_2),])
    df.net <- with(df.net, df.net[order(interaction_name_2),])
    df.net$interaction_name_2 <- factor(df.net$interaction_name_2, levels = unique(df.net$interaction_name_2))
    #df.net$interaction_name_2 <- factor(df.net$interaction_name_2, levels = unique(df.net$interaction_name_2)
    cells.order <- group.names
    df.net$source.target <- factor(df.net$source.target, levels = cells.order) 
    df <- df.net
  } else {
    dataset.name <- names(object@net)
    df.net.all <- subsetCommunication(object, slot.name = "net",
                                      sources.use = sources.use, targets.use = targets.use,
                                      signaling = signaling,
                                      pairLR.use = pairLR.use,
                                      thresh = thresh)
    df.all <- data.frame()
    for (ii in 1:length(comparison)) {
      cells.level <- levels(object@idents[[comparison[ii]]])
      if (is.numeric(sources.use)) {
        sources.use <- cells.level[sources.use]
      }
      if (is.numeric(targets.use)) {
        targets.use <- cells.level[targets.use]
      }

      df.net <- df.net.all[[comparison[ii]]]
      df.net$interaction_name_2 <- as.character(df.net$interaction_name_2)
      df.net$source.target <- paste(df.net$source, df.net$target, sep = " -> ")
      source.target <- paste(rep(sources.use, each = length(targets.use)), targets.use, sep = " -> ")
      source.target.isolate <- setdiff(source.target, unique(df.net$source.target))
      if (length(source.target.isolate) > 0) {
        df.net.isolate <- as.data.frame(matrix(NA, nrow = length(source.target.isolate), ncol = ncol(df.net)))
        colnames(df.net.isolate) <- colnames(df.net)
        df.net.isolate$source.target <- source.target.isolate
        df.net.isolate$interaction_name_2 <- df.net$interaction_name_2[1]
        df.net.isolate$pval <- 1
        a <- stringr::str_split(df.net.isolate$source.target, " -> ", simplify = T)
        df.net.isolate$source <- as.character(a[, 1])
        df.net.isolate$target <- as.character(a[, 2])
        df.net <- rbind(df.net, df.net.isolate)
      }

      df.net$source <- factor(df.net$source, levels = cells.level[cells.level %in% unique(df.net$source)])
      df.net$target <- factor(df.net$target, levels = cells.level[cells.level %in% unique(df.net$target)])
      group.names <- paste(rep(levels(df.net$source), each = length(levels(df.net$target))), levels(df.net$target), sep = " -> ")
      group.names0 <- group.names
      group.names <- paste0(group.names0, " (", dataset.name[comparison[ii]], ")")

      if (nrow(df.net) > 0) {
        df.net$pval[df.net$pval > 0.05] = 1
        df.net$pval[df.net$pval > 0.01 & df.net$pval <= 0.05] = 2
        df.net$pval[df.net$pval <= 0.01] = 3
        df.net$prob[df.net$prob == 0] <- NA
        df.net$prob.original <- df.net$prob
        df.net$prob <- -1/log(df.net$prob)
      } else {
        df.net <- as.data.frame(matrix(NA, nrow = length(group.names), ncol = 5))
        colnames(df.net) <- c("interaction_name_2","source.target","prob","pval","prob.original")
        df.net$source.target <- group.names0
      }
      # df.net$group.names <- sub(paste0(' \\(',dataset.name[comparison[ii]],'\\)'),'',as.character(df.net$source.target))
      df.net$group.names <- as.character(df.net$source.target)
      df.net$source.target <- paste0(df.net$source.target, " (", dataset.name[comparison[ii]], ")")
      df.net$dataset <- dataset.name[comparison[ii]]
      df.all <- rbind(df.all, df.net)
    }
    if (nrow(df.all) == 0) {
      stop("No interactions are detected. Please consider changing the cell groups for analysis. ")
    }

    idx1 <- which(is.infinite(df.all$prob) | df.all$prob < 0)
    if (sum(idx1) > 0) {
      values.assign <- seq(max(df.all$prob, na.rm = T)*1.1, max(df.all$prob, na.rm = T)*1.5, length.out = length(idx1))
      position <- sort(df.all$prob.original[idx1], index.return = TRUE)$ix
      df.all$prob[idx1] <- values.assign[match(1:length(idx1), position)]
    }

    df.all$interaction_name_2[is.na(df.all$interaction_name_2)] <- df.all$interaction_name_2[!is.na(df.all$interaction_name_2)][1]

    df <- df.all
    df <- order(df, decreasing = TRUE)
    df <- with(df, df[order(interaction_name_2),])
    #df <- with(df, df[interaction_name_2,])
    df$interaction_name_2 <- factor(df$interaction_name_2, levels = unique(df$interaction_name_2))

    cells.order <- c()
    dataset.name.order <- c()
    for (i in 1:length(group.names0)) {
      for (j in 1:length(comparison)) {
        cells.order <- c(cells.order, paste0(group.names0[i], " (", dataset.name[comparison[j]], ")"))
        dataset.name.order <- c(dataset.name.order, dataset.name[comparison[j]])
      }
    }
    df$source.target <- factor(df$source.target, levels = cells.order)
  }

  min.cutoff <- quantile(df$prob, min.quantile,na.rm= T)
  max.cutoff <- quantile(df$prob, max.quantile,na.rm= T)
  df$prob[df$prob < min.cutoff] <- min.cutoff
  df$prob[df$prob > max.cutoff] <- max.cutoff


  if (remove.isolate) {
    df <- df[!is.na(df$prob), ]
    line.on <- FALSE
  }
  if (!is.null(max.dataset)) {
    # line.on <- FALSE
    # df <- df[!is.na(df$prob),]
    signaling <- as.character(sort(unique(df$interaction_name_2))) 
    for (i in signaling) {
      df.i <- df[df$interaction_name_2 == i, ,drop = FALSE]
      cell <- as.character(sort(unique(df.i$group.names))) 
      for (j in cell) {
        df.i.j <- df.i[df.i$group.names == j, , drop = FALSE]
        values <- df.i.j$prob
        idx.max <- which(values == max(values, na.rm = T))
        idx.min <- which(values == min(values, na.rm = T))
        #idx.na <- c(which(is.na(values)), which(!(dataset.name[comparison] %in% df.i.j$dataset)))
        dataset.na <- c(df.i.j$dataset[is.na(values)], setdiff(dataset.name[comparison], df.i.j$dataset))
        if (length(idx.max) > 0) {
          if (!(df.i.j$dataset[idx.max] %in% dataset.name[max.dataset])) {
            df.i.j$prob <- NA
          } else if ((idx.max != idx.min) & !is.null(min.dataset)) {
            if (!(df.i.j$dataset[idx.min] %in% dataset.name[min.dataset])) {
              df.i.j$prob <- NA
            } else if (length(dataset.na) > 0 & sum(!(dataset.name[min.dataset] %in% dataset.na)) > 0) {
              df.i.j$prob <- NA
            }
          }
        }
        df.i[df.i$group.names == j, "prob"] <- df.i.j$prob
      }
      df[df$interaction_name_2 == i, "prob"] <- df.i$prob
    }
    #df <- df[!is.na(df$prob), ]
  }
  if (remove.isolate) {
    df <- df[!is.na(df$prob), ]
    line.on <- FALSE
  }
  if (nrow(df) == 0) {
    stop("No interactions are detected. Please consider changing the cell groups for analysis. ")
  }
  df$interaction_name_2 <- factor(df$interaction_name_2, levels = unique(df$interaction_name_2))
  #df$interaction_name_2 <- factor(sort(df$interaction_name_2, levels = unique(df$interaction_name_2))) 
  df$source.target = droplevels(df$source.target, exclude = setdiff(levels(df$source.target),unique(df$source.target)))

  g <- ggplot(df, aes(x = source.target, y = interaction_name_2, color = prob, size = pval)) +
    geom_point(pch = 16) +
    theme_linedraw() + theme(panel.grid.major = element_blank()) +
    theme(axis.text.x = element_text(angle = angle.x, hjust= hjust.x, vjust = vjust.x),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    scale_x_discrete(position = "bottom")

  values <- c(1,2,3); names(values) <- c("p > 0.05", "0.01 < p < 0.05","p < 0.001")
  g <- g + scale_radius(range = c(min(df$pval), max(df$pval)), breaks = sort(unique(df$pval)),labels = names(values)[values %in% sort(unique(df$pval))], name = "p-value")
  #g <- g + scale_radius(range = c(1,3), breaks = values,labels = names(values), name = "p-value")
  if (min(df$prob, na.rm = T) != max(df$prob, na.rm = T)) {
    g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99), na.value = "white", limits=c(quantile(df$prob, 0,na.rm= T), quantile(df$prob, 1,na.rm= T)),
                                    breaks = c(quantile(df$prob, 0,na.rm= T), quantile(df$prob, 1,na.rm= T)), labels = c("min","max")) +
      guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob."))
  } else {
    g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99), na.value = "white") +
      guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob."))
  }

  g <- g + theme(text = element_text(size = font.size),plot.title = element_text(size=font.size.title)) +
    theme(legend.title = element_text(size = 8), legend.text = element_text(size = 6))

  if (grid.on) {
    if (length(unique(df$source.target)) > 1) {
      g <- g + geom_vline(xintercept=seq(1.5, length(unique(df$source.target))-0.5, 1),lwd=0.1,colour=color.grid)
    }
    if (length(unique(df$interaction_name_2)) > 1) {
      g <- g + geom_hline(yintercept=seq(1.5, length(unique(df$interaction_name_2))-0.5, 1),lwd=0.1,colour=color.grid)
    }
  }
  if (!is.null(title.name)) {
    g <- g + ggtitle(title.name) + theme(plot.title = element_text(hjust = 0.5))
  }

  if (!is.null(comparison)) {
    if (line.on) {
      xintercept = seq(0.5+length(dataset.name[comparison]), length(group.names0)*length(dataset.name[comparison]), by = length(dataset.name[comparison]))
      g <- g + geom_vline(xintercept = xintercept, linetype="dashed", color = "grey60", size = line.size)
    }
    if (color.text.use) {
      if (is.null(group)) {
        group <- 1:length(comparison)
        names(group) <- dataset.name[comparison]
      }
      if (is.null(color.text)) {
        color <- ggPalette(length(unique(group)))
      } else {
        color <- color.text
      }
      names(color) <- names(group[!duplicated(group)])
      color <- color[group]
      #names(color) <- dataset.name[comparison]
      dataset.name.order <- levels(df$source.target)
      dataset.name.order <- stringr::str_match(dataset.name.order, "\\(.*\\)")
      dataset.name.order <- stringr::str_sub(dataset.name.order, 2, stringr::str_length(dataset.name.order)-1)
      xtick.color <- color[dataset.name.order]
      g <- g + theme(axis.text.x = element_text(colour = xtick.color))
    }
  }
  if (!show.legend) {
    g <- g + theme(legend.position = "none")
  }
  if (return.data) {
    return(list(communication = df, gg.obj = g))
  } else {
    return(g)
  }

}




netVisual_bubble(data, sources.use = 'vEVT', targets.use= c('dFB', 'dVEC', 'dLEC', 'DSC_1', 'DSC_2', 'dMSC'), 
                 remove.isolate= TRUE, thresh = 0.001, color.heatmap = "Spectral", max.quantile=0.8,
                signaling = c("FN1", "COLLAGEN", "LAMININ", "DESMOSOME", 
                                                       "VEGF", "VISFATIN", "PTN",  "CALCR", 
                                                      "ADGRE5", "APP", "CDH5", "EPHA", "NECTIN"))


