setwd("C:/Users/kaaru/Desktop/Assignment/Internship/Nick")
options(future.globals.maxSize = 100000 * 1024^2)
library(dplyr)
library(Seurat)
library(cowplot)
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(fgsea)
library(msigdbr)
library(presto)
library(tibble)
library(patchwork)
library(SingleR)
library(celldex)
library(viridis)
library(EnhancedVolcano)
# Create data directory
data_dir <- file.path(getwd(), "Data")
dir.create(data_dir, showWarnings = FALSE)

# Define data paths
data1 <- file.path(data_dir, "F2_PBS_out", "sample_filtered_feature_bc_matrix")
data2 <- file.path(data_dir, "F3_CUE_out", "sample_filtered_feature_bc_matrix")
data3 <- file.path(data_dir, "F4_PBS_out", "sample_filtered_feature_bc_matrix")
data4 <- file.path(data_dir, "F9_CUE_out", "sample_filtered_feature_bc_matrix")
data5 <- file.path(data_dir, "PBS1")
data6 <- file.path(data_dir, "CUE1")
data7 <- file.path(data_dir, "PBS2")
data8 <- file.path(data_dir, "CUE2")

# Load data
dataPBS <- Read10X(data.dir = data1)
dataCUE <- Read10X(data.dir = data2)
dataPBS0 <- Read10X(data.dir = data3)
dataCUE0 <- Read10X(data.dir = data4)
dataPBS1 <- Read10X(data.dir = data5)
dataCUE1 <- Read10X(data.dir = data6)
dataPBS2 <- Read10X(data.dir = data7)
dataCUE2 <- Read10X(data.dir = data8)

F2.PBS <- CreateSeuratObject(counts = dataPBS, project = "F2_PBS", min.cells = 5)
F2.PBS$treatment <- "PBS"
F2.PBS$sample <- "F2"
F2.PBS[["percent.mt.F2PBS"]]<-PercentageFeatureSet(F2.PBS, pattern="^mt-")
VlnPlot(F2.PBS, features=c("nFeature_RNA", "nCount_RNA", "percent.mt.F2PBS"))
F2.PBS<-subset(F2.PBS, subset = nFeature_RNA>500 & percent.mt.F2PBS<10 & nCount_RNA<20000)

F3.CUE <- CreateSeuratObject(counts = dataCUE, project = "F3_CUE", min.cells = 5)
F3.CUE$treatment <- "CUE"
F3.CUE$sample <- "F3"
F3.CUE[["percent.mt.F3CUE"]]<-PercentageFeatureSet(F3.CUE, pattern="^mt-")
VlnPlot(F3.CUE, features=c("nFeature_RNA", "nCount_RNA", "percent.mt.F3CUE"))
F3.CUE<-subset(F3.CUE, subset = nFeature_RNA>500 & percent.mt.F3CUE<10 & nCount_RNA<20000)

F4.PBS <- CreateSeuratObject(counts = dataPBS0, project = "F4_PBS", min.cells = 5)
F4.PBS$treatment <- "PBS"
F4.PBS$sample <- "F4"
F4.PBS[["percent.mt.F4PBS"]]<-PercentageFeatureSet(F4.PBS, pattern="^mt-")
VlnPlot(F4.PBS, features=c("nFeature_RNA", "nCount_RNA", "percent.mt.F4PBS"))
F4.PBS<-subset(F4.PBS, subset = nFeature_RNA>500 & percent.mt.F4PBS<10 & nCount_RNA<20000)

F9.CUE <- CreateSeuratObject(counts = dataCUE0, project = "F9_CUE", min.cells = 5)
F9.CUE$treatment <- "CUE"
F9.CUE$sample <- "F9"
F9.CUE[["percent.mt.F9CUE"]]<-PercentageFeatureSet(F9.CUE, pattern="^mt-")
VlnPlot(F9.CUE, features=c("nFeature_RNA", "nCount_RNA", "percent.mt.F9CUE"))
F9.CUE<-subset(F9.CUE, subset = nFeature_RNA>500 & percent.mt.F9CUE<10 & nCount_RNA<20000)

D1.PBS <- CreateSeuratObject(counts = dataPBS1$`Gene Expression`, project = "D1_PBS", min.cells = 5)
D1.PBS$treatment <- "PBS"
D1.PBS$sample <- "D1"
D1.PBS[["percent.mt.D1PBS"]]<-PercentageFeatureSet(D1.PBS, pattern="^mt-")
VlnPlot(D1.PBS, features=c("nFeature_RNA", "nCount_RNA", "percent.mt.D1PBS"))
D1.PBS<-subset(D1.PBS, subset = nFeature_RNA>500 & percent.mt.D1PBS<10 & nCount_RNA<20000)

D2.CUE <- CreateSeuratObject(counts = dataCUE1$`Gene Expression`, project = "D2_CUE", min.cells = 5)
D2.CUE$treatment <- "CUE"
D2.CUE$sample <- "D2"
D2.CUE[["percent.mt.D2CUE"]]<-PercentageFeatureSet(D2.CUE, pattern="^mt-")
VlnPlot(D2.CUE, features=c("nFeature_RNA", "nCount_RNA", "percent.mt.D2CUE"))
D2.CUE<-subset(D2.CUE, subset = nFeature_RNA>500 & percent.mt.D2CUE<10 & nCount_RNA<20000)

D3.PBS <- CreateSeuratObject(counts = dataPBS2$`Gene Expression`, project = "D3_PBS", min.cells = 5)
D3.PBS$treatment <- "PBS"
D3.PBS$sample <- "D3"
D3.PBS[["percent.mt.D3PBS"]]<-PercentageFeatureSet(D3.PBS, pattern="^mt-")
VlnPlot(D3.PBS, features=c("nFeature_RNA", "nCount_RNA", "percent.mt.D3PBS"))
D3.PBS<-subset(D3.PBS, subset = nFeature_RNA>500 & percent.mt.D3PBS<10 & nCount_RNA<20000)

D4.CUE <- CreateSeuratObject(counts = dataCUE2$`Gene Expression`, project = "D4_CUE", min.cells = 5)
D4.CUE$treatment <- "CUE"
D4.CUE$sample <- "D4"
D4.CUE[["percent.mt.D4CUE"]]<-PercentageFeatureSet(D4.CUE, pattern="^mt-")
VlnPlot(D4.CUE, features=c("nFeature_RNA", "nCount_RNA", "percent.mt.D4CUE"))
D4.CUE<-subset(D4.CUE, subset = nFeature_RNA>500 & percent.mt.D4CUE<10 & nCount_RNA<20000)

F2.PBS <- RenameCells(object = F2.PBS, add.cell.id = "F2_PBS")
F3.CUE <- RenameCells(object = F3.CUE, add.cell.id = "F3_CUE")
F4.PBS <- RenameCells(object = F4.PBS, add.cell.id = "F4_PBS")
F9.CUE <- RenameCells(object = F9.CUE, add.cell.id = "F9_CUE")
D1.PBS <- RenameCells(object = D1.PBS, add.cell.id = "D1_PBS")
D2.CUE <- RenameCells(object = D2.CUE, add.cell.id = "D2_CUE")
D3.PBS <- RenameCells(object = D3.PBS, add.cell.id = "D3_PBS")
D4.CUE <- RenameCells(object = D4.CUE, add.cell.id = "D4_CUE")


obj.combined <- merge(F2.PBS, y = c(F3.CUE,F4.PBS,F9.CUE,D1.PBS,D2.CUE,D3.PBS,D4.CUE), project = "Cue_treatment")
obj.combined <- NormalizeData(obj.combined, normalization.method = "LogNormalize", scale.factor = 10000)
obj.combined <- FindVariableFeatures(obj.combined, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(obj.combined)
obj.combined <- ScaleData(obj.combined, features = all.genes)
obj.combined <- RunPCA(obj.combined, features = VariableFeatures(object = obj.combined), nfeatures.print = 5, npcs = 60)
print(obj.combined[["pca"]], dims = 1, nfeatures = 10)
integrated.combined <- IntegrateLayers(
  object = obj.combined, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = FALSE
)
integrated.combined[["RNA"]] <- JoinLayers(integrated.combined[["RNA"]])
ElbowPlot(obj.combined, ndims = 60)
obj.combined <- FindNeighbors(integrated.combined,reduction = "integrated.rpca", dims = 1:30)
obj.combined <- FindClusters(obj.combined, resolution = 2, algorithm = 2)
obj.combined <- RunUMAP(obj.combined, dims = 1:30,reduction = "integrated.rpca")
DimPlot(obj.combined, label = TRUE, reduction = "umap", split.by = "treatment")

saveRDS(obj.combined, file = "./reanalysis_clustered.rds")

# Marker identification and clustering
markers <- FindAllMarkers(obj.combined, only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
write.csv(markers, './reanalysis_markers.csv')

# Cell annotation using singleR
ref_immgen <- celldex::ImmGenData()
predictions_main = SingleR(test = GetAssayData(obj.combined), 
                           ref = ref_immgen,
                           labels = ref_immgen$label.main)
predictions_fine = SingleR(test = GetAssayData(obj.combined), 
                           ref = ref_immgen,
                           labels = ref_immgen$label.fine)
plotScoreHeatmap(predictions_main)

obj.combined[['immgen_singler_main']] = rep('NA', ncol(obj.combined))
obj.combined$immgen_singler_main[rownames(predictions_main)] = predictions_main$labels

obj.combined[['immgen_singler_fine']] = rep('NA', ncol(obj.combined))
obj.combined$immgen_singler_fine[rownames(predictions_fine)] = predictions_fine$labels

ggplot(obj.combined[[]], aes(x = orig.ident, fill = immgen_singler_main)) + geom_bar(position = "fill") + scale_fill_viridis(discrete = TRUE)
DimPlot(obj.combined, group.by = c("immgen_singler_main"), label = TRUE, split.by = "treatment")


DimPlot(obj.combined, group.by = "seurat_clusters", label = TRUE) + 
  DimPlot(obj.combined, group.by = "immgen_singler_main", label = TRUE) +
  FeaturePlot(obj.combined, features = "Ptprc")

# Cluster markers
cluster8.markers <- FindMarkers(obj.combined, ident.1 = 8, logfc.threshold = 0.6, test.use = "wilcox", only.pos = TRUE, min.pct = 0.25)
cluster16.markers <- FindMarkers(obj.combined, ident.1 = 16, logfc.threshold = 0.6, test.use = "wilcox", only.pos = TRUE, min.pct = 0.25)
cluster31.markers <- FindMarkers(obj.combined, ident.1 = 31, logfc.threshold = 0.6, test.use = "wilcox", only.pos = TRUE, min.pct = 0.25)
cluster32.markers <- FindMarkers(obj.combined, ident.1 = 32, logfc.threshold = 0.6, test.use = "wilcox", only.pos = TRUE, min.pct = 0.25)

# DE Analysis
# 8
eight<-subset(obj.combined, idents = 8)
Idents(eight)<-"treatment"
avg.eight <- log1p(AverageExpression(eight, verbose = FALSE)$RNA)
avg.eight <-data.frame(avg.eight)
avg.eight$gene <- rownames(avg.eight)
eight_comp_p2<-ggplot(avg.eight, aes(PBS, CUE)) + geom_point() + geom_smooth(method = "lm", se = FALSE, color = "gray") + ggtitle("Cluster 8 DE Analysis") + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), plot.title = element_text(hjust = 0.5)) + xlab("PBS Treated") + ylab("CUE Treated")
eight_comp_p2<-LabelPoints(plot = eight_comp_p2, points = c("Trbv3","Ecm1","Pparg","Plcl1","Rnf128","Ly6a","Trbv13-3","Uros","Gpr15","Trbv19"), color="blue", repel = TRUE, xnudge = 0, ynudge = 0)
plot_grid(eight_comp_p2)

clus8_diff<-FindMarkers(eight, ident.1="PBS",ident.2="CUE", verbose = FALSE, logfc.threshold = 0.6, min.pct = 0.1)

DoHeatmap(eight, features = c(c("Ly6a","Trbv13-3","Uros","Gpr15","Trbv19"),c("Trbv3","Ecm1","Pparg","Plcl1","Rnf128"))) +
  scale_fill_gradientn(colors = c("blue", "white", "red"))

# 16
sixteen<-subset(obj.combined, idents = 16)
Idents(sixteen)<-"treatment"
avg.sxt <- log1p(AverageExpression(sixteen, verbose = FALSE)$RNA)
avg.sxt <-data.frame(avg.sxt)
avg.sxt$gene <- rownames(avg.sxt)
sxt_comp_p2<-ggplot(avg.sxt, aes(PBS, CUE)) + geom_point() + geom_smooth(method = "lm", se = FALSE, color = "gray") + ggtitle("Cluster 16 DE Analysis") + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), plot.title = element_text(hjust = 0.5)) + xlab("PBS Treated") + ylab("CUE Treated")
sxt_comp_p2<-LabelPoints(plot = sxt_comp_p2, points = c("Ifit3","Ifit3b","Ly6a","A330040F15Rik","Parp11","Ccl3","Ccl4","Adgre5","Nr4a2","Hip1"), color="blue", repel = TRUE, xnudge = 0, ynudge = 0)
plot_grid(sxt_comp_p2)

clus16_diff<-FindMarkers(sixteen, ident.1="PBS",ident.2="CUE", verbose = FALSE, logfc.threshold = 0.6, min.pct = 0.1)

DoHeatmap(sixteen, features = c("Ifit3","Ifit3b","Ly6a","A330040F15Rik","Parp11","Ccl3","Ccl4","Adgre5","Nr4a2","Hip1")) +
  scale_fill_gradientn(colors = c("blue", "white", "red"))

# 31
treg<-subset(obj.combined, idents = 31)
treg <- FindVariableFeatures(treg)
Idents(treg)<-"treatment"
avg.treg <- log1p(AverageExpression(treg, verbose = FALSE)$RNA)
avg.treg <-data.frame(avg.treg)
avg.treg$gene <- rownames(avg.treg)
treg_comp_p2<-ggplot(avg.treg, aes(PBS, CUE)) + geom_point() + geom_smooth(method = "lm", se = FALSE, color = "gray") + ggtitle("Cluster 31 DE Analysis") + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), plot.title = element_text(hjust = 0.5)) + xlab("PBS Treated") + ylab("CUE Treated")
treg_comp_p2<-LabelPoints(plot = treg_comp_p2, points = c("Hsph1","Rtp4","Emb","Tkt","Ifi27l2a","Cytip","Tnfrsf4","Tnfrsf9","Lrrc8c","Il2ra","Ccnd2","Arl5a","Gzmb","Cish","Socs1","Socs2","Hivep3","Nckap5","Gadd45g"), color="blue", repel = TRUE, xnudge = 0, ynudge = 0)
plot_grid(treg_comp_p2)

clus31_diff<-FindMarkers(treg, ident.1="PBS",ident.2="CUE", verbose = FALSE, logfc.threshold = 0.6, min.pct = 0.1)
write.csv(clus31_diff, './clus31.csv')

DoHeatmap(treg, features = c("Hsph1","Rtp4","Emb","Tkt","Ifi27l2a","Cytip","Tnfrsf4","Tnfrsf9","Lrrc8c","Il2ra","Ccnd2","Arl5a","Gzmb","Cish","Socs1","Socs2","Hivep3","Nckap5","Gadd45g")) +
  scale_fill_gradientn(colors = c("blue", "white", "red"))

# 32
tpro<-subset(obj.combined, idents = 32)
Idents(tpro)<-"treatment"
avg.tpro <- log1p(AverageExpression(tpro, verbose = FALSE)$RNA)
avg.tpro <-data.frame(avg.tpro)
avg.tpro$gene <- rownames(avg.tpro)
tpro_comp_p2<-ggplot(avg.tpro, aes(PBS, CUE)) + geom_point() + geom_smooth(method = "lm", se = FALSE, color = "gray") + ggtitle("Cluster 32 DE Analysis") + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), plot.title = element_text(hjust = 0.5)) + xlab("PBS Treated") + ylab("CUE Treated")
tpro_comp_p2<-LabelPoints(plot = tpro_comp_p2, points = c("Ccl4","Trav12-3","B3gnt5","Hbb-bs","H2-Aa","Nckap5","Nccrp1","Rln3","Tspoap1","Foxp3"), color="blue", repel = TRUE, xnudge = 0, ynudge = 0)
plot_grid(tpro_comp_p2)

clus32_diff<-FindMarkers(tpro, ident.1="PBS",ident.2="CUE", verbose = FALSE, logfc.threshold = 0.6, min.pct = 0.1)

DoHeatmap(tpro, features = c("Ccl4","Trav12-3","B3gnt5","Hbb-bs","H2-Aa","Nckap5","Nccrp1","Rln3","Tspoap1","Foxp3","Il1rl1")) +
  scale_fill_gradientn(colors = c("blue", "white", "red"))

subset_tcell <- subset(obj.combined, idents=c(8,16,31,32))

# Generate violin plots
genes <- c("Tgfbr1","Tgfbr2","Tgfbr3")
plots <- list()
for (gene in genes) {
  p <- VlnPlot(subset_tcell, features = gene, split.by = "treatment", cols = c("red","blue")) + ylim(0, 4)
  plots[[gene]] <- p
}
combined_plot <- wrap_plots(plots, ncol = 2)
print(combined_plot)

## Average expression of transcripts from violin plots
subset_tcell$seurat_clusters.treatment<-paste(Idents(subset_tcell), subset_tcell$treatment, sep="_")
subset_tcell$seurat_clusters<- Idents(subset_tcell)
Idents(subset_tcell)<-"seurat_clusters.treatment"

AverageExpression(subset_tcell, features = "Ifng")


# B Cell analysis
bcell <- subset(obj.combined, idents = c(0,1,3,6,10,11))
bcell.markers <- FindAllMarkers(bcell,only.pos = TRUE)

Idents(bcell)<-"treatment"
avg.bcell <- log1p(AverageExpression(bcell, verbose = FALSE)$RNA)
avg.bcell <-data.frame(avg.bcell)
avg.bcell$gene <- rownames(avg.bcell)
b_comp_p2<-ggplot(avg.bcell, aes(PBS, CUE)) + geom_point() + geom_smooth(method = "lm", se = FALSE, color = "gray") + ggtitle("B cell DE Analysis") + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), plot.title = element_text(hjust = 0.5)) + xlab("PBS Treated") + ylab("CUE Treated")
b_comp_p2<-LabelPoints(plot = b_comp_p2, points = c("Hbb-bs","Gm9750","Cfp","Fam129c","Gm10282","Wwox","Adamts13","Il31ra","Hist2h2aa1","Tmtc2"), color="blue", repel = TRUE, xnudge = 0, ynudge = 0)
plot_grid(b_comp_p2)

b_diff<-FindMarkers(bcell, ident.1="PBS",ident.2="CUE", verbose = FALSE, logfc.threshold = 0.6, min.pct = 0.1)

b_de_sig <- b_diff[b_diff$p_val_adj < 0.05,]
b_de_sig %>%
  top_n(n = 5, wt = -avg_log2FC)
b_de_sig %>%
  top_n(n = 5, wt = avg_log2FC)

EnhancedVolcano(b_diff,
                lab = rownames(b_diff),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'B Cell Differential Expression Analysis',
                pCutoff = 0.01,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 3.0,
                colAlpha = 1,
                col = c('black','pink','purple','red3'),
                boxedLabels = TRUE,
                parseLabels = TRUE,
                drawConnectors = TRUE,
                legendPosition = 'bottom',
) 

DoHeatmap(bcell, features = c("Hbb-bs","Gm9750","Cfp","Fam129c","Gm10282","Wwox","Adamts13","Il31ra","Hist2h2aa1","Tmtc2")) +
  scale_fill_gradientn(colors = c("blue", "white", "red"))


# Generate violin plots
genes <- c("Wwox","Adamts13","Il31ra","Hist2h2aa1","Tmtc2","Vpreb3")
plots <- list()
for (gene in genes) {
  p <- VlnPlot(bcell, features = gene, split.by = "treatment", cols = c("red","blue")) + ylim(0, 4)
  plots[[gene]] <- p
}
combined_plot <- wrap_plots(plots, ncol = 2)
print(combined_plot)

# NK Cell analysis
nkcell <- subset(obj.combined, idents = c(4,7,27))
nkcell.markers <- FindAllMarkers(nkcell,only.pos = TRUE)

Idents(nkcell)<-"treatment"
avg.nkcell <- log1p(AverageExpression(nkcell, verbose = FALSE)$RNA)
avg.nkcell <-data.frame(avg.nkcell)
avg.nkcell$gene <- rownames(avg.nkcell)
nk_comp_p2<-ggplot(avg.nkcell, aes(PBS, CUE)) + geom_point() + geom_smooth(method = "lm", se = FALSE, color = "gray") + ggtitle("NK cell DE Analysis") + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), plot.title = element_text(hjust = 0.5)) + xlab("PBS Treated") + ylab("CUE Treated")
nk_comp_p2<-LabelPoints(plot = nk_comp_p2, points = c("Cx3cr1","Steap3","Serpinb9b","Klhl30","Ly6c2","Tcrg-C4","Serpina3g","Dtx1","Dscam","Podnl1"), color="blue", repel = TRUE, xnudge = 0, ynudge = 0)
plot_grid(nk_comp_p2)

nk_diff<-FindMarkers(nkcell, ident.1="PBS",ident.2="CUE", verbose = FALSE, logfc.threshold = 0.6, min.pct = 0.1)

EnhancedVolcano(nk_diff,
                lab = rownames(nk_diff),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'NK Cell Differential Expression Analysis',
                pCutoff = 0.001,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 3.0,
                colAlpha = 1,
                col = c('black','pink','purple','red3'),
                boxedLabels = TRUE,
                parseLabels = TRUE,
                drawConnectors = TRUE,
                legendPosition = 'bottom',
) 

DoHeatmap(nkcell, features = c("Cx3cr1","Steap3","Serpinb9b","Klhl30","Ly6c2","Tcrg-C4","Serpina3g","Dtx1","Dscam","Podnl1")) +
  scale_fill_gradientn(colors = c("blue", "white", "red"))

# Generate violin plots
genes <- c("Tcrg-C4","Serpina3g","Dtx1","Dscam","Podnl1")
plots <- list()
for (gene in genes) {
  p <- VlnPlot(nkcell, features = gene, split.by = "treatment", cols = c("red","blue")) + ylim(0, 4)
  plots[[gene]] <- p
}
combined_plot <- wrap_plots(plots, ncol = 2)
print(combined_plot)

# Mast Cell Analysis
mastcell <- subset(obj.combined, idents = c(37))

Idents(mastcell)<-"treatment"
avg.mastcell <- log1p(AverageExpression(mastcell, verbose = FALSE)$RNA)
avg.mastcell <-data.frame(avg.mastcell)
avg.mastcell$gene <- rownames(avg.mastcell)
mast_comp_p2<-ggplot(avg.mastcell, aes(PBS, CUE)) + geom_point() + geom_smooth(method = "lm", se = FALSE, color = "gray") + ggtitle("Mast cell DE Analysis") + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), plot.title = element_text(hjust = 0.5)) + xlab("PBS Treated") + ylab("CUE Treated")
mast_comp_p2<-LabelPoints(plot = mast_comp_p2, points = c("Dusp6","Ifitm3","Ifitm1","H2-K1","Prss34","Rps28","Plin3","Micu2","Rpl36a","Nme3"), color="blue", repel = TRUE, xnudge = 0, ynudge = 0)
plot_grid(mast_comp_p2)

mast_diff<-FindMarkers(mastcell, ident.1="PBS",ident.2="CUE", verbose = FALSE, logfc.threshold = 0.6, min.pct = 0.1)

EnhancedVolcano(mast_diff,
                lab = rownames(mast_diff),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'Mast Cell Differential Expression Analysis',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 3.0,
                colAlpha = 1,
                col = c('black','pink','purple','red3'),
                boxedLabels = TRUE,
                parseLabels = TRUE,
                drawConnectors = TRUE,
                legendPosition = 'bottom',
) 

DoHeatmap(mastcell, features = c("Dusp6","Ifitm3","Ifitm1","H2-K1","Prss34","Rps28","Plin3","Micu2","Rpl36a","Nme3")) +
  scale_fill_gradientn(colors = c("blue", "white", "red"))

# Macrophage Analysis
macrophage <- subset(obj.combined, idents = c(14, 29))

Idents(macrophage)<-"treatment"
avg.macrophage <- log1p(AverageExpression(macrophage, verbose = FALSE)$RNA)
avg.macrophage <-data.frame(avg.macrophage)
avg.macrophage$gene <- rownames(avg.macrophage)
macro_comp_p2<-ggplot(avg.macrophage, aes(PBS, CUE)) + geom_point() + geom_smooth(method = "lm", se = FALSE, color = "gray") + ggtitle("Macrophage DE Analysis") + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), plot.title = element_text(hjust = 0.5)) + xlab("PBS Treated") + ylab("CUE Treated")
macro_comp_p2<-LabelPoints(plot = macro_comp_p2, points = c("Ly6a","Acod1","C5ar1","Slfn4","Ifit2","Tns1","Trerf1","Mras","Lag3","Cd63"), color="blue", repel = TRUE, xnudge = 0, ynudge = 0)
plot_grid(macro_comp_p2)

macro_diff<-FindMarkers(macrophage, ident.1="PBS",ident.2="CUE", verbose = FALSE, logfc.threshold = 0.6, min.pct = 0.1)
write.csv(macro_diff, "./macrodiff.csv")
EnhancedVolcano(macro_diff,
                lab = rownames(macro_diff),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'Macrophage Differential Expression Analysis',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 3.0,
                colAlpha = 1,
                col = c('black','pink','purple','red3'),
                boxedLabels = TRUE,
                parseLabels = TRUE,
                drawConnectors = TRUE,
                legendPosition = 'bottom',
) 

DoHeatmap(macrophage, features = c("Ly6a","Acod1","C5ar1","Slfn4","Ifit2","Tns1","Trerf1","Mras","Lag3","Cd63")) +
  scale_fill_gradientn(colors = c("blue", "white", "red"))

macrophage$treatment <- factor(macrophage$treatment, levels = c("PBS", "CUE"))
VlnPlot(macrophage, features = c("Ly6a","Acod1","C5ar1","Slfn4","Ifit2","Tns1","Trerf1","Mras","Lag3","Cd63"), split.by = "treatment")

# Neutrophil Analysis
neutrophil <- subset(obj.combined, idents = c(2, 25, 42))

Idents(neutrophil)<-"treatment"
avg.neutro <- log1p(AverageExpression(neutrophil, verbose = FALSE)$RNA)
avg.neutro <-data.frame(avg.neutro)
avg.neutro$gene <- rownames(avg.neutro)
neutro_comp_p2<-ggplot(avg.neutro, aes(PBS, CUE)) + geom_point() + geom_smooth(method = "lm", se = FALSE, color = "gray") + ggtitle("Neutrophil DE Analysis") + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), plot.title = element_text(hjust = 0.5)) + xlab("PBS Treated") + ylab("CUE Treated")
neutro_comp_p2<-LabelPoints(plot = neutro_comp_p2, points = c("Gbp4","Igfbp6","Trim30c","Scimp","Ggct","Orm1","Tmem45a2","Astn2","Lhfpl2","Ly6g5b"), color="blue", repel = TRUE, xnudge = 0, ynudge = 0)
plot_grid(neutro_comp_p2)

neutro_diff<-FindMarkers(neutrophil, ident.1="PBS",ident.2="CUE", verbose = FALSE, logfc.threshold = 0.6, min.pct = 0.1)
write.csv(macro_diff, "./macrodiff.csv")
EnhancedVolcano(neutro_diff,
                lab = rownames(neutro_diff),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'Neutrophil Differential Expression Analysis',
                pCutoff = 0.001,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 3.0,
                colAlpha = 1,
                col = c('black','pink','purple','red3'),
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                legendPosition = 'bottom',
) 

DoHeatmap(neutrophil, features = c("Gbp4","Igfbp6","Trim30c","Scimp","Ggct","Orm1","Tmem45a2","Astn2","Lhfpl2","Ly6g5b")) +
  scale_fill_gradientn(colors = c("blue", "white", "red"))

macrophage$treatment <- factor(macrophage$treatment, levels = c("PBS", "CUE"))
VlnPlot(macrophage, features = c("Ly6a","Acod1","C5ar1","Slfn4","Ifit2","Tns1","Trerf1","Mras","Lag3","Cd63"), split.by = "treatment")

# Generate violin plots
genes <- c("Orm1","Tmem45a2","Astn2","Lhfpl2","Ly6g5b")
plots <- list()
for (gene in genes) {
  p <- VlnPlot(neutrophil, features = gene, split.by = "treatment", cols = c("red","blue")) + ylim(0, 4)
  plots[[gene]] <- p
}
combined_plot <- wrap_plots(plots, ncol = 2)
print(combined_plot)

# DC Analysis
dc <- subset(obj.combined, idents = c(40,39))

Idents(dc)<-"treatment"
avg.dc <- log1p(AverageExpression(dc, verbose = FALSE)$RNA)
avg.dc <-data.frame(avg.dc)
avg.dc$gene <- rownames(avg.dc)
dc_comp_p2<-ggplot(avg.dc, aes(PBS, CUE)) + geom_point() + geom_smooth(method = "lm", se = FALSE, color = "gray") + ggtitle("Dendritic Cell DE Analysis") + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), plot.title = element_text(hjust = 0.5)) + xlab("PBS Treated") + ylab("CUE Treated")
dc_comp_p2<-LabelPoints(plot = dc_comp_p2, points = c("Serpinb1a","Slc30a6","Cxcl9","Nup188","Ifit2","Cetn4","Egr3","Tmem158","Pla2g4f","E2f5"), color="blue", repel = TRUE, xnudge = 0, ynudge = 0)
plot_grid(dc_comp_p2)

dc_diff<-FindMarkers(dc, ident.1="PBS",ident.2="CUE", verbose = FALSE, logfc.threshold = 0.6, min.pct = 0.1)
write.csv(macro_diff, "./macrodiff.csv")
EnhancedVolcano(dc_diff,
                lab = rownames(dc_diff),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'Dendritic Cell Differential Expression Analysis',
                pCutoff = 0.001,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 3.0,
                colAlpha = 1,
                col = c('black','pink','purple','red3'),
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                legendPosition = 'bottom',
) 

DoHeatmap(dc, features = c("Serpinb1a","Slc30a6","Cxcl9","Nup188","Ifit2","Cetn4","Egr3","Tmem158","Pla2g4f","E2f5")) +
  scale_fill_gradientn(colors = c("blue", "white", "red"))

# ILC Analysis
ilc <- subset(obj.combined, idents = c(9, 33, 34))

Idents(ilc)<-"treatment"
avg.ilc <- log1p(AverageExpression(ilc, verbose = FALSE)$RNA)
avg.ilc <-data.frame(avg.ilc)
avg.ilc$gene <- rownames(avg.ilc)
ilc_comp_p2<-ggplot(avg.ilc, aes(PBS, CUE)) + geom_point() + geom_smooth(method = "lm", se = FALSE, color = "gray") + ggtitle("ILC DE Analysis") + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), plot.title = element_text(hjust = 0.5)) + xlab("PBS Treated") + ylab("CUE Treated")
ilc_comp_p2<-LabelPoints(plot = ilc_comp_p2, points = c("Ly6a","Art2b","Hspa1b","Emp1","Tuba8","Gm36723","Gm47283","Ttn","Ccl5","Nr4a3"), color="blue", repel = TRUE, xnudge = 0, ynudge = 0)
plot_grid(ilc_comp_p2)

ilc_diff<-FindMarkers(ilc, ident.1="PBS",ident.2="CUE", verbose = FALSE, logfc.threshold = 0.6, min.pct = 0.1)
write.csv(ilc_diff, "./CSV/ILCdiff.csv")
EnhancedVolcano(ilc_diff,
                lab = rownames(ilc_diff),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'ILC Differential Expression Analysis',
                pCutoff = 0.001,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 3.0,
                colAlpha = 1,
                col = c('black','pink','purple','red3'),
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                legendPosition = 'bottom',
) 

DoHeatmap(ilc, features = c("Ly6a","Art2b","Hspa1b","Emp1","Tuba8","Gm36723","Gm47283","Ttn","Ccl5","Nr4a3")) +
  scale_fill_gradientn(colors = c("blue", "white", "red"))































