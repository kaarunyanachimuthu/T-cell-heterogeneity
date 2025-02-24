setwd("C:/Users/kaaru/Desktop/Assignment/Internship/Nick")
options(future.globals.maxSize = 100000 * 1024^2)

library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(patchwork)
library(dplyr)
library(ggplot2)

set.seed(1)

Cd4 <- readRDS("./cd4-1.rds")
Cd4.cue <- subset(Cd4, subset = treatment == "CUE")
DefaultAssay(Cd4.cue) <- "RNA"
t.list <- SplitObject(Cd4.cue, split.by = "orig.ident")

# Run the standard workflow for visualization and clustering
for (i in 1:length(t.list)) {
  t.list[[i]] <- NormalizeData(t.list[[i]], verbose = FALSE)
  t.list[[i]] <- FindVariableFeatures(t.list[[i]], selection.method = "vst",
                                      nfeatures = 2000, verbose = FALSE)
}
features <- SelectIntegrationFeatures(object.list = t.list)
t.list <- lapply(X = t.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
anchors <- FindIntegrationAnchors(object.list = t.list, reduction = "rpca",
                                  dims = 1:10)
Cd4.cue <- IntegrateData(anchorset = anchors, dims = 1:10, k.weight = 50)
Cd4.cue[["RNA"]] <- JoinLayers(Cd4.cue[["RNA"]])
Cd4.cue <- ScaleData(Cd4.cue, features = rownames(Cd4.cue), verbose = FALSE)
Cd4.cue <- RunPCA(Cd4.cue)
ElbowPlot(Cd4.cue)
Cd4.cue <- FindNeighbors(Cd4.cue, reduction = "pca", dims = 1:10)
Cd4.cue <- FindClusters(Cd4.cue, resolution = 0.5)
Cd4.cue <- RunUMAP(Cd4.cue, reduction = "pca", dims = 1:10)

DimPlot(Cd4.cue, reduction = "umap", label=T, pt.size = 1.5) 
itreg.markers <- FindMarkers(Cd4.cue, ident.1 = 6, only.pos = TRUE)
write.csv(itreg.markers, file = "itreg markers.csv")



cds <- as.cell_data_set(Cd4.cue)

fData(cds)
rownames(fData(cds))[1:10]

fData(cds)$gene_short_name <- rownames(fData(cds))

cd4.partition <- c(rep(1,length(cds@colData@rownames)))
names(cd4.partition) <- cds@colData@rownames
cd4.partition <- as.factor(cd4.partition)
cds@clusters$UMAP$partitions <- cd4.partition

cds@int_colData@listData$reducedDims$UMAP <- Cd4.cue@reductions$umap@cell.embeddings

p1 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE, label_groups_by_cluster = FALSE, group_label_size = 5)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE, label_groups_by_cluster = FALSE, group_label_size = 5)
wrap_plots(p1, p2)

integrated.sub <- subset(as.Seurat(cds, assay = NULL), monocle3_partitions == 1)
cds <- as.cell_data_set(integrated.sub)

cds <- learn_graph(cds, use_partition = TRUE, verbose = FALSE)

plot_cells(cds,
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           group_label_size = 4)

cds <- order_cells(cds, reduction_method = "UMAP", root_cells = colnames(cds[,clusters(cds) == 3]))
plot_cells(cds,
           color_cells_by = "pseudotime",
           group_cells_by = "cluster",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           trajectory_graph_color = "grey60",
           trajectory_graph_segment_size = 1.5,
           cell_size = 1)

# ordering cells by monocle3 pseudotime
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))
ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime, median), fill = seurat_clusters)) +
  geom_boxplot()

integrated.sub <- as.Seurat(cds, assay = NULL)
FeaturePlot(integrated.sub, "monocle3_pseudotime", label = TRUE)

cds_graph_test_results <- graph_test(cds,
                                     neighbor_graph = "principal_graph",
                                     cores = 8)
rowData(cds)$gene_short_name <- row.names(rowData(cds))

head(cds_graph_test_results, error=FALSE, message=FALSE, warning=FALSE)

deg_ids <- rownames(subset(cds_graph_test_results[order(cds_graph_test_results$morans_I, decreasing = TRUE),], q_value < 0.05))

plot_cells(cds,
           genes=c(head(deg_ids)),
           show_trajectory_graph = FALSE,
           label_cell_groups = FALSE,
           label_leaves = FALSE, cell_size = 1)
