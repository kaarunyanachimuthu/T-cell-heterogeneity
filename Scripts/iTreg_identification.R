library(Seurat)
library(magrittr)
library(ggplot2)

load("./Data/signaturelistGSE109742.Rda")
# Extract signatures that you want in the heatmap
gene_sets <- list(
  iTreg_Sig_4 = signaturelist[["GSE19512_NAUTRAL_VS_INDUCED_TREG_down"]],  
  iTreg_Sig_3 =  signaturelist[["GSE40443_INDUCED_VS_TOTAL_TREG_up"]],  
  iTreg_Sig_2 = signaturelist[["GSE14308_INDUCED_VS_NATURAL_TREG_up"]],  
  iTreg_Sig_1 = signaturelist[["GSE14415_INDUCED_VS_NATURAL_TREG_up"]],
  eTreg_Sig = signaturelist[["GSE19512_NAUTRAL_VS_INDUCED_TREG_up"]]
)

intersect(signaturelist[["GSE14308_INDUCED_VS_NATURAL_TREG_up"]], signaturelist[["GSE14415_INDUCED_VS_NATURAL_TREG_up"]])

main <- subset(obj.combined, seurat_clusters == 31)
main <- main[rowSums(main) != 0, ] # Remove genes with zero expression everywhere
Treg$seurat_clusters.treatment<-paste(Idents(Treg), Treg$treatment, sep="_")
Treg$seurat_clusters<- Idents(Treg)
Idents(Treg)<-"seurat_clusters.treatment"

library(Seurat)
library(pheatmap)

# Compute module scores for each gene signature
Treg <- AddModuleScore(Treg, features = gene_sets, name = "Signature")

# Extract module scores from metadata
score_matrix <- as.data.frame(Treg@meta.data[, grep("^Signature", colnames(Treg@meta.data))])
colnames(score_matrix) <- names(gene_sets)  # Rename columns to match signature names
score_matrix$Cluster <- Treg$seurat_clusters.treatment  # Add cluster information

# Compute average scores per cluster
cluster_avg_scores <- aggregate(. ~ Cluster, data = score_matrix, FUN = mean)
rownames(cluster_avg_scores) <- cluster_avg_scores$Cluster
cluster_avg_scores$Cluster <- NULL  # Remove redundant cluster column

# Scale scores across clusters
cluster_avg_scores <- t(scale(t(cluster_avg_scores)))  # Scale per signature

# Generate heatmap
pheatmap(cluster_avg_scores, cluster_cols = FALSE, color = colorRampPalette(c("blue", "white", "red"))(100), cluster_rows = TRUE)

# Define your custom cluster order (modify as needed)
custom_order <- c("0_PBS", "1_PBS", "2_PBS", "0_CUE","1_CUE","2_CUE")  # Example: Order clusters as 2 → 0 → 3 → 1

# Reorder rows based on your custom order
cluster_avg_scores <- cluster_avg_scores[match(custom_order, rownames(cluster_avg_scores)), , drop = FALSE]

# Generate heatmap with manual ordering
pheatmap(cluster_avg_scores, cluster_rows = FALSE, cluster_cols = FALSE, 
         color = colorRampPalette(c("blue", "white", "red"))(100))



























