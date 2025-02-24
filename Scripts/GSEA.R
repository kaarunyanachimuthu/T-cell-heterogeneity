library(Seurat)
library(fgsea)
library(msigdbr)
library(dplyr)
library(tibble)

# Step 1: Perform DE analysis between conditions
Cd4 <- readRDS('./Cd4_pbs_vs_cue.RDS')
six <- subset(Cd4, idents = 6)
six$seurat_clusters.treatment<-paste(Idents(six), six$treatment, sep="_")
six$seurat_clusters<- Idents(six)
Idents(six)<-"seurat_clusters.treatment"

de_results <- FindMarkers(six, ident.1 = "6_PBS", ident.2 = "6_CUE", min.pct = 0.25, only.pos = TRUE)
ranked_genes <- de_results %>%
  rownames_to_column(var = "gene") %>%
  arrange(desc(avg_log2FC)) %>%
  pull(avg_log2FC, name = "gene")

# Step 2: Get Hallmark gene sets from MSigDB
msigdb_sets <- msigdbr(species = "Mus musculus", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)
pathways <- split(msigdb_sets$gene_symbol, msigdb_sets$gs_name)

# Step 3: Run fgsea
fgsea_results <- fgsea(pathways = pathways, 
                       stats = ranked_genes, 
                       minSize = 15, 
                       maxSize = 500, 
                       nperm = 1000)

# Step 4: Sort by NES and view top results
fgsea_results <- fgsea_results %>% arrange(desc(NES))
head(fgsea_results)

# Step 5: Plot the top enriched pathways
topPathways <- fgsea_results %>% filter(padj < 1) %>% head(10) %>% pull(pathway)
plotGseaTable(pathways = pathways[topPathways],
              stats = ranked_genes,
              fgseaRes = fgsea_results,
              gseaParam = 0.5)
