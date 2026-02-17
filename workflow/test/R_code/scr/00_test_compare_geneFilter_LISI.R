# libraries 
library(Seurat)
library(cluster)
library(lisi)
library(ggplot2)
library(dplyr)
library(patchwork)

# 1. Load data
test_filter <- readRDS("out/object/pbmc_test/connect_5k_pbmc_harmony_filter.rds")
test_all <- readRDS("out/object/pbmc_test/connect_5k_pbmc_harmony_all.rds")

# 2. Silhouette Analysis (Cluster Quality)
# Measures how similar a cell is to its own cluster vs other clusters.
calc_silhouette <- function(obj, reduction = "harmony", res = "RNA_snn_res.0.3") {
  # Calculate distance matrix in Harmony space
  dist_matrix <- dist(Embeddings(obj, reduction))
  clusters <- as.numeric(obj@meta.data[[res]])
  
  # Calculate silhouette
  sil <- silhouette(clusters, dist_matrix)
  return(data.frame(
    avg_sil = mean(sil[, 3]),
    sil_widths = sil[, 3],
    cluster = as.factor(clusters)
  ))
}

sil_all <- calc_silhouette(test_all) %>% mutate(test = "All Genes")
sil_filt <- calc_silhouette(test_filter) %>% mutate(test = "Filtered")

# 3. LISI Analysis (Batch Mixing vs Cell Type Separation)
# iLISI: Integration LISI (Batch mixing). Closer to # of batches is better.
# cLISI: Cell-type LISI. Closer to 1 is better (less mixing of different types).
calc_lisi_metrics <- function(obj, reduction = "harmony") {
  coords <- Embeddings(obj, reduction)
  meta <- obj@meta.data[, c("orig.ident", "RNA_snn_res.0.3")]
  
  res_lisi <- compute_lisi(coords, meta, c("orig.ident", "RNA_snn_res.0.3"))
  colnames(res_lisi) <- c("iLISI", "cLISI")
  return(res_lisi)
}

lisi_all <- calc_lisi_metrics(test_all) %>% mutate(test = "All Genes")
lisi_filt <- calc_lisi_metrics(test_filter) %>% mutate(test = "Filtered")

# 4. Visualization & Reporting
comparison_df <- bind_rows(
  sil_all %>% select(val = avg_sil, test) %>% mutate(metric = "Avg Silhouette"),
  sil_filt %>% select(val = avg_sil, test) %>% mutate(metric = "Avg Silhouette"),
  data.frame(val = mean(lisi_all$iLISI), test = "All Genes", metric = "iLISI (Batch Mixing)"),
  data.frame(val = mean(lisi_filt$iLISI), test = "Filtered", metric = "iLISI (Batch Mixing)"),
  data.frame(val = mean(lisi_all$cLISI), test = "All Genes", metric = "cLISI (Type Separation)"),
  data.frame(val = mean(lisi_filt$cLISI), test = "Filtered", metric = "cLISI (Type Separation)")
)

# Plot Results
p_sil <- ggplot(bind_rows(sil_all, sil_filt), aes(x=test, y=sil_widths, fill=test)) +
  geom_violin(alpha=0.7) + geom_boxplot(width=0.1, outlier.shape=NA) +
  theme_minimal() + ggtitle("Cluster Stability (Silhouette)")

p_lisi <- ggplot(bind_rows(lisi_all, lisi_filt), aes(x=test, y=iLISI, fill=test)) +
  geom_boxplot() + theme_minimal() + 
  ggtitle("Batch Mixing (iLISI)") + 
  labs(subtitle = paste0("Ideal iLISI: ", length(unique(test_all$orig.ident))))

p_sil / p_lisi