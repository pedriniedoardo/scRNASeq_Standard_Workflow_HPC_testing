# ======================================================================
# == load libraries ==
# ======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

# ======================================================================
# == custom functions ==
# ======================================================================

# specify the version of Seurat Assay -------------------------------------
# set seurat compatible with seurat4 workflow
options(Seurat.object.assay.version = "v5")

# ======================================================================
# == parameters ==
# ======================================================================

# Seurat:
#   embedding:
#   n_dims: 30
# resolutions: [.3, .6, 1, 1.2]
# markers: 
#   min_pct: .1
# 
# report:
#   figure_extension: ['pdf', 'png']
# qc_features: ['nFeature_RNA','nCount_RNA','percent.mt','percent.ribo','percent.globin','scDblFinder.score']

# config["report"]["qc_features"]
# qc_features <- snakemake@params$qc_features
qc_features <- c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo","percent.globin","scDblFinder.score") 

# config["report"]["figure_extension"]
# figure_extension <- snakemake@params$figure_extension
figure_extension <- c("png","pdf")

# config["out_location"]
# out_location <- snakemake@params$out_location
out_location <- "out/plot/"

# config["Seurat"]["embedding"]["resolutions"]
# resolutions <- snakemake@params$resolutions
resolutions <- c(0.3, 0.6, 1, 1.2)

# cellbender_tag <- snakemake@params$cellbender_tag
cellbender_tag <- "default"

# ======================================================================
# == load merged and integrated object ==
# ======================================================================

# merged_rds <- snakemake@input$merged_rds
# integrated_rds <- snakemake@input$integrated_rds

# --- inputs ---
merged_rds     <- "../../../results/Seurat/object/merged_obj_default.rds"
integrated_rds <- "../../../results/Seurat/object/integrated_obj_default.rds"

merged <- readRDS(merged_rds)
integrated <- readRDS(integrated_rds)

n_samples <- length(unique(integrated$orig.ident))
n_res <- length(resolutions)

# ======================================================================
# == set panel rows and cols ===
# ======================================================================

panel_layout <- function(n) {
  ncol <- ceiling(sqrt(n))
  nrow <- ceiling(n / ncol)
  return(c(nrow, ncol))
}

n_rows_clusters <- panel_layout(n_res)[1]
n_cols_clusters <- panel_layout(n_res)[2]

n_rows_samples  <- panel_layout(n_samples)[1]
n_cols_samples  <- panel_layout(n_samples)[2]

# ======================================================================
# == QC violin plots ==
# ======================================================================

for (feat in qc_features) {
  for (ext in figure_extension) {
    ggsave(
      paste0(out_location, "Seurat/plot/vln_", feat, "_", cellbender_tag, ".", ext),
      VlnPlot(integrated, features = feat, group.by = 'orig.ident'),
      height = 7, width = 5 + (1.5 * n_samples))
  }
}

# ======================================================================
# == QC dim plots ==
# ======================================================================

for (feat in qc_features) {
  for (ext in figure_extension) {
    ggsave(
      paste0(out_location, "Seurat/plot/dim_", feat, "_", cellbender_tag, ".", ext),
      FeaturePlot(integrated, features = feat,order = T),
      height = 6, width = 8
    )
  }
}

# ======================================================================
# == idents dim plots ==
# ======================================================================

for (object_type in c('integrated', 'merged')) {
  object <- switch(object_type,
                   'integrated' = integrated,
                   'merged'     = merged)
  for (ext in figure_extension) {
    ggsave(
      paste0(out_location, "Seurat/plot/dim_", object_type, "_ident_", cellbender_tag, ".", ext),
      DimPlot(object, group.by = 'orig.ident',order = T),
      height = 6, width = 8
    )
  }
}

# ======================================================================
# == clusters dim plots ==
# ======================================================================

for (object_type in c('integrated', 'merged')) {
  object <- switch(object_type,
                   'integrated' = integrated,
                   'merged'     = merged)
  for (ext in figure_extension) {
    gg_list <- lapply(resolutions, function(res) {
      DimPlot(object, group.by = paste0('RNA_snn_res.', res))
    })
    ggsave(
      paste0(out_location, "Seurat/plot/dim_", object_type, "_cluster_", cellbender_tag, ".", ext),
      wrap_plots(gg_list),
      height = 4 + (2.5 * n_rows_clusters), width = 4 + (2.5 * n_cols_clusters)
    )
  }
}
