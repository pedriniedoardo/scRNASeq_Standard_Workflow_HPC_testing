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

qc_features      = snakemake@params$qc_features
figure_extension = snakemake@params$figure_extension
out_location     = snakemake@params$out_location
resolutions      = snakemake@params$resolutions

# ======================================================================
# == load merged and integrated object ==
# ======================================================================

merged_rds     = snakemake@input$merged_rds
integrated_rds = snakemake@input$integrated_rds

merged         = readRDS(merged_rds)
integrated     = readRDS(integrated_rds)
n_samples      = length(unique(integrated$orig.ident))
n_res          = length(resolutions)

# ======================================================================
# == set panel rows and cols ===
# ======================================================================

panel_layout <- function(n) {
  ncol <- ceiling(sqrt(n))
  nrow <- ceiling(n / ncol)
  return(c(nrow, ncol))
}

n_rows_clusters = panel_layout(n_res)[1]
n_cols_clusters = panel_layout(n_res)[2]

n_rows_samples  = panel_layout(n_samples)[1]
n_cols_samples  = panel_layout(n_samples)[2]

# ======================================================================
# == QC violin plots ==
# ======================================================================

for(feat in qc_features) {
  for(ext in figure_extension) {
    ggsave(
      paste0(out_location, "Seurat/plot/vln_", feat, ".", ext),
      VlnPlot(integrated, features = feat, group.by = 'orig.ident'), 
      height = 5, width = 3 + (1.5 * n_samples))
  }
}

# ======================================================================
# == QC dim plots ==
# ======================================================================

for(feat in qc_features) {
  for(ext in figure_extension) {
    ggsave(
      paste0(out_location, "Seurat/plot/dim_", feat, ".", ext),
      FeaturePlot(integrated, features = feat), 
      height = 6, width = 8
      )
  }
}

# ======================================================================
# == idents dim plots ==
# ======================================================================

# merged and integrated

for(object_type in c('integrated','merged')) {
  object = switch(object_type, 
                  'integrated' = integrated, 
                  'merged'     = merged
                  )
  for(ext in figure_extension) {
    ggsave(
      paste0(out_location, "Seurat/plot/dim_", object_type, "_ident.", ext),
      DimPlot(object, group.by = 'orig.ident'),
      height = 6, width = 8
    )
  }
}

# ======================================================================
# == clusters dim plots ==
# ======================================================================

# merged and integrated

for(object_type in c('integrated','merged')) {
  object = switch(object_type, 
                  'integrated' = integrated, 
                  'merged'     = merged
  )
  for(ext in figure_extension) {
    gg_list = lapply(resolutions, function(res) {
      DimPlot(object, group.by = paste0('RNA_snn_res.', res))
      })
    ggsave(
      paste0(out_location, "Seurat/plot/dim_", object_type, "_cluster.", ext),
      wrap_plots(gg_list),
      height = 4 + (2.5 * n_rows_clusters), width = 4 + (2.5 * n_cols_clusters)
      )
  }
}


