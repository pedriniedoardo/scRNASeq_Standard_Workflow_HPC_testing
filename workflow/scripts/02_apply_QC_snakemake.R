# ======================================================================
# == load libraries ==
# ======================================================================

library(scater)
library(Seurat)
library(tidyverse)
library(robustbase)
library(patchwork)
library(scDblFinder)

# ======================================================================
# == custom functions ==
# ======================================================================

# specify the version of Seurat Assay -------------------------------------
# set seurat compatible with seurat4 workflow
options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 1000 * 1024^2)

# ======================================================================
# == Snakemake integation ==
# ======================================================================

# read in the rds file
rds <- snakemake@input$rds

# sample id
id_sample <- snakemake@wildcards$sample_name

# read in the LUT_QC
LUT_QC <- snakemake@input$LUT_QC

message("input rds: ", rds)
message("sample id: ", id_sample)
message("table LUT QC: ", LUT_QC)

# define the output
# out_id_object <-  "results/Seurat/object/01_connect_5k_pbmc_NGSC3_ch1_gex_1_obj_preQC.rds"
# out_id_meta <- "results/Seurat/table/01_connect_5k_pbmc_NGSC3_ch1_gex_1_meta_preQC.rds"
out_id_object <- snakemake@output$rds
out_id_meta <- snakemake@output$meta

message("output rds: ", out_id_object)
message("output meta: ", out_id_meta)


# ======================================================================
# == load input ==
# ======================================================================

df_LUT <- read_csv(LUT_QC) %>%
  filter(sample_name == id_sample)

# 1000
featureLow_thr <- df_LUT$featureLow_thr

# 6000
featureHigh_thr <- df_LUT$featureHigh_thr

# 15
mito_thr <- df_LUT$mito_thr

# "01000_06000_15_V5"
# build the label
label <- paste(featureLow_thr,featureHigh_thr,mito_thr,"V5",sep = "_")

# ======================================================================
# == run the standard processing ==
# ======================================================================

scobj <- readRDS(rds)

# add the filtering label
scobj$label <- label

# add the filtering variable based on the fixed threshold
scobj$discard_threshold <- scobj@meta.data %>%
  mutate(test = percent.mt > mito_thr | nFeature_RNA < featureLow_thr | nFeature_RNA > featureHigh_thr) %>% 
  pull(test)

# preprocess the dataset before the doublet identification as recommended in:
# https://bioconductor.org/packages/release/bioc/vignettes/scDblFinder/inst/doc/scDblFinder.html
# 1.5.11
# according to the documentation the doublet identification should be run before the fine filtering
# remove low coverage cells and proprocess to generata clusters needed for the doublet identification
scobj <- subset(scobj,subset = nCount_RNA > 500) %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters() %>%
  # do I run the UMAP ? I do not need it for the doublet identification, but can be useful in case someone wants to explore an individual sample
  RunUMAP(dims = 1:30)

# run scDblFinder after filtering the low coverage cells
sce_scobj <- scDblFinder(GetAssayData(scobj, layer="counts"), clusters=Idents(scobj))

# port the resulting scores back to the Seurat object:
scobj$scDblFinder.score <- sce_scobj$scDblFinder.score
scobj$scDblFinder.class <- sce_scobj$scDblFinder.class

# perform the filtering based on the fixed threshold defined
scobj_filter <- subset(scobj, subset = discard_threshold == 0)

# preprocess data after filtering 
scobj_filter <- scobj_filter %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters() %>%
  RunUMAP(dims = 1:30) 

# ======================================================================
# == save output ==
# ======================================================================

saveRDS(scobj_filter,out_id_object)
write_tsv(scobj_filter@meta.data %>% rownames_to_column("barcodes"),out_id_meta)