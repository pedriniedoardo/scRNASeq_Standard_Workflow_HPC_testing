# ======================================================================
# == load libraries ==
# ======================================================================

# add the installtion step of the GitHub package.
# all the dependencies have been tracked in the yaml file

# Bootstrap CyteTypeR
if (!requireNamespace("CyteTypeR", quietly = TRUE)) {
  message("Installing CyteTypeR from GitHub...")
  # fix the version
  # remotes::install_github("NygenAnalytics/CyteTypeR", ref = "v0.4.2", upgrade = "never")
  
  # lates version
  remotes::install_github("NygenAnalytics/CyteTypeR", upgrade = "never")
}

library(CyteTypeR)
library(tidyverse)
library(patchwork)
library(Matrix)
library(Seurat)
library(presto)
library(duckdb)

# ======================================================================
# == Snakemake integation ==
# ======================================================================

# input:
rds <- snakemake@input$rds
markers <- snakemake@input$markers
metadata <- snakemake@input$metadata
cellbender_tag <- snakemake@params$cellbender_tag

# rds <- "../../../results/Seurat/object/integrated_obj_cellbender.rds"
# markers <- "../../../results/Seurat/table/integrated_markers_cellbender.tsv"
# metadata <- "../../../config/meta_annotation.json"

message("input rds: ", rds)
message("input metadata: ", metadata)
message("cellbender tag: ", cellbender_tag)

# output:
out_id_rds <- snakemake@output$rds
out_id_varsh5 <- snakemake@output$vars_h5_path
out_id_duckdb <- snakemake@output$obs_duckdb_path
out_id_filename <- snakemake@output$query_filename

# out_id_rds <- "out/object/integrated_obj_cellbender_CyteTypeR.rds"
# out_id_varsh5 <- "out/object/vars.h5"
# out_id_duckdb <- "out/object/obs.duckdb"
# out_id_filename <- "out/object/query.json"

message("output CyteTypeR rds: ", out_id_rds)
message("output CyteTypeR vars.h5: ", out_id_varsh5)
message("output CyteTypeR obs.duckdb: ", out_id_duckdb)
message("output CyteTypeR query.json: ", out_id_filename)

# params
# token <- snakemake@params$token
token <- Sys.getenv("CYTETYPER_TOKEN")
# this should be handled in the snakemake rule
if (token == "") {
  stop("CYTETYPER_TOKEN environment variable is not set or is empty. Please export it before running Snakemake.")
}

cov_markers <- snakemake@params$cov_markers

# token <- ""
# cov_markers <- "RNA_snn_res.0.3"

# ======================================================================
# == load input ==
# ======================================================================

# read in the full object
scobj <- readRDS(rds)

# read in the makrers
df_markers <- read_tsv(markers)

# read in the metadata
list_metadata <- yaml::read_yaml(metadata)

# modify the title to incorporate the cellbender tag
list_metadata$title <- paste(list_metadata$title,cellbender_tag,sep = "_")

# ======================================================================
# == run the standard processing ==
# ======================================================================

# reshape the marker output to be accepted
# This one assumes that the order of the genes is maintained from the output of Seurat's FindAllMarkers() function during the integration step
df_markers_fix <- df_markers %>%
  # mutate(cluster = as.numeric(cluster)) %>%
  filter(resolution == str_remove(cov_markers,pattern = "RNA_snn_res.")) %>%
  select(cluster,gene,avg_log2FC)

# run the query preparation
prepped_data <- PrepareCyteTypeR(obj = scobj,
                                 marker_table = df_markers_fix,
                                 n_top_genes = 10,
                                 group_key = cov_markers,
                                 aggregate_metadata = TRUE,
                                 coordinates_key = "umap",
                                 vars_h5_path = out_id_varsh5,
                                 obs_duckdb_path = out_id_duckdb)

# run the annotation
results <- CyteTypeR(obj = scobj,
                     auth_token = token,
                     prepped_data = prepped_data, 
                     study_context = list_metadata$study_context, 
                     metadata = list_metadata,
                     query_filename = out_id_filename)

# ======================================================================
# == save output ==
# ======================================================================

saveRDS(results,out_id_rds)
