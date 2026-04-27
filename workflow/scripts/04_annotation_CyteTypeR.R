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
# results <- CyteTypeR(obj = scobj,
#                      auth_token = token,
#                      prepped_data = prepped_data, 
#                      study_context = list_metadata$study_context, 
#                      metadata = list_metadata,
#                      query_filename = out_id_filename)

# --- new implementation to recover output data ---
results  <- NULL

# handle the console output to collect the job_id generated
msg_log <- tempfile()
msg_con <- file(msg_log, open = "wt")
sink(msg_con, type = "message")

# run the CyteTypeR
results <- tryCatch(
  {
    message("Submitting query to CyteTypeR server...")
    CyteTypeR(obj = scobj,
              auth_token = token,
              prepped_data = prepped_data,
              study_context = list_metadata$study_context,
              metadata = list_metadata,
              query_filename = out_id_filename)
  },
  error = function(e) {
    message("CyteTypeR failed: ", conditionMessage(e))
    NULL
  },
  finally = {
    sink(type = "message")
    close(msg_con)
  }
)

# replay captured log so it is printed in the console
logged <- readLines(msg_log)
cat(paste(logged, collapse = "\n"), "\n")

# extract job_id from captured output
all_log <- paste(logged, collapse = " ")
job_id <- regmatches(
  all_log,
  regexpr("[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}", all_log)
)
if (length(job_id) == 0) job_id <- NULL

# if the run has failed because of a dropped connection, pull the meta_resutls
meta_results <- NULL
if (is.null(results)) {
  if (is.null(job_id)) {
    stop("CyteTypeR failed before a job_id was obtained. Cannot retrieve results.")
  }
  # Reconstruct job details in scobj@misc so GetResults uses the same pipeline as CyteTypeR (normalize + correct endpoint), instead of the bare job_id path.
  scobj@misc[["cytetype_jobDetails"]] <- list(
    job_id = job_id,
    api_url = "https://prod.cytetype.nygen.io",
    group_key = cov_markers,
    cluster_labels = prepped_data$clusterLabels
  )
  message("Retrieving results via GetResults (job_id: ", job_id, ")...")
  meta_results <- GetResults(obj = scobj,
                             auth_token = token,
                             results_prefix = "cytetype")
}

# testing
# test <- GetResults(obj = scobj,
#                    auth_token = token,
#                    results_prefix = "cytetype")
#
# test2 <- GetResults(auth_token = token,
#                     results_prefix = "cytetype",
#                     job_id = job_id)
# test3 <- CyteTypeR:::.transform_results_seurat(
#   meta_results,
#   cluster_map = prepped_data$clusterLabels)

# GetResults returns the normalized result list, not an annotated Seurat object.
# Transform and apply annotations to scobj to mirror what CyteTypeR does internally. cluster_map handles the sequential ID -> original label translation internally.
if (!is.null(meta_results)) {
  transformed_results <- CyteTypeR:::.transform_results_seurat(
    meta_results,
    cluster_map = prepped_data$clusterLabels
  )

  # add the main annotations from the collected results
  ann_colname <- paste("cytetype", "annotation",cov_markers, sep = "_")
  onto_colname <- paste("cytetype", "cellOntologyTerm", cov_markers, sep = "_")
  ontoID_colname <- paste("cytetype", "cellOntologyTermID", cov_markers, sep = "_")
  state_colname <- paste("cytetype", "cellState", cov_markers, sep = "_")
  # this is missing from the regular run of CyteTypeR output
  granular_colname <- paste("cytetype", "granularAnnotation", cov_markers, sep = "_")

  cells_group <- as.character(scobj@meta.data[[cov_markers]])
  scobj@meta.data[[ann_colname]] <- factor(setNames(transformed_results$annotation, transformed_results$clusterId)[cells_group])
  scobj@meta.data[[onto_colname]] <- factor(setNames(transformed_results$ontologyTerm, transformed_results$clusterId)[cells_group])
  scobj@meta.data[[ontoID_colname]] <- factor(setNames(transformed_results$ontologyTermID, transformed_results$clusterId)[cells_group])
  scobj@meta.data[[state_colname]] <- factor(setNames(transformed_results$cellState, transformed_results$clusterId)[cells_group])
  # this is missing from the regular run of CyteTypeR output
  scobj@meta.data[[granular_colname]] <- factor(setNames(transformed_results$granularAnnotation, transformed_results$clusterId)[cells_group])

  # for some reason in the default run of CyteTypeR there are also some other covariates added to the metadata but are redundant
  # cytetype_RNA_snn_res.x -> annoation
  # cytetype_ontologyTerm_RNA_snn_res.x -> cellOntologyTerm
  # cytetype_ontologyID_RNA_snn_res.x -> cellOntologyTermID

  # I have deceded to skip them

  scobj@misc$cytetype_results <- transformed_results

  results <- scobj
}

# table(results@meta.data$RNA_snn_res.0.5,results@meta.data$cytetype_cellState_RNA_snn_res.0.5)

# ======================================================================
# == save output ==
# ======================================================================

saveRDS(results,out_id_rds)
