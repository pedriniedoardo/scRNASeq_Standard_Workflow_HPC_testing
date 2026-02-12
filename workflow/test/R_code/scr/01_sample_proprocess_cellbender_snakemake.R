# ======================================================================
# == load libraries ==
# ======================================================================

library(scater)
library(Seurat)
library(tidyverse)
library(robustbase)
library(patchwork)
library(scuttle)
library(hdf5r)
library(scCustomize)

# ======================================================================
# == custom functions ==
# ======================================================================

# specify the version of Seurat Assay -------------------------------------
# set seurat compatible with seurat4 workflow
options(Seurat.object.assay.version = "v5")

# ======================================================================
# == Snakemake integation ==
# ======================================================================

# # define the project folder
# # folder = directory(config["out_location"] + "cellranger/merged/{sample_name}")
# result_folder <- "/idle/ric.cosr/ric.cosr/pedrini.edoardo/sc-rna-seq-cosr-standard-workflow-test2/workflow/test/cellbender/test_connect_5k_pbmc"
# # result_folder <- snakemake@input$result_folder
# 
# # define the data to be loaded
# # in_id_data <- "cellbender_output_filtered.h5"
# in_id_data <- "filtered_feature_bc_matrix.h5"
# 
# # define ogganisim
# in_id_org <- "9606"
# # in_id_org <- snakemake@params$id_org

# -------------------------------------------------------------------------
# new implementation after introducing the cellbender module
# Get the input object (this is either the default h5 from cellranger or the one from cellbender)
# input_target <- snakemake@input[["target"]]
input_target <- "../../../results/cellbender/merged/connect_5k_pbmc_NGSC3_ch1_gex_1/cellbender_out_filtered.h5"
input_target <- "../../../results/cellranger/merged/connect_5k_pbmc_NGSC3_ch1_gex_1"

message("input target: ", input_target)

# Get the sample id
# folder = directory(config["out_location"] + "cellranger/merged/{sample_name}")
# in_id_sample <- snakemake@wildcards$sample_name
in_id_sample <- "connect_5k_pbmc_NGSC3_ch1_gex_1"

message("sample id: ", in_id_sample)

# Get the boolean flag
# use_cellbender <- snakemake@params[["use_cellbender"]]
use_cellbender <- T

# define ogganisim
in_id_org <- "9606"
# in_id_org <- snakemake@params$id_org
message("id organism: ", in_id_org)

# 3. Determine the final h5 path logic
if (use_cellbender) {
  # If CellBender ran, the input_target IS the h5 file
  h5_file <- input_target
  message(paste("Loading CellBender output from:", h5_file))
  
} else {
  # If CellRanger ran, the input_target IS the folder
  # We must construct the path to the file inside it
  h5_file <- file.path(input_target, "outs", "filtered_feature_bc_matrix.h5")
  message(paste("Loading CellRanger output from:", h5_file))
}

# 4. Check if file exists (good practice for debugging)
if (!file.exists(h5_file)) {
  stop(paste("The h5 file was not found at:", h5_file))
}

# -------------------------------------------------------------------------

# define the output
# out_id_object <-  "results/Seurat/object/01_connect_5k_pbmc_NGSC3_ch1_gex_1_obj_preQC.rds"
# out_id_meta <- "results/Seurat/table/01_connect_5k_pbmc_NGSC3_ch1_gex_1_meta_preQC.rds"
out_id_object <-  "out/object/01_cellbender_test_connect_5k_pbmc_obj_preQC.rds"
out_id_meta <- "out/table/01_cellbender_test_connect_5k_pbmc_meta_preQC.tsv"
# out_id_object <- snakemake@output$rds
# out_id_meta <- snakemake@output$meta

message("output rds: ", out_id_object)
message("output meta: ", out_id_meta)

# ======================================================================
# == load input ==
# ======================================================================

# read in the matrix
# change the implementaiton by reading in the h5 rather than the folder
# data <- Read10X_h5(filename = file.path(result_folder,in_id_data))
# this function is needed as the output of cellbender cannot be read directly from Read10X_h5
data <- Read_CellBender_h5_Mat(file_name = h5_file)

# ======================================================================
# == run the standard processing ==
# ======================================================================

# crete the object
datasc <- CreateSeuratObject(counts = data, project = in_id_sample,
                             # potentially parametrize the values below
                             # do not filter any gene bofore the integration
                             min.cells = 0,
                             min.features = 200)

# add the mefadata accordin to the organism
if(in_id_org == "9606"){
  # add the metadata
  datasc$percent.mt <- PercentageFeatureSet(datasc, pattern = "^MT-")
  datasc$percent.ribo <- PercentageFeatureSet(datasc, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
  datasc$percent.globin <- Seurat::PercentageFeatureSet(datasc,pattern = "^HB[^(P)]")
} else if(in_id_org == "10090"){
  # add the metadata
  datasc$percent.mt <- PercentageFeatureSet(datasc, pattern = "^mt-")
  datasc$percent.ribo <- PercentageFeatureSet(datasc, pattern = "^Rp[sl][[:digit:]]|^Rplp[[:digit:]]|^Rpsa")
  datasc$percent.globin <- Seurat::PercentageFeatureSet(datasc,pattern = "^Hb[^(p)]")
} else {
  # decide whether to run something or nothing
}

# -------------------------------------------------------------------------
# add the filtering variable based on the adaptive threshold multivalue
stats <- cbind(log10(datasc@meta.data$nCount_RNA),
               log10(datasc@meta.data$nFeature_RNA),
               datasc@meta.data$percent.mt)

outlying <- adjOutlyingness(stats, only.outlyingness = TRUE)
multi.outlier <- isOutlier(outlying, type = "higher")
datasc$discard_multi <- as.vector(multi.outlier)

# add the filtering variable based on the adaptive threshold single values
high_QC_mito <- isOutlier(datasc@meta.data$percent.mt, type="high", log=TRUE)
QC_features <- isOutlier(datasc@meta.data$nFeature_RNA, type="both", log=TRUE)

datasc$discard_single <- high_QC_mito | QC_features
# -------------------------------------------------------------------------

# ======================================================================
# == save output ==
# ======================================================================

saveRDS(datasc,out_id_object)
write_tsv(datasc@meta.data %>% rownames_to_column("barcodes"),out_id_meta)
