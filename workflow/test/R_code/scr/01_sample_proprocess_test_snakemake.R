# ======================================================================
# == load libraries ==
# ======================================================================

library(scater)
library(Seurat)
library(tidyverse)
library(robustbase)
library(patchwork)
library(scuttle)
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

# define the project folder
# folder = directory(config["out_location"] + "cellranger/merged/{sample_name}")
# result_folder <- "/beegfs/scratch/ric.cosr/pedrini.edoardo/test/test_seurat_COSR/data/connect_5k_pbmc_NGSC3_ch1_gex_1"
result_folder <-"../../../results/cellranger/merged/connect_5k_pbmc_NGSC3_ch1_gex_1"
# result_folder <- snakemake@input$result_folder

# define the sample id
# folder = directory(config["out_location"] + "cellranger/merged/{sample_name}")
in_id_sample <- "connect_5k_pbmc_NGSC3_ch1_gex_1"
# in_id_sample <- snakemake@wildcards$sample_name

# define the data to be loaded
in_id_data <- "outs/filtered_feature_bc_matrix"

# define ogganisim
in_id_org <- "9606"
# in_id_org <- snakemake@params$id_org

message("input folder: ", result_folder)
message("sample id: ", in_id_sample)
message("id organism: ", in_id_org)

# define the output
# out_id_object <-  "results/Seurat/object/01_connect_5k_pbmc_NGSC3_ch1_gex_1_obj_preQC.rds"
# out_id_meta <- "results/Seurat/table/01_connect_5k_pbmc_NGSC3_ch1_gex_1_meta_preQC.rds"
out_id_object <-  "out/object/01_connect_5k_pbmc_NGSC3_ch1_gex_1_obj_preQC.rds"
out_id_meta <- "out/table/01_connect_5k_pbmc_NGSC3_ch1_gex_1_meta_preQC.tsv"
# out_id_object <- snakemake@output$rds
# out_id_meta <- snakemake@output$meta

message("output rds: ", out_id_object)
message("output meta: ", out_id_meta)

# ======================================================================
# == load input ==
# ======================================================================

# read in the matrix
# change the implementaiton by reading in the h5 rather than the folder
data_folder <- Read10X(data.dir = file.path(result_folder,in_id_data))

# is id different from loading the h5 file?
data_h5 <- Read10X_h5(filename = file.path(result_folder,"outs/filtered_feature_bc_matrix.h5"))

# is id different using the funciton needed for the cellbender output?
data_h5_test <- Read_CellBender_h5_Mat(file_name = file.path(result_folder,"outs/filtered_feature_bc_matrix.h5"))

dim(data_folder)
dim(data_h5)
dim(data_h5_test)

# are the three version the same?
all.equal(data_folder,data_h5)
all.equal(data_folder,data_h5_test)
