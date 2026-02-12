# ======================================================================
# == load libraries ==
# ======================================================================

library(scater)
library(Seurat)
library(tidyverse)
library(robustbase)
library(patchwork)
library(scuttle)

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
result_folder <- "/beegfs/scratch/ric.cosr/pedrini.edoardo/test/test_seurat_COSR/data/connect_5k_pbmc_NGSC3_ch1_gex_1"
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
data <- Read10X(data.dir = file.path(result_folder,in_id_data))

# ======================================================================
# == run the standard processing ==
# ======================================================================

# crete the object
datasc <- CreateSeuratObject(counts = data, project = in_id_sample,
                             # potentially parametrize the values below
                             # do not filter any gene before the integration
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