# -------------------------------------------------------------------------
test_cellbender <- Read_CellBender_h5_Mat(file_name = "../../../results/cellbender/merged/connect_5k_pbmc_NGSC3_ch1_gex_1/cellbender_out_filtered.h5")
test_default <- Read_CellBender_h5_Mat(file_name = "../../../results/cellranger/merged/connect_5k_pbmc_NGSC3_ch1_gex_1/outs/filtered_feature_bc_matrix.h5")

dim(test_cellbender)
dim(test_default)

# -------------------------------------------------------------------------
test_cellbender <- readRDS("../../../results/Seurat/object/integrated_obj_cellbender.rds")
test_default <- readRDS("../../../results/Seurat/object/integrated_obj_default.rds")

DimPlot(test_cellbender)
DimPlot(test_default)

dim(test_cellbender)
dim(test_default)

test_cellbender@meta.data
test_default@meta.data


# -------------------------------------------------------------------------
test <- readRDS("../../../../sc-rna-seq-cosr-standard-workflow-test/results/Seurat/object/integrated_obj.rds")

test2 <- readRDS("/beegfs/scratch/ric.cosr/pedrini.edoardo/project_edoardo/240531_scRNAseq_MSSpinal_Absinta/out/object/analysis_R44/26_sobj_integrated_cleanup_manualAnnotation.rds")
