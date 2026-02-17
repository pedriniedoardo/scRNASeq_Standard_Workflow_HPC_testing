# AIM ---------------------------------------------------------------------
# test the integration approach without any filtering for the genes

# renv integration --------------------------------------------------------
# to load the packages
# source(".Rprofile")

# libraries ---------------------------------------------------------------
library(scater)
library(Seurat)
library(tidyverse)
library(robustbase)
library(patchwork)
library(scDblFinder)
library(presto)
library(harmony)

# Snakemake integation ----------------------------------------------------
# define the inputs

# define the input spatial file
# rds <- c("out/object/test_human_default/connect_5k_pbmc_NGSC3_ch1_gex_listTest.rds")
rds <- snakemake@input$spe_input
message("input spe object: ", rds)

# define the output file
# out_object <- "out/object/test_human_default/connect_5k_pbmc_NGSC3_ch1_gex_all.rds"
out_object <- snakemake@output$spe_output
message("output object: ", out_object)

# specify the version of Seurat Assay -------------------------------------
# set seurat compatible with seurat4 workflow
options(Seurat.object.assay.version = "v5")

# parameters --------------------------------------------------------------
n_dims = 30
resolutions = c(.3, .6, 1, 1.2)
min_pct = 0.1

# load the data -----------------------------------------------------------
# load objects 
all_objects <- readRDS(rds)

# rename cells (to avoid duplicated barcodes)
# fix the metadata to avoid confusion with resolutions names
lapply(all_objects, function(obj) {
  obj = RenameCells(obj, add.cell.id = Project(obj))
  
  # remove the clustering information to avoid confusions after integration
  meta_in <- obj@meta.data
  meta_fix <- meta_in %>%
    select(-contains("RNA_snn_res")) %>%
    select(-"seurat_clusters")
  
  # swap the metadata
  obj@meta.data <- meta_fix
  
  # remove the scaled.data layer
  obj[["RNA"]]$scale.data <- NULL
  
  return(obj)
}) -> all_objects


# merge and processing ----------------------------------------------------
n_cells <- sum(sapply(all_objects,ncol))
options(future.globals.maxSize = length(all_objects) * n_cells*.1 * 1024^2)

# merge
merged_full <- merge(all_objects[[1]], all_objects[-1])

# -------------------------------------------------------------------------
# apply the gene filtering after merging
# manual implementation
counts <- GetAssayData(JoinLayers(merged_full[["RNA"]]), assay = "RNA", layer = "counts")
genes_to_keep <- rowSums(counts > 0) >= 20
merged <- merged_full[genes_to_keep, ]
# -------------------------------------------------------------------------

# standard pre-processing
# (-defaults- arguments exposed)
merged <- merged %>% 
  NormalizeData(
    assay = NULL,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    margin = 1,
    verbose = TRUE
  ) %>% 
  FindVariableFeatures(
    assay = NULL,
    selection.method = "vst",
    loess.span = 0.3,
    clip.max = "auto",
    # mean.function = FastExpMean,
    # dispersion.function = FastLogVMR,
    num.bin = 20,
    binning.method = "equal_width",
    nfeatures = 2000,
    mean.cutoff = c(0.1, 8),
    dispersion.cutoff = c(1, Inf),
    verbose = TRUE
  ) %>%
  ScaleData(
    features = NULL, # Default is variable features.
    assay = NULL,
    vars.to.regress = NULL,
    split.by = NULL,
    model.use = "linear",
    use.umi = FALSE,
    do.scale = TRUE,
    do.center = TRUE,
    scale.max = 10,
    block.size = 1000,
    min.cells.to.block = 3000,
    verbose = TRUE
  ) %>%
  RunPCA(
    assay = NULL,
    features = NULL, # If features=NULL, PCA will be run using the variable features for the Assay
    npcs = 50,
    rev.pca = FALSE,
    weight.by.var = TRUE,
    verbose = TRUE,
    ndims.print = 1:5,
    nfeatures.print = 30,
    reduction.name = "pca",
    reduction.key = "PC_",
    seed.use = 42
  )

# integration -------------------------------------------------------------
integrated = IntegrateLayers(
  merged, 
  method = HarmonyIntegration,
  orig.reduction = "pca",
  assay = NULL,
  features = NULL, # Ignored
  scale.layer = "scale.data", # Ignored
  new.reduction = "harmony",
  layers = NULL, # Ignored
  npcs = 50L,
  key = "harmony_",
  theta = NULL, # Diversity clustering penalty parameter
  lambda = NULL, # Ridge regression penalty parameter
  sigma = 0.1, # Width of soft kmeans clusters
  nclust = NULL, # Number of clusters in model
  tau = 0, # Protection against overclustering small datasets with large ones
  block.size = 0.05, # What proportion of cells to update during clustering
  max.iter.harmony = 10L, # Maximum number of rounds to run Harmony
  max.iter.cluster = 20L, # Maximum number of rounds to run clustering at each round of Harmony
  epsilon.cluster = 1e-05, # Convergence tolerance for clustering round of Harmony
  epsilon.harmony = 1e-04, # Convergence tolerance for Harmony
  verbose = TRUE
)

# embedding and clustering ------------------------------------------------
integrated <- integrated %>% 
  FindNeighbors(
    reduction = "harmony", # SdP: cannot be changed
    dims = 1:n_dims, 
    assay = NULL,
    features = NULL, # Features to use as input for building the (S)NN; used only when dims is NULL
    k.param = 20, # Defines k for the k-nearest neighbor algorithm (SdP: change to make it compliant with RunUMAP)
    return.neighbor = FALSE,
    # compute.SNN = !return.neighbor,
    prune.SNN = 1/15, # Sets the cutoff for acceptable Jaccard index when computing the neighborhood overlap for the SNN construction. Any edges with values less than or equal to this will be set to 0 and removed from the SNN graph
    nn.method = "annoy", # Method for nearest neighbor finding. Options include: rann, annoy
    n.trees = 50, # More trees gives higher precision when using annoy approximate nearest neighbor search
    annoy.metric = "euclidean", # Distance metric for annoy. Options include: euclidean, cosine, manhattan, and hamming
    nn.eps = 0, # Error bound when performing nearest neighbor seach using RANN; default of 0.0 implies exact nearest neighbor search
    verbose = TRUE,
    do.plot = FALSE,
    graph.name = NULL, # Default is assay.name_(s)nn
    l2.norm = FALSE, # Take L2Norm of the data
    cache.index = FALSE
  ) %>% 
  RunUMAP(
    dims = 1:n_dims, 
    reduction = "harmony", # SdP: cannot be changed
    features = NULL, # dims must be NULL to run on features
    graph = NULL, # Name of graph on which to run UMAP
    assay = DefaultAssay(object = object),
    nn.name = NULL, # Name of knn output on which to run UMAP
    slot = "data",
    umap.method = "uwot",
    reduction.model = NULL, # DimReduc object that contains the umap model
    return.model = FALSE,
    n.neighbors = 30L, # SdP: same as in FindNeighbors?
    n.components = 2L,
    metric = "cosine",
    n.epochs = NULL,
    learning.rate = 1,
    min.dist = 0.3,
    spread = 1,
    set.op.mix.ratio = 1,
    local.connectivity = 1L,
    repulsion.strength = 1,
    negative.sample.rate = 5L,
    a = NULL,
    b = NULL,
    uwot.sgd = FALSE,
    seed.use = 42L,
    metric.kwds = NULL,
    angular.rp.forest = FALSE,
    densmap = FALSE,
    dens.lambda = 2,
    dens.frac = 0.3,
    dens.var.shift = 0.1,
    verbose = TRUE,
    reduction.name = "umap",
    reduction.key = NULL
  ) %>% 
  FindClusters(
    graph.name = NULL,
    cluster.name = NULL,
    modularity.fxn = 1,
    initial.membership = NULL,
    node.sizes = NULL,
    resolution = resolutions, # SdP: multiple resolutions?
    # argument is now depracated
    # method = "matrix",
    algorithm = 1,
    n.start = 10,
    n.iter = 10,
    random.seed = 0,
    group.singletons = TRUE,
    temp.file.location = NULL,
    edge.file.name = NULL,
    verbose = TRUE
  )

# join layers -------------------------------------------------------------
integrated <- JoinLayers(integrated)

# save output -------------------------------------------------------------
saveRDS(integrated,out_object)
