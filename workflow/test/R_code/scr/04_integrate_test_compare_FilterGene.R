# AIM ---------------------------------------------------------------------
# the aim is to compare runnign the pipeline with or without filtering genes at the level of single samples before integration.
# does it affect the selection of the HVG?

# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(UpSetR)
library(ComplexUpset)
library(ComplexHeatmap)

# custom functions --------------------------------------------------------
# define the jaccard score function
jaccard <- function(a, b) {
  intersection <- length(intersect(a, b))
  union <- length(a) + length(b) - intersection
  return (intersection/union)
}

# test
a <- c('potato', 'tomotto', 'chips', 'baloon')
b <- c('car', 'chips', 'bird', 'salt')

jaccard(a, b)

# parameters --------------------------------------------------------------
# define the resolution parameter for the comparison
res <- "RNA_snn_res.0.3"

spe_filter <- "out/object/test_human_default/01_connect_5k_pbmc_NGSC3_ch1_gex_filter.rds"
spe_all <- "out/object/test_human_default/01_connect_5k_pbmc_NGSC3_ch1_gex_all.rds"

nm <- "connect_5k_pbmc_NGSC3_ch1_gex"

p03_out <- paste0("out/plot/04_connect_5k_pbmc_NGSC3_ch1_gex_UMAP_",res,".pdf")
p04_out <- "out/plot/04_connect_5k_pbmc_NGSC3_ch1_gex_UpsetPlot.pdf"
p05_out <- "out/plot/04_connect_5k_pbmc_NGSC3_ch1_gex_BarGeneExp.pdf"
ht02_out <- paste0("out/plot/04_connect_5k_pbmc_NGSC3_ch1_gex_Jaccard_",res,".pdf")
p06_out <- paste0("out/plot/04_connect_5k_pbmc_NGSC3_ch1_gex_BarPropHVG_",res,".pdf")

tab01_all <- paste0("out/table/04_connect_5k_pbmc_NGSC3_ch1_gex_markersAll_",res,".tsv")
tab01_filter <- paste0("out/table/04_connect_5k_pbmc_NGSC3_ch1_gex_markersFilter_",res,".tsv")

# read in the data --------------------------------------------------------
test_filter <- readRDS(spe_filter)
p01 <- DimPlot(test_filter, group.by = res) + ggtitle("test_filter")

test_all <- readRDS(spe_all)
p02 <- DimPlot(test_all, group.by = res) + ggtitle("test_all")

# make the plots togetegr
p03 <- (p01+p02) + plot_annotation(nm)
ggsave(plot = p03,filename = p03_out,width = 10,height = 4)

# exploration -------------------------------------------------------------
# compare the HVG between the two implementations
HVG_filter <- VariableFeatures(test_filter)
HVG_all <- VariableFeatures(test_all)

list_dataset <- list(HVG_filter=HVG_filter,HVG_all=HVG_all)

p04 <- ComplexUpset::upset(fromList(list_dataset),colnames(fromList(list_dataset)),wrap=T) + ggtitle(paste("HVG",nm))
ggsave(plot = p04,filename = p04_out,width = 6,height = 4)

# pull the different ones
# pull all the genes from each term
df1 <- lapply(list_dataset,function(x){
  data.frame(gene = x)
}) %>% 
  bind_rows(.id = "path")

#
head(df1)

# pull all the unique genes
df2 <- data.frame(gene=unique(unlist(list_dataset)))

#
head(df2)

# generate the intersections
df_int <- lapply(df2$gene,function(x){
  # pull the name of the intersections
  intersection <- df1 %>% 
    dplyr::filter(gene==x) %>% 
    arrange(path) %>% 
    pull("path") %>% 
    paste0(collapse = "|")
  
  # build the dataframe
  data.frame(gene = x,int = intersection)
}) %>% 
  bind_rows()

#
head(df_int,n=20)

# show the intersections
df_int %>% 
  group_by(int) %>% 
  summarise(n=n()) %>% 
  arrange(desc(n))

# check the expression of the HVG_all unique
# add a covariate for the whols sample
test_all$test <- "test_all"
df_exp_all_all <- AverageExpression(test_all,features = df_int %>% filter(int == "HVG_all") %>% pull(gene),layer = "counts",group.by = "test")
df_exp_all_filter <- AverageExpression(test_all,features = df_int %>% filter(int == "HVG_filter") %>% pull(gene),layer = "counts",group.by = "test")
df_exp_all_common <- AverageExpression(test_all,features = df_int %>% filter(int == "HVG_all|HVG_filter") %>% pull(gene),layer = "counts",group.by = "test")

test_all@assays$RNA$counts[1:10,1:10]

# compare the average expression from all the intersections
# bind_rows(
#   df_exp_all_all$RNA %>%
#     data.frame(avg_count = .) %>%
#     rownames_to_column("gene") %>%
#     # pivot_longer(names_to = "gene",values_to = "count",-barcode) %>%
#     mutate(test= "all_unique"),
#   
#   df_exp_all_filter$RNA %>%
#     data.frame(avg_count = .) %>%
#     rownames_to_column("gene") %>%
#     # pivot_longer(names_to = "gene",values_to = "count",-barcode) %>%
#     mutate(test= "filter_unique")
# ) %>%
#   ggplot(aes(x=test,y=avg_count)) + geom_boxplot(outlier.shape = NA) + scale_y_continuous(trans = "log1p") +
#   coord_cartesian(ylim = c(0,0.5)) +
#   theme_bw() +
#   geom_point(position = position_jitter(width = 0.2),alpha=0.1)

# add also the expression from the common HVG
df_exp_all_HVG <- bind_rows(
  df_exp_all_all$RNA %>%
    data.frame(avg_count = .) %>%
    rownames_to_column("gene") %>%
    # pivot_longer(names_to = "gene",values_to = "count",-barcode) %>%
    mutate(test= "all_unique"),
  
  df_exp_all_filter$RNA %>%
    data.frame(avg_count = .) %>%
    rownames_to_column("gene") %>%
    # pivot_longer(names_to = "gene",values_to = "count",-barcode) %>%
    mutate(test= "filter_unique"),
  
  df_exp_all_common$RNA %>%
    data.frame(avg_count = .) %>%
    rownames_to_column("gene") %>%
    # pivot_longer(names_to = "gene",values_to = "count",-barcode) %>%
    mutate(test= "common")
)

# define the limits for the plotting
thr_min <- 0
thr_max <- round(summary(df_exp_all_HVG$avg_count)[5]*4,digits = 1)

# make the plot
p05 <- df_exp_all_HVG %>%
  mutate(test = factor(test,levels = c("common","filter_unique","all_unique"))) %>%
  ggplot(aes(x=test,y=avg_count)) + geom_boxplot(outlier.shape = NA) + scale_y_continuous(trans = "log1p") +
  coord_cartesian(ylim = c(thr_min,thr_max)) +
  theme_bw() +
  geom_point(position = position_jitter(width = 0.2),alpha=0.1) +
  ggtitle(nm)
ggsave(plot = p05,filename = p05_out,width = 5,height = 5)

# test cluster similarity -------------------------------------------------
# check how similar are the clusters using Jaccard score

# pull the meta
meta_all <- test_all@meta.data %>%
  rownames_to_column("barcode") %>%
  mutate(cluster_test = paste0("all-",RNA_snn_res.0.3))

meta_filter <- test_filter@meta.data %>%
  rownames_to_column("barcode") %>%
  mutate(cluster_test = paste0("filt-",RNA_snn_res.0.3))

# build the dataset for the correlatino plot
df_crossing <- crossing(seurat_clusters_all = unique(meta_all$cluster_test),
                        seurat_clusters_filter = unique(meta_filter$cluster_test))

# build the scatter plot
df_jaccard_score <- pmap(list(all = df_crossing$seurat_clusters_all,
                              filt = df_crossing$seurat_clusters_filter), function(all,filt){
                                
                                # calculate the jaccard score
                                a <- meta_all %>%
                                  filter(cluster_test == all) %>% pull(barcode)
                                b <- meta_filter %>%
                                  filter(cluster_test == filt) %>% pull(barcode)
                                jaccard_score <- jaccard(a,b)
                                
                                # build a data.frame
                                df <- data.frame(all = all,
                                                 filt = filt,
                                                 jaccard_score = jaccard_score)
                                return(df)
                              }) %>%
  bind_rows()

# shape it as a matrix
mat_jaccard_score <- df_jaccard_score %>%
  pivot_wider(names_from = all,values_from = jaccard_score) %>%
  column_to_rownames("filt")

ht02 <- Heatmap(mat_jaccard_score,
                 name = "Jaccard\nscore",
                 column_title = nm,
                 # col = colorRamp2(c(-1, 0, 1), colors = c("blue", "white", "red")),
                 col = viridis::viridis(option = "turbo",n = 20),
                 row_names_side = "right",
                 row_names_gp = gpar(fontsize = 8),
                 column_names_side = "bottom",
                 column_names_gp = gpar(fontsize = 8),
                 row_dend_reorder = FALSE,
                 column_dend_reorder = FALSE,
                 row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                 column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                 show_column_names = T,
                 show_row_names = T)


pdf(file = ht02_out,width = 7,height = 6)
draw(ht02)
dev.off()

# test markers ------------------------------------------------------------
# stefano suggested to calculate the common marker per cluster and compare them with the set of HVG
# if there are more HVG in the markers genes, it is better

Idents(test_all) <- res
markers.all <- FindAllMarkers(test_all)

Idents(test_filter) <- res
markers.filter <- FindAllMarkers(test_filter)

# save the marlers
markers.all %>%
  write_tsv(tab01_all)

markers.filter %>%
  write_tsv(tab01_filter)

# pull the top 100 markers and compare how many are in the variable genes
markers.all_fix <- markers.all %>%
  group_by(cluster) %>%
  dplyr::slice(1:100) %>%
  split(f = .$cluster) %>%
  lapply(function(x){
    x %>%
      mutate(isHVG = gene %in% HVG_all)
  }) %>%
  bind_rows() %>%
  ungroup()

# do the same for the filtered dataset
markers.filter_fix <- markers.filter %>%
  group_by(cluster) %>%
  dplyr::slice(1:100) %>%
  split(f = .$cluster) %>%
  lapply(function(x){
    x %>%
      mutate(isHVG = gene %in% HVG_filter)
  }) %>%
  bind_rows() %>%
  ungroup()

# make a general assessment
markers.all_fix %>%
  summarise(avg_isHVG = mean(isHVG))
markers.filter_fix %>%
  summarise(avg_isHVG = mean(isHVG))

# plot the proportion per cluster
markers.all_summary <- markers.all_fix %>%
  group_by(cluster) %>%
  summarise(prop_isHVG = mean(isHVG)) %>%
  mutate(test = "all")

markers.filter_summary <- markers.filter_fix %>%
  group_by(cluster) %>%
  summarise(prop_isHVG = mean(isHVG)) %>%
  mutate(test = "filter")

# plot the different proportions
p06 <- bind_rows(markers.filter_summary,
          markers.all_summary) %>%
  ggplot(aes(x=test,y=prop_isHVG)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2),shape = 1) + theme_bw() +
  ggtitle(nm)
ggsave(plot = p06,filename = p06_out,width = 5,height = 5)
