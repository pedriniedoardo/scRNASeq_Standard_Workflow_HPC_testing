# AIM ---------------------------------------------------------------------
# the aim is to compare runnign the pipeline with or without filtering genes at the level of single samples before integration.
# does it affect the selection of the HVG?

# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(UpSetR)
library(ComplexUpset)

# read in the data --------------------------------------------------------
test_filter <- readRDS("out/object/pbmc_test/connect_5k_pbmc_harmony_filter.rds")
p1 <- DimPlot(test_filter)

test_all <- readRDS("out/object/pbmc_test/connect_5k_pbmc_harmony_all.rds")
p2 <- DimPlot(test_all)

p1+p2

# exploration -------------------------------------------------------------
# compare the HVG between the two implementations
HVG_filter <- VariableFeatures(test_filter)
HVG_all <- VariableFeatures(test_all)

list_dataset <- list(HVG_filter=HVG_filter,HVG_all=HVG_all)

ComplexUpset::upset(fromList(list_dataset),colnames(fromList(list_dataset)),wrap=T) + ggtitle("HVG")

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
bind_rows(
  df_exp_all_all$RNA %>%
    data.frame(avg_count = .) %>%
    rownames_to_column("gene") %>%
    # pivot_longer(names_to = "gene",values_to = "count",-barcode) %>%
    mutate(test= "all"),
  
  df_exp_all_filter$RNA %>%
    data.frame(avg_count = .) %>%
    rownames_to_column("gene") %>%
    # pivot_longer(names_to = "gene",values_to = "count",-barcode) %>%
    mutate(test= "filter")
) %>%
  ggplot(aes(x=test,y=avg_count)) + geom_boxplot(outlier.shape = NA) + scale_y_continuous(trans = "log1p") +
  coord_cartesian(ylim = c(0,0.5)) +
  theme_bw() +
  geom_point(position = position_jitter(width = 0.2),alpha=0.1)

# add also the expression from the common HVG
bind_rows(
  df_exp_all_all$RNA %>%
    data.frame(avg_count = .) %>%
    rownames_to_column("gene") %>%
    # pivot_longer(names_to = "gene",values_to = "count",-barcode) %>%
    mutate(test= "all_unique"),
  
  df_exp_all_filter$RNA %>%
    data.frame(avg_count = .) %>%
    rownames_to_column("gene") %>%
    # pivot_longer(names_to = "gene",values_to = "count",-barcode) %>%
    mutate(test= "filter_unique")
) %>%
  ggplot(aes(x=test,y=avg_count)) + geom_boxplot(outlier.shape = NA) + scale_y_continuous(trans = "log1p") +
  coord_cartesian(ylim = c(0,0.5)) +
  theme_bw() +
  geom_point(position = position_jitter(width = 0.2),alpha=0.1)

# test --------------------------------------------------------------------
# stefano suggested to calculate the common marker per cluster and compare them with the set of HVG
# if there are more HVG in the markers genes, it is better


