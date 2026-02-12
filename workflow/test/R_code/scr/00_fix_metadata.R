# -------------------------------------------------------------------------
# add dataset_id

library(tidyverse)

df <- read_csv("data/test_dataset_10x.csv")

df_full <- df %>%
  mutate(dataset_id = paste(species,source,tissue,product,sep = "_"))

df_full %>%
  write_csv("data/test_dataset_10x_full.csv")

# -------------------------------------------------------------------------
# pull the paths
df <- read_tsv("data/path.txt",col_names = F) %>%
  dplyr::rename(full_path = X1)

df_full <- df %>%
  mutate(file = basename(full_path)) %>%
  mutate(path = dirname(full_path))

# make it minimal
df_full %>%
  group_by(path) %>%
  summarise() %>%
  write_csv("data/path_summary.txt")
