# AIM ---------------------------------------------------------------------
# for benchmakrs purposes, try to generate fake samples from individual samples.
# save the object as a list

# libraries ---------------------------------------------------------------
library(tidyverse)
library(Seurat)

# processing --------------------------------------------------------------
# define the folders
folder_input <- "data/sample_obj/test_human_default/"
folder_output <- "out/object/test_human_default/"

# define the samples.
# notice all of them are single samples
sample_id <- dir(folder_input) %>%
  str_subset(pattern = "postQC.rds") %>%
  str_subset(pattern = "connect_5k_pbmc_NGSC3",negate = T)

# generate the replicates
lapply(sample_id,function(file){
  # track the progress
  print(file)
  
  # generate fake replicates for testing the integration procedure
  # read in the dataset
  test <- readRDS(paste0(folder_input,file))
  
  # create the fake sample id
  test$sample <- test@meta.data %>%
    mutate(sample = sample(c("sample_0","sample_1","sample_2"),replace = T,size = nrow(.))) %>%
    pull(sample)
  
  # pull the labels of the ids
  sample_id <- test$sample %>% unique()
  
  # subse the daset
  list_sobj <- lapply(sample_id,function(id){
    test <- subset(test,sample == id)
    return(test)
  })
  
  # define the name of the list of object
  list_id <- str_replace_all(file,pattern = "default_obj_postQC",replacement = "listTest")
  
  # save the list
  saveRDS(list_sobj,paste0(folder_output,list_id))
})

# -------------------------------------------------------------------------
# handle the real replicated from the PBMCs dataset generating list from the two
# define the folders
folder_input <- "data/sample_obj/test_human_default/"
folder_output <- "out/object/test_human_default/"

# define the samples.
# notice all of them are single samples
sample_id <- dir(folder_input) %>%
  str_subset(pattern = "postQC.rds") %>%
  str_subset(pattern = "connect_5k_pbmc_NGSC3",negate = F)

# generate the replicates
list_pbmc <- lapply(sample_id,function(file){
  # track the progress
  print(file)
  
  # generate fake replicates for testing the integration procedure
  # read in the dataset
  test <- readRDS(paste0(folder_input,file))
  
  return(test)
})

# save the list
saveRDS(list_pbmc,paste0(folder_output,"connect_5k_pbmc_NGSC3_ch1_gex_listTest.rds"))

# -------------------------------------------------------------------------
