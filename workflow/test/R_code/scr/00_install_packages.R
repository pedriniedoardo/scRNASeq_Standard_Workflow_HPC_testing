# match the same version set in the yaml file
renv::install(c("scater@1.34.1",
                "tidyverse@2.0.0",
                "robustbase@0.99-6",
                "patchwork@1.3.2",
                "scuttle@1.16.0",
                "immunogenomics/presto@1.0.0",
                "Seurat@5.3.1",
                "harmony@1.2.4"))

# cannot find the specific version definded in the yaml (scdblfinder=1.23.4).
# I will install the latest available for this version of bioconductor 3.20
renv::install("scDblFinder")
