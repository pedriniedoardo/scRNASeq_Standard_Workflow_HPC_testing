rule reportImages:
  '''
  This rule is intended for producing images on the embedded Seurat's object
  '''
  input:
    merged_rds     = rules.integrateSamples.output.merged_rds,
    integrated_rds = rules.integrateSamples.output.integrated_rds
  output:
    vln_qc = expand(
        config["out_location"] + "Seurat/plot/vln_{feature}.{extension}",
        feature   = config["report"]["qc_features"],
        extension = config["report"]["figure_extension"]
    ),
    dim_qc = expand(
        config["out_location"] + "Seurat/plot/dim_{feature}.{extension}",
        feature   = config["report"]["qc_features"],
        extension = config["report"]["figure_extension"]
    ),
    dim_ident = expand(
      config["out_location"] + "Seurat/plot/dim_{object_type}_ident.{extension}",
      object_type = ["merged", "integrated"],
      extension   = config["report"]["figure_extension"]
      ),
    dim_clusters = expand(
      config["out_location"] + "Seurat/plot/dim_{object_type}_cluster.{extension}",
      object_type = ["merged", "integrated"],
      extension   = config["report"]["figure_extension"]
      )
  params:
    out_location     = config["out_location"],
    qc_features      = config["report"]["qc_features"],
    figure_extension = config["report"]["figure_extension"],
    resolutions      = config["Seurat"]["embedding"]["resolutions"]
  conda: config["env_seurat"]
  log:
      'logs/report/reportImages.log'
  benchmark:
      'benchmarks/report/reportImages.txt'
  script:
      "../scripts/05_report_images.R"

    
