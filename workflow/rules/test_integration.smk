rule runTestSeuratAll:
    '''
    this rule is running the integration without filtering
    '''
    input:
        spe_input = config["location_benchmark_obj"] + "{sample_benchmark}_listTest.rds"
    output:
        spe_output = config["location_benchmark_obj"] + "01_{sample_benchmark}_all.rds"
    conda:
        config["env_seurat"]
    log:
        'logs/test_benchmark/{sample_benchmark}/runTestSeuratAll.log'
    benchmark:
        'benchmarks/test_benchmark/{sample_benchmark}/runTestSeuratAll.txt'
    script:
        "../test/R_code/scr/snakemake/01_integrate_geneFilter_all_snakemake.R"

rule runTestSeuratFilter:
    '''
    this rule is running the integration with filtering
    '''
    input:
        spe_input = config["location_benchmark_obj"] + "{sample_benchmark}_listTest.rds"
    output:
        spe_output = config["location_benchmark_obj"] + "01_{sample_benchmark}_filter.rds"
    conda:
        config["env_seurat"]
    log:
        'logs/test_benchmark/{sample_benchmark}/runTestSeuratFilter.log'
    benchmark:
        'benchmarks/test_benchmark/{sample_benchmark}/runTestSeuratFilter.txt'
    script:
        "../test/R_code/scr/snakemake/01_integrate_geneFilter_filter_snakemake.R"

rule runTestSeuratCompare:
    '''
    Compare pipeline results with and without gene filtering:
    Evaluates HVG overlap, Jaccard similarity of clusters, and marker gene retention.
    '''
    input:
        spe_all    = rules.runTestSeuratAll.output.spe_output,
        spe_filter = rules.runTestSeuratFilter.output.spe_output
    output:
        plot_umap = config["location_benchmark_plot"] + "04_{sample_benchmark}_UMAP_{res}.pdf",
        plot_upset = config["location_benchmark_plot"] + "04_{sample_benchmark}_UpsetPlot_{res}.pdf",
        plot_exp = config["location_benchmark_plot"] + "04_{sample_benchmark}_BoxplotGeneExp_{res}.pdf",
        plot_jaccard = config["location_benchmark_plot"] + "04_{sample_benchmark}_Jaccard_{res}.pdf",
        plot_prop_hvg = config["location_benchmark_plot"] + "04_{sample_benchmark}_BoxplotPropHVG_{res}.pdf",
        tab_markers_all = config["location_benchmark_table"] + "04_{sample_benchmark}_markersAll_{res}.tsv",
        tab_markers_filt = config["location_benchmark_table"] + "04_{sample_benchmark}_markersFilter_{res}.tsv"
    conda:
        config["env_seurat"]
    params:
        res = "{res}",
        nm = "{sample_benchmark}"
    log:
        'logs/test_benchmark/{sample_benchmark}/runTestSeuratCompare_{res}.log'
    benchmark:
        'benchmarks/test_benchmark/{sample_benchmark}/runTestSeuratCompare_{res}.txt'
    script:
        # Assuming you save the R snippet to this location
        "../test/R_code/scr/snakemake/02_integrate_geneFilter_compare_snakemake.R"