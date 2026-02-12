rule runCellbender:
    '''
    This rule run cellbender on the output of cellranger to remove the soup
    '''
    input:
        result_folder = rules.runCellrangerMultiRun.output.folder
    output:
        full_h5 = config["out_location"] + "cellbender/merged/{sample_name}/cellbender_out.h5",
        filter_h5 = config["out_location"] + "cellbender/merged/{sample_name}/cellbender_out_filtered.h5",
        summary_plots = config["out_location"] + "cellbender/merged/{sample_name}/cellbender_out.pdf",
        barcodes = config["out_location"] + "cellbender/merged/{sample_name}/cellbender_out_cell_barcodes.csv",
        metrics = config["out_location"] + "cellbender/merged/{sample_name}/cellbender_out_metrics.csv"
    resources:
        slurm_partition="cuda",
        slurm_extra="--gres=gpu:1"
    conda:
        config["env_cellbender"]
    shadow: "minimal"
    params:
        raw_h5 = config["out_location"] + "cellranger/merged/{sample_name}/outs/raw_feature_bc_matrix.h5"
    shell:
        '''
        cellbender remove-background \
        --input {params.raw_h5} \
        --output {output.full_h5} \
        --cuda
        '''