rule runCellbenderTestFast:
    input:
        raw_mtx = "/idle/ric.cosr/ric.cosr/pedrini.edoardo/sc-rna-seq-cosr-standard-workflow-test2/workflow/test/cellbender/example/tiny_raw_feature_bc_matrix.h5ad"
    output:
        h5 = "/idle/ric.cosr/ric.cosr/pedrini.edoardo/sc-rna-seq-cosr-standard-workflow-test2/workflow/test/cellbender/example/out_GPU/tiny_output.h5"
    resources:
        slurm_partition="cuda",
        slurm_extra="--gres=gpu:1"
    conda:
        config["env_cellbender"]
    params:
    shell:
        '''
        cellbender remove-background \
        --input {input.raw_mtx} \
        --output {output.h5} \
        --expected-cells 500 \
        --total-droplets-included 2000 \
        --cuda
        '''

rule runCellbenderTestFast2:
    input:
        raw_mtx = "/idle/ric.cosr/ric.cosr/pedrini.edoardo/sc-rna-seq-cosr-standard-workflow-test2/workflow/test/cellbender/example/tiny_raw_feature_bc_matrix.h5ad"
    output:
        h5 = "/idle/ric.cosr/ric.cosr/pedrini.edoardo/sc-rna-seq-cosr-standard-workflow-test2/workflow/test/cellbender/example/out_GPU_minimal/tiny_output.h5"
    resources:
        slurm_partition="cuda",
        slurm_extra="--gres=gpu:1"
    conda:
        config["env_cellbender"]
    shadow: "minimal"
    params:
    shell:
        '''
        cellbender remove-background \
        --input {input.raw_mtx} \
        --output {output.h5} \
        --expected-cells 500 \
        --total-droplets-included 2000 \
        --cuda
        '''

rule runCellbenderTest:
    input:
        raw_mtx = "/idle/ric.cosr/ric.cosr/pedrini.edoardo/sc-rna-seq-cosr-standard-workflow-test2/workflow/test/cellbender/test_connect_5k_pbmc/raw_feature_bc_matrix.h5"
    output:
        h5 = "/idle/ric.cosr/ric.cosr/pedrini.edoardo/sc-rna-seq-cosr-standard-workflow-test2/workflow/test/cellbender/test_connect_5k_pbmc/cellbender_output.h5"
    resources:
        slurm_partition="cuda",
        slurm_extra="--gres=gpu:1"
    conda:
        config["env_cellbender"]
    shadow: "minimal"
    params:
    shell:
        '''
        cellbender remove-background \
        --input {input.raw_mtx} \
        --output {output.h5} \
        --cuda
        '''
