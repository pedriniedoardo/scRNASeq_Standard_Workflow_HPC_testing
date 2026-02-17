rule runVelocytoTest:
    '''
    run velocyto on the cellanger output
    '''
    input:
        cellranger_folder = "/idle/ric.cosr/ric.cosr/pedrini.edoardo/sc-rna-seq-cosr-standard-workflow-test2/workflow/test/velocyto/connect_5k_pbmc_NGSC3_ch1_gex_1",
        transcriptome_gtf = config["ref"]["reference"] + "/genes/genes.gtf",
        mask_gtf = config["ref"]["repeat_mask_gtf"]
    output:
        loom = "/idle/ric.cosr/ric.cosr/pedrini.edoardo/sc-rna-seq-cosr-standard-workflow-test2/workflow/test/velocyto/connect_5k_pbmc_NGSC3_ch1_gex_1/velocyto/connect_5k_pbmc_NGSC3_ch1_gex_1.loom"
    conda:
        config["env_velocyto"]
    log:
        "logs/velocyto/connect_5k_pbmc_NGSC3_ch1_gex_1.log"
    params:
        cpus = config["set-threads"]["runVelocyto"],
        # RAM = config["set-resources"]["runVelocyto"]["mem_mb"]
    shell:
        '''
        echo "Starting velocyto..." > {log}
        
        velocyto run10x \
        -m {input.mask_gtf} \
        --samtools-threads {params.cpus} \
        --samtools-memory {params.RAM} \
        {input.cellranger_folder} \
        {input.transcriptome_gtf} \
        >> {log} 2>&1
        '''
