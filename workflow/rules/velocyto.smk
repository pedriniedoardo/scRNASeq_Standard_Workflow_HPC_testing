rule runVelocyto:
    '''
    run velocyto on the cellanger output
    '''
    input:
        cellranger_folder = rules.runCellrangerMultiRun.output.folder,
        transcriptome_gtf = config["ref"]["reference"] + "/genes/genes.gtf",
        mask_gtf = config["ref"]["repeat_mask_gtf"]
    output:
        loom = config["out_location"] + "velocyto/merged/{sample_name}/{sample_name}.loom"
    conda:
        config["env_velocyto"]
    log:
        "logs/velocyto/{sample_name}.log"
    threads:
        config["set-threads"]["runVelocyto"]
    resources:
        # Starts at base memory, doubles or adds for each retry
        mem_mb = lambda wildcards, attempt: config["set-resources"]["runVelocyto"]["base_mb"] * attempt
    params:
        # FIX: Check if mem_mb is an integer (real run) or a TBDString (dry run)
        mem_thread = lambda wildcards, resources, threads: int((resources.mem_mb - 2000) / threads) if isinstance(resources.mem_mb, int) else 2000
    shadow: 
        # use shallow to make visible the whole snakemake root
        # use minimal to make visible only the specified inputs
        "minimal"
    shell:
        '''
        echo "Starting velocyto..." > {log}
        
        velocyto run10x \
        -m {input.mask_gtf} \
        --samtools-threads {threads} \
        --samtools-memory {params.mem_thread} \
        {input.cellranger_folder} \
        {input.transcriptome_gtf} \
        >> {log} 2>&1

        # move the output to the destination
        mkdir -p $(dirname {output.loom})
        mv {input.cellranger_folder}/velocyto/*.loom {output.loom}
        '''