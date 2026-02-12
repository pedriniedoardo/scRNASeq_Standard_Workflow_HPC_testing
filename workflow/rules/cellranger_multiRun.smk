rule runCellrangerMultiRun:
    '''
    This is the rule to run cellranger from multi run samples
    The way this is implemented now I cannot validate the presence of the input files
    '''
    input:
        transcriptome = config["ref"]["reference"],
        # this way I am testing for the existance of the actual input files
        fastqs = lambda wildcards: SAMPLES_merge[wildcards.sample_name]["sample_path_list"],
    output:
        output = config["out_location"] + "cellranger/merged/{sample_name}_finished.log",
        folder = directory(config["out_location"] + "cellranger/merged/{sample_name}")
    conda: config["env_cellranger"]
    log:
        'logs/merged/{sample_name}/runCellrangerMultiRun.log'
    benchmark:
        'benchmarks/merged/{sample_name}/runCellrangerMultiRun.txt'
    resources:
        # mem_mb = 16000,
        cpus = config["set-threads"]["runCellranger"]
    threads:
        config["set-threads"]["runCellranger"]
    params:
        fastqs = lambda wildcards: SAMPLES_merge[wildcards.sample_name]["sample_path_string"],
        INTRON_flag = config['INTRON_flag'],
        BAM_flag = config['BAM_flag'],
        cpus = config["set-threads"]["runCellranger"],
        RAM = config["set-resources"]["runCellranger"]["mem_gb"],
        INTRON_param = config['INTRON_param'],
        outs = directory(config["out_location"] + "cellranger/merged/{sample_name}")
    shell:
        '''
        echo "generate the output folder" >> {log}
        mkdir -p {params.outs} >> {log}
        
        echo "start cellranger run for <{wildcards.sample_name}>" >> {log}

        cellranger count --id={wildcards.sample_name} \
        --transcriptome={input.transcriptome} \
        --fastqs={params.fastqs} \
        --sample={wildcards.sample_name} \
        --localcores={params.cpus} \
        --create-bam={params.BAM_flag} \
        {params.INTRON_param}={params.INTRON_flag} \
        --localmem={params.RAM} \
        --output-dir={params.outs}
        
        # generate the output file
        touch {output.output}

        echo "completed run cellranger for <{wildcards.sample_name}>" >> {log}
        '''

rule moveCellrangerSummaryMultiRun:
    '''
    This is the rule to move the outputs of interest.
    '''
    input:
        cellranger_summary = rules.runCellrangerMultiRun.output.output
    output:
        summary = config["out_location"] + "web_summaries/merged/{sample_name}_web_summary.html"
    # conda: config["env_bioinfo"]
    log:
        'logs/merged/{sample_name}/moveCellrangerSummaryMultiRun.log'
    benchmark:
        'benchmarks/merged/{sample_name}/moveCellrangerSummaryMultiRun.txt'
    resources:
        mem_mb = 250,
        cpus = 1
    threads: 1
    params:
        wd_summary = config["out_location"] + 'cellranger/merged/{sample_name}/outs/web_summary.html',
    shell:
        '''
        # copy
        echo "copy summary <{wildcards.sample_name}>" >> {log}

        # copy the web summary from the output in an individula folder
        cp {params.wd_summary} {output.summary} >> {log}

        echo "summary copied <{wildcards.sample_name}>" >> {log}
        '''

rule runMultiqc1MultiRun:
    '''
    This rule allow to run the multiqc on the final folder.
    Currently it is set to run after cellranger. this is the only step to grep some outputs.
    '''
    input:
        # this is needed to trigger it after the generation of the outputs on all the outputs
        input_file=expand(rules.moveCellrangerSummaryMultiRun.output.summary,
        sample_name=SAMPLES_merge.keys())
    output:
        html = config["out_location"] + "multiQC/merged/multiqc_report.html"
    conda: config["env_bioinfo"]
    log:
        'logs/merged/runMultiqc1MultiRun.log'
    benchmark:
        'benchmarks/merged/runMultiqc1MultiRun.txt'
    resources:
        mem_gb = 4,
        cpus = 2
    threads: 2
    params:
        folder_in = config["out_location"] + "cellranger/merged/",
        folder_out = config["out_location"] + "multiQC/merged/"
    shell:
        '''
        echo "start multiqc" >> {log}
        multiqc --force {params.folder_in} -o {params.folder_out}
        echo "end multiqc" >> {log}
        '''
