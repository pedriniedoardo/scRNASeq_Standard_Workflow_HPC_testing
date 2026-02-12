rule runCellranger:
    '''
    This is the rule to run cellranger.
    '''
    input:
        fastqs = lambda wildcards: SAMPLES_single[(wildcards.sample_name, wildcards.sample_run)]["sample_path_list"],
        transcriptome = config["ref"]["reference"],
    output:
        output = config["out_location"] + "cellranger/single/{sample_run}/{sample_name}_finished.log",
        folder = directory(config["out_location"] + "cellranger/single/{sample_run}/{sample_name}")
    conda: config["env_cellranger"]
    log:
        'logs/single/{sample_run}/{sample_name}/runCellranger.log'
    benchmark:
        'benchmarks/single/{sample_run}/{sample_name}/runCellranger.txt'
    resources:
        # mem_mb = 16000,
        cpus = config["set-threads"]["runCellranger"]
    threads:
        config["set-threads"]["runCellranger"]
    params:
        fastqs = lambda wildcards: SAMPLES_single[(wildcards.sample_name, wildcards.sample_run)]["sample_path_string"],
        sample_name_unique = lambda wildcards: SAMPLES_single[(wildcards.sample_name, wildcards.sample_run)]["sample_name_unique"],
        INTRON_flag = config['INTRON_flag'],
        BAM_flag = config['BAM_flag'],
        cpus = config["set-threads"]["runCellranger"],
        RAM = config["set-resources"]["runCellranger"]["mem_gb"],
        INTRON_param = config['INTRON_param'],
        outs = directory(config["out_location"] + "cellranger/single/{sample_run}/{sample_name}")
    shell:
        '''
        echo "generate the output folder" >> {log}
        mkdir -p {params.outs} >> {log}
        
        echo "start cellranger run for <{wildcards.sample_name}>" >> {log}

        cellranger count --id={params.sample_name_unique} \
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

rule moveCellrangerSummary:
    '''
    This is the rule to move the outputs of interest.
    '''
    input:
        cellranger_summary = rules.runCellranger.output.output
    output:
        summary = config["out_location"] + "web_summaries/single/{sample_run}/{sample_name}_web_summary.html"
    # conda: config["env_bioinfo"]
    log:
        'logs/single/{sample_run}/{sample_name}/moveCellrangerSummary.log'
    benchmark:
        'benchmarks/single/{sample_run}/{sample_name}/moveCellrangerSummary.txt'
    resources:
        mem_mb = 250,
        cpus = 1
    threads: 1
    params:
        wd_summary = config["out_location"] + 'cellranger/single/{sample_run}/{sample_name}/outs/web_summary.html',
    shell:
        '''
        # copy
        echo "copy summary <{wildcards.sample_name}>" >> {log}

        # copy the web summary from the output in an individula folder
        cp {params.wd_summary} {output.summary} >> {log}

        echo "summary copied <{wildcards.sample_name}>" >> {log}
        '''

# rule runMultiqc1:
#     '''
#     This rule allow to run the multiqc on the final folder.
#     Currently it is set to run after cellranger. this is the only step to grep some outputs.
#     This implementation allow to have one report per run.
#     '''
#     input:
#         # this is needed to trigger it after the generation of the outputs on all the outputs
#         input_file=expand(rules.moveCellrangerSummary.output.summary,
#                             zip, # Use zip mode
#                             sample_name=[k[0] for k in SAMPLES_single.keys()], # List of first elements (names)
#                             sample_run=[k[1] for k in SAMPLES_single.keys()])
#     output:
#         html = config["out_location"] + "multiQC/single/{sample_run}/multiqc_report.html"
#     conda: config["env_bioinfo"]
#     log:
#         'logs/single/{sample_run}/runMultiqc1.log'
#     benchmark:
#         'benchmarks/single/{sample_run}/runMultiqc1.txt'
#     resources:
#         mem_gb = 4,
#         cpus = 2
#     threads: 2
#     params:
#         folder_in = config["out_location"] + "cellranger/single/{sample_run}/",
#         folder_out = config["out_location"] + "multiQC/single/{sample_run}/"
#     shell:
#         '''
#         echo "start multiqc" >> {log}
#         multiqc {params.folder_in} -o {params.folder_out}
#         echo "end multiqc" >> {log}
#         '''

rule runMultiqc1:
    '''
    This rule allow to run the multiqc on the final folder.
    Currently it is set to run after cellranger. this is the only step to grep some outputs.
    This implementation allow to have all the runs in the same report.
    '''
    input:
        # this is needed to trigger it after the generation of the outputs on all the outputs
        input_file=expand(rules.moveCellrangerSummary.output.summary,
                            zip, # Use zip mode
                            sample_name=[k[0] for k in SAMPLES_single.keys()], # List of first elements (names)
                            sample_run=[k[1] for k in SAMPLES_single.keys()])
    output:
        html = config["out_location"] + "multiQC/single/multiqc_report.html"
    conda: config["env_bioinfo"]
    log:
        'logs/single/runMultiqc1.log'
    benchmark:
        'benchmarks/single/runMultiqc1.txt'
    resources:
        mem_gb = 4,
        cpus = 2
    threads: 2
    params:
        folder_in = config["out_location"] + "cellranger/single/",
        folder_out = config["out_location"] + "multiQC/single/"
    shell:
        '''
        echo "start multiqc" >> {log}
        multiqc --force {params.folder_in} -o {params.folder_out}
        echo "end multiqc" >> {log}
        '''
