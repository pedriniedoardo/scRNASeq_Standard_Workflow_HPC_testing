rule cellranger_mkref:
    """
    This rule builds a Cell Ranger reference by combining transcriptome FASTA/GTF
    with a transgene FASTA (and auto-generating a minimal GTF).
    """
    input:
        transcriptome_fasta = config["ref"]["reference_fasta"],
        transcriptome_gtf = config["ref"]["reference_gtf"],
        transgene_fasta = config["trans"]["transgene_fasta"]
    output:
        ref_dir = directory(
            f"reference/{config['ref']['build']}_{config['trans']['transgene_name']}_cellranger_ref"
        )
    conda:
        config["env_cellranger"]
    log:
        f"logs/reference/{config['ref']['build']}_{config['trans']['transgene_name']}/cellranger_mkref.log"
    benchmark:
        f"benchmarks/reference/{config['ref']['build']}_{config['trans']['transgene_name']}/cellranger_mkref.txt"
    resources:
        cpus = config["set-threads"]["cellranger_mkref"]
    threads:
        config["set-threads"]["cellranger_mkref"]
    params:
        combined_fa = f"reference/{config['ref']['build']}_{config['trans']['transgene_name']}.fa",
        combined_gtf = f"reference/{config['ref']['build']}_{config['trans']['transgene_name']}.gtf",
        transgene_gtf = f"reference/{config['trans']['transgene_name']}.gtf",
        genome_name = f"{config['ref']['build']}_{config['trans']['transgene_name']}",
        RAM = config["set-resources"]["cellranger_mkref"]["mem_gb"],
        cpus = config["set-threads"]["cellranger_mkref"]
    shell:
        r'''
        echo "===== Building Cell Ranger Reference =====" >> {log}
        mkdir -p reference >> {log}

        echo "Merging transcriptome and transgene FASTA..." >> {log}
        cat {input.transcriptome_fasta} {input.transgene_fasta} > {params.combined_fa}

        echo "Auto-generating transgene GTF..." >> {log}
        transgene_id=$(grep ">" {input.transgene_fasta} | sed 's/>//; s/ .*//')
        seq_length=$(grep -v ">" {input.transgene_fasta} | tr -d '\n' | wc -c)
        echo -e "${{transgene_id}}\ttransgene\texon\t1\t${{seq_length}}\t.\t+\t.\tgene_id \"${{transgene_id}}\"; transcript_id \"${{transgene_id}}\"; gene_name \"${{transgene_id}}\";" > {params.transgene_gtf}

        echo "Merging transcriptome and transgene GTFs..." >> {log}
        awk '$1 !~ /^#/' {input.transcriptome_gtf} > {params.combined_gtf}
        awk '$1 !~ /^#/' {params.transgene_gtf} >> {params.combined_gtf}

        echo "Running cellranger mkref..." >> {log}
        cellranger mkref \
            --genome={params.genome_name} \
            --fasta={params.combined_fa} \
            --genes={params.combined_gtf} \
            --nthreads={params.cpus} \
            --memgb={params.RAM} \
            --output-dir={output.ref_dir}

        echo "===== Reference built successfully at {output.ref_dir} =====" >> {log}
        '''

