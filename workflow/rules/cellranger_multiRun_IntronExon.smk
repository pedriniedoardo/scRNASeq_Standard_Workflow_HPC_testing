rule moveCellrangerIntronExonMultiRun:
    '''
    This is the rule to move the outputs of interest.
    '''
    input:
        cellranger_summary = rules.runCellrangerMultiRun.output.output
    output:
        bam_file = config["out_location"] + "PureIntronExon/merged/{sample_name}/possorted_genome_bam.bam"
    # conda: config["env_bioinfo"]
    log:
        'logs/merged/{sample_name}/moveCellrangerIntronExonMultiRun.log'
    benchmark:
        'benchmarks/merged/{sample_name}/moveCellrangerIntronExonMultiRun.txt'
    resources:
        mem_mb = 250,
        cpus = 1
    threads: 1
    params:
        bam = config["out_location"] + "cellranger/merged/{sample_name}/outs/possorted_genome_bam.bam"
    shell:
        '''
        # copy
        echo "move the BAM file <{wildcards.sample_name}>" >> {log}

        # move the bam file to the folder of interest for the next step
        mv {params.bam} {output.bam_file} >> {log}

        echo "moved BAM file <{wildcards.sample_name}>" >> {log}
        '''

rule copyGtfMultiRun:
    '''
    thi rule allow to copy the gtf file in the pure intron exon folder to allow the correct anntation of the features
    '''
    input:
        transcriptome_gtf = config["ref"]["reference"] + "/genes/genes.gtf"
    output:
        gtf = config["out_location"] + "PureIntronExon/merged/genes.gtf"
    # conda: config["env_bioinfo"]
    log:
        'logs/merged/copyGtfMultiRun.log'
    benchmark:
        'benchmarks/merged/copyGtfMultiRun.txt'
    resources:
        mem_mb = 250,
        cpus = 1
    threads: 1
    params:
    shell:
        '''
        echo "start copy the gtf file" >> {log}
        cp {input.transcriptome_gtf} {output.gtf} >> {log}
        echo "completed copy of the gtf file" >> {log}
        '''

rule pureExonMultiRun:
    '''
    This is the rule to perform the extraction of the pure exonic reads from the dataset.
    '''
    input:
        bam_file = rules.moveCellrangerIntronExonMultiRun.output.bam_file
    output:
        exon_bam = config["out_location"] + "PureIntronExon/merged/{sample_name}/possorted_genome_bam_exon_filterd.bam",
        exon_feature = config["out_location"] + "PureIntronExon/merged/{sample_name}/feature_test_exon_fullBarcodes.txt.gz",
        exon_count = config["out_location"] + "PureIntronExon/merged/{sample_name}/summarized_feature_test_exon_fullBarcodes.txt"
    conda: config["env_bioinfo"]
    log:
        'logs/merged/{sample_name}/pureExonMultiRun.log'
    benchmark:
        'benchmarks/merged/{sample_name}/pureExonMultiRun.txt'
    resources:
        # mem_mb = lambda wildcards, attempt: config['filtering_RAM'] * attempt,
        # mem_mb = 16000,
        mem_mb = config["set-resources"]["pureExon"]["mem_mb"],
        cpus = config["set-threads"]["pureExon"]
    threads:
        config["set-threads"]["pureExon"]
    params:
        outs_folder = directory(config["out_location"] + "PureIntronExon/merged/{sample_name}")
    shell:
        '''
        echo "start pure exon filtering for <{wildcards.sample_name}>" >> {log}
         
        # Filter BAM file for records with 'xf:i:25' 
        samtools view {input.bam_file} | LC_ALL=C grep "xf:i:25" | grep "RE:A:E" > {params.outs_folder}/body_filtered_exon_sam
        echo "exon filtering for <{wildcards.sample_name}> done" >> {log}
        
        # Extract the BAM header and write to header_filted_sam
        samtools view -H {input.bam_file} > {params.outs_folder}/header_filted_exon_sam
        echo "header exon filtering for <{wildcards.sample_name}> done" >> {log}
        
        # Combine the header and extracted records
        cat {params.outs_folder}/header_filted_exon_sam {params.outs_folder}/body_filtered_exon_sam > {params.outs_folder}/possorted_genome_bam_exon_filterd.sam
        echo "generated sam for <{wildcards.sample_name}>" >> {log}
        
        # Convert SAM to BAM
        samtools view -b {params.outs_folder}/possorted_genome_bam_exon_filterd.sam > {output.exon_bam}
        # remove the files that are not needed any more
        rm -fr {params.outs_folder}/body_filtered_exon_sam {params.outs_folder}/possorted_genome_bam_exon_filterd.sam {params.outs_folder}/header_filted_exon_sam
        echo "SAM converted to BAM for <{wildcards.sample_name}>" >> {log}

        # extract the whole feature per barcode
        samtools view {params.outs_folder}/possorted_genome_bam_exon_filterd.bam | cut -f 12- | gzip > {params.outs_folder}/feature_test_exon_fullBarcodes.txt.gz
        echo "extract the pure exonic features per barcode for <{wildcards.sample_name}>" >> {log}
        
        # count the feature per barcode
        zcat {params.outs_folder}/feature_test_exon_fullBarcodes.txt.gz | \
        awk -F'\t' '
        {{
            gene = ""; barcode = "";
            for (i = 1; i <= NF; i++) {{
                if ($i ~ /^GX:Z:/) gene = $i;
                if ($i ~ /^CB:Z:/) barcode = $i;
            }}
            
            if (gene != "" && barcode != "") {{
                print gene, barcode;
            }}
        }}' | \
        sort | uniq -c | sort -k1,1nr -k2,2 -k3,3 > {params.outs_folder}/summarized_feature_test_exon_fullBarcodes.txt
        echo "count the pure exonic features per barcode for <{wildcards.sample_name}>" >> {log}
        '''

rule pureIntronMultiRun:
    '''
    This is the rule to perform the extraction of the pure intronic reads from the dataset.
    '''
    input:
        bam_file = rules.moveCellrangerIntronExonMultiRun.output.bam_file
    output:
        intron_bam = config["out_location"] + "PureIntronExon/merged/{sample_name}/possorted_genome_bam_intron_filterd.bam",
        intron_feature = config["out_location"] + "PureIntronExon/merged/{sample_name}/feature_test_intron_fullBarcodes.txt.gz",
        intron_count = config["out_location"] + "PureIntronExon/merged/{sample_name}/summarized_feature_test_intron_fullBarcodes.txt"
    conda: config['env_bioinfo']
    log:
        'logs/merged/{sample_name}/pureIntronMultiRun.log'
    benchmark:
        'benchmarks/merged/{sample_name}/pureIntronMultiRun.txt'
    resources:
        # mem_mb = lambda wildcards, attempt: config['filtering_RAM'] * attempt,
        # mem_mb = 16000,
        mem_mb = config["set-resources"]["pureIntron"]["mem_mb"],
        cpus = config["set-threads"]["pureIntron"]
    threads:
        config["set-threads"]["pureIntron"]
    params:
        outs_folder = directory(config["out_location"] + "PureIntronExon/merged/{sample_name}")
    shell:
        '''
        echo "start pure intron filtering for <{wildcards.sample_name}>" >> {log}
    
        # Filter BAM file for records with 'xf:i:25' 
        samtools view {input.bam_file} | LC_ALL=C grep "xf:i:25" | grep "RE:A:N" > {params.outs_folder}/body_filtered_intron_sam
        echo "intron filtering for <{wildcards.sample_name}> done" >> {log}
        
        # Extract the BAM header and write to header_filted_sam
        samtools view -H {input.bam_file} > {params.outs_folder}/header_filted_intron_sam
        echo "header intron filtering for <{wildcards.sample_name}> done" >> {log}
        
        # Combine the header and extracted records
        cat {params.outs_folder}/header_filted_intron_sam {params.outs_folder}/body_filtered_intron_sam > {params.outs_folder}/possorted_genome_bam_intron_filterd.sam
        echo "generated sam for <{wildcards.sample_name}>" >> {log}
        
        # Convert SAM to BAM
        samtools view -b {params.outs_folder}/possorted_genome_bam_intron_filterd.sam > {output.intron_bam}
    
        # remove the files that are not needed any more
        rm -fr {params.outs_folder}/body_filtered_intron_sam {params.outs_folder}/possorted_genome_bam_intron_filterd.sam {params.outs_folder}/header_filted_intron_sam
        echo "SAM converted to BAM for <{wildcards.sample_name}>" >> {log}

        # extract the whole feature per barcode
        samtools view {params.outs_folder}/possorted_genome_bam_intron_filterd.bam | cut -f 12- | gzip > {params.outs_folder}/feature_test_intron_fullBarcodes.txt.gz
        echo "extract the pure intronic features per barcode for <{wildcards.sample_name}>" >> {log}
        
        # count the feature per barcode
        zcat {params.outs_folder}/feature_test_intron_fullBarcodes.txt.gz | \
        awk -F'\t' '
        {{
            gene = ""; barcode = "";
            for (i = 1; i <= NF; i++) {{
                if ($i ~ /^GX:Z:/) gene = $i;
                if ($i ~ /^CB:Z:/) barcode = $i;
            }}
            
            if (gene != "" && barcode != "") {{
                print gene, barcode;
            }}
        }}' | \
        sort | uniq -c | sort -k1,1nr -k2,2 -k3,3 > {params.outs_folder}/summarized_feature_test_intron_fullBarcodes.txt
        echo "count the pure intronic features per barcode for <{wildcards.sample_name}>" >> {log}
        '''

rule pureAllMultiRun:
    '''
    This is the rule to perform the extraction of the pure intronic and exonic reads from the dataset.
    '''
    input:
        bam_file = rules.moveCellrangerIntronExonMultiRun.output.bam_file
    output:
        all_bam = config["out_location"] + "PureIntronExon/merged/{sample_name}/possorted_genome_bam_all_filterd.bam",
        all_feature = config["out_location"] + "PureIntronExon/merged/{sample_name}/feature_test_all_fullBarcodes.txt.gz",
        all_count = config["out_location"] + "PureIntronExon/merged/{sample_name}/summarized_feature_test_all_fullBarcodes.txt"
    conda: config['env_bioinfo']
    log:
        'logs/merged/{sample_name}/pureAllMultiRun.log'
    benchmark:
        'benchmarks/merged/{sample_name}/pureAllMultiRun.txt'
    resources:
        # mem_mb = lambda wildcards, attempt: config['filtering_RAM'] * attempt,
        # mem_mb = 16000,
        mem_mb = config["set-resources"]["pureAll"]["mem_mb"],
        cpus = config["set-threads"]["pureAll"]
    threads:
        config["set-threads"]["pureAll"]
    params:
        outs_folder = directory(config["out_location"] + "PureIntronExon/merged/{sample_name}")
    shell:
        '''
        echo "start pure all filtering for <{wildcards.sample_name}>" >> {log}
        
        # Filter BAM file for records with 'xf:i:25' 
        samtools view {input.bam_file} | LC_ALL=C grep "xf:i:25" > {params.outs_folder}/body_filtered_all_sam
        echo "all filtering for <{wildcards.sample_name}> done" >> {log}
        
        # Extract the BAM header and write to header_filted_sam
        samtools view -H {input.bam_file} > {params.outs_folder}/header_filted_all_sam
        echo "header all filtering for <{wildcards.sample_name}> done" >> {log}
        
        # Combine the header and extracted records
        cat {params.outs_folder}/header_filted_all_sam {params.outs_folder}/body_filtered_all_sam > {params.outs_folder}/possorted_genome_bam_all_filterd.sam
        echo "generated sam for <{wildcards.sample_name}>" >> {log}
        
        # Convert SAM to BAM
        samtools view -b {params.outs_folder}/possorted_genome_bam_all_filterd.sam > {output.all_bam}
    
        # remove the files that are not needed any more
        rm -fr {params.outs_folder}/body_filtered_all_sam {params.outs_folder}/possorted_genome_bam_all_filterd.sam {params.outs_folder}/header_filted_all_sam
        echo "SAM converted to BAM for <{wildcards.sample_name}>" >> {log}

        # extract the whole feature per barcode
        samtools view {params.outs_folder}/possorted_genome_bam_all_filterd.bam | cut -f 12- | gzip > {params.outs_folder}/feature_test_all_fullBarcodes.txt.gz
        echo "extract the pure intronic and exonic features per barcode for <{wildcards.sample_name}>" >> {log}
        
        # count the feature per barcode
        zcat {params.outs_folder}/feature_test_all_fullBarcodes.txt.gz | \
        awk -F'\t' '
        {{
            gene = ""; barcode = "";
            for (i = 1; i <= NF; i++) {{
                if ($i ~ /^GX:Z:/) gene = $i;
                if ($i ~ /^CB:Z:/) barcode = $i;
            }}
            
            if (gene != "" && barcode != "") {{
                print gene, barcode;
            }}
        }}' | \
        sort | uniq -c | sort -k1,1nr -k2,2 -k3,3 > {params.outs_folder}/summarized_feature_test_all_fullBarcodes.txt
        echo "count the pure all features per barcode for <{wildcards.sample_name}>" >> {log}
        '''
