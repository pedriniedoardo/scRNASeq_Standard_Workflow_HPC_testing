rule testRenv:
    '''
    A simple dummy rule to test if the Renv environment activates correctly.
    '''
    output:
        "results/test/renv_test_passed.txt"
    conda: config["env_testR"]
    script:
        "../scripts/00_test_renv.R"

rule testAnnotationCyteTypeR:
    '''
    This rule run the automatic annotation using CyteTypeR
    This is a test rule
    '''
    input:
        # rds = rules.integrateSamples.output.rds,
        # markers = rules.integrateSamples.output.markers,
        # metadata = config["meta_annotation"]
        rds = "/idle/ric.cosr/ric.cosr/pedrini.edoardo/sc-rna-seq-cosr-standard-workflow-test2/results/Seurat/object/integrated_obj_cellbender.rds",
        markers = "/idle/ric.cosr/ric.cosr/pedrini.edoardo/sc-rna-seq-cosr-standard-workflow-test2/results/Seurat/table/integrated_markers_cellbender.tsv",
        metadata = config["meta_annotation"]
    output:
        # rds = config["out_location"] + "Seurat/object/integrated_obj_" + CELLBENDER_TAG + "_CyteTypeR.rds",
        # vars_h5_path = config["out_location"] + "Seurat/object/vars_" + CELLBENDER_TAG + ".h5",
        # obs_duckdb_path = config["out_location"] + "Seurat/object/obs_" + CELLBENDER_TAG + ".duckdb",
        # query_filename = config["out_location"] + "Seurat/object/query_" + CELLBENDER_TAG + ".json"
        rds = "/idle/ric.cosr/ric.cosr/pedrini.edoardo/sc-rna-seq-cosr-standard-workflow-test2/workflow/test/cytetypeR/test_out/scobj_CyteTypeR.rds",
        vars_h5_path = "/idle/ric.cosr/ric.cosr/pedrini.edoardo/sc-rna-seq-cosr-standard-workflow-test2/workflow/test/cytetypeR/test_out/vars.h5",
        obs_duckdb_path =  "/idle/ric.cosr/ric.cosr/pedrini.edoardo/sc-rna-seq-cosr-standard-workflow-test2/workflow/test/cytetypeR/test_out/obs.duckdb",
        query_filename = "/idle/ric.cosr/ric.cosr/pedrini.edoardo/sc-rna-seq-cosr-standard-workflow-test2/workflow/test/cytetypeR/test_out/query.json"
    conda: config["env_annotationCyteTypeR"]
    log:
        'logs/Seurat/integrated/testAnnotationCyteTypeR.log'
    benchmark:
        'benchmarks/Seurat/integrated/testAnnotationCyteTypeR.txt'
    params:
        # id_org = config["ref"]["organism_id"]
        # token = config["token_CyteTypeR"],
        cov_markers = config["cov_markers_CyteTypeR"]
    script:
        "../scripts/04_annotation_CyteTypeR.R"
