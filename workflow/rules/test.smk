rule testRenv:
    '''
    A simple dummy rule to test if the Renv environment activates correctly.
    '''
    output:
        "results/test/renv_test_passed.txt"
    conda: config["env_testR"]
    script:
        "../scripts/00_test_renv.R"