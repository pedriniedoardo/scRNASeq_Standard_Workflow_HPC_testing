# make sure to ask for the gpu:
# srun -p cuda --gpus v100:1 --cpus-per-gpu=1 --mem=4GB --pty bash
. /home/pedrini.edoardo/miniconda3/bin/activate;
conda activate env_cellbender_pip

# cellbender -v

cellbender remove-background \
     --input /idle/ric.cosr/ric.cosr/pedrini.edoardo/sc-rna-seq-cosr-standard-workflow-test2/workflow/test/cellbender/example/tiny_raw_feature_bc_matrix.h5ad \
     --output /idle/ric.cosr/ric.cosr/pedrini.edoardo/sc-rna-seq-cosr-standard-workflow-test2/workflow/test/cellbender/example/out_GPU/tiny_output.h5 \
     --expected-cells 500 \
     --total-droplets-included 2000 \
     --cuda