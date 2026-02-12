# --- Global Variables ---

# Determine the filename tag based on the config flag. this can be used to tailor the name of the output files
if config.get("run_cellbender", False):
    CELLBENDER_TAG = "cellbender"
else:
    CELLBENDER_TAG = "default"

# --- Helper functions ---

# Helper function: allow to switch the input to the runSamplePreprocessing in case we run Cellbender
# config.get It looks for the key; if it doesn't find it, it uses the fallback value (the second argument).
# If the key exists: It uses the value from the config.
def get_preprocessing_input(wildcards):
    if config.get("run_cellbender", False):
        return rules.runCellbender.output.filter_h5
    else:
        return rules.runCellrangerMultiRun.output.folder

