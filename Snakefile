# Snakefile

import pandas as pd

# --- Load config and variables ---
configfile: "configs/config.yaml"
samples = pd.read_csv(config["METAFILE"], sep=',', header=0)['Group']
end = config["end"]
compression = config["COMPRESSION_TYPE"]
final_path = f"{config['FINALOUTPUT']}/{config['PROJECT']}/genome"
input_path = config["INPUTPATH"]

# Include rule files
#include: "workflow/rules/trim.snakefile"
include: "workflow/rules/uncompress.snakefile"
include: "workflow/rules/quality_control.snakefile"
#include: "workflow/rules/align_HiSat2.snakefile"
#include: "workflow/rules/align_kallisto.snakefile"
#include: "workflow/rules/spladder_run.snakefile"
#include: "workflow/rules/bisbee_run.snakefile"
#include: "workflow/rules/interproscan_run.snakefile"

sys.path.append("workflow/scripts")
from uncompress_helpers import *

rule all:
    input:
        report = final_path + "/fastqc/report_quality_control.html"
