# Snakefile

# Load main config
configfile: "configs/config.yaml"

# Include rule files
include: "workflow/trim.snakefile"
include: "workflow/quality_control.snakefile"
include: "workflow/align_HiSat2.snakefile"
include: "workflow/align_kallisto.snakefile"
include: "workflow/spladder_run.snakefile"
include: "workflow/bisbee_run.snakefile"
include: "workflow/interproscan_run.snakefile"

rule all:
    input:
        expand("FINALOUTPUT/PROJECT/genome/bisbee/{sample}.grepped.filtered.fasta", sample=SAMPLES),
        expand("FINALOUTPUT/PROJECT/genome/InterProScan/{sample}.tsv", sample=SAMPLES)

