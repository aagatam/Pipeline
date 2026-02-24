########################################
# Trimming reads
########################################

TRIM_THREADS = config["trimming"]["threads"]
LAYOUT = config["sequencing"]["layout"]

rule trim:
    input:
        r1 = f"{OUTDIR}/{PROJECT}/reads/{{sample}}_R1.fastq",
        r2 = f"{OUTDIR}/{PROJECT}/reads/{{sample}}_R2.fastq" if LAYOUT == "paired" else None
    output:
        r1 = f"{OUTDIR}/{PROJECT}/trim/{{sample}}_R1.trimmed.fastq.gz",
        r2 = f"{OUTDIR}/{PROJECT}/trim/{{sample}}_R2.trimmed.fastq.gz" if LAYOUT == "paired" else None
    log:
        f"{OUTDIR}/{PROJECT}/logs/trim/{{sample}}.log"
    threads: TRIM_THREADS
    shell:
        """
        mkdir -p {OUTDIR}/{PROJECT}/trim
        mkdir -p {OUTDIR}/{PROJECT}/logs/trim

        if [ "{LAYOUT}" = "paired" ]; then
            fastp \
                -i {input.r1} \
                -I {input.r2} \
                -o {output.r1} \
                -O {output.r2} \
                --thread {threads} \
                > {log} 2>&1
        else
            fastp \
                -i {input.r1} \
                -o {output.r1} \
                --thread {threads} \
                > {log} 2>&1
        fi
        """