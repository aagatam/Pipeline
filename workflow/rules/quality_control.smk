########################################
# Quality Control with FastQC + MultiQC
########################################

LAYOUT = config["sequencing"]["layout"]

rule fastqc:
    input:
        r1 = f"{OUTDIR}/{PROJECT}/trim/{{sample}}_R1.trimmed.fastq.gz",
        r2 = f"{OUTDIR}/{PROJECT}/trim/{{sample}}_R2.trimmed.fastq.gz" if LAYOUT == "paired" else None
    output:
        r1_html = f"{OUTDIR}/{PROJECT}/qc/fastqc/{{sample}}_R1_fastqc.html",
        r2_html = f"{OUTDIR}/{PROJECT}/qc/fastqc/{{sample}}_R2_fastqc.html" if LAYOUT == "paired" else None
    log:
        f"{OUTDIR}/{PROJECT}/logs/fastqc/{{sample}}.log"
    threads: 2
    shell:
        """
        mkdir -p {OUTDIR}/{PROJECT}/qc/fastqc
        mkdir -p {OUTDIR}/{PROJECT}/logs/fastqc

        if [ "{LAYOUT}" = "paired" ]; then
            fastqc {input.r1} {input.r2} \
                --outdir {OUTDIR}/{PROJECT}/qc/fastqc \
                > {log} 2>&1
        else
            fastqc {input.r1} \
                --outdir {OUTDIR}/{PROJECT}/qc/fastqc \
                > {log} 2>&1
        fi
        """
        

rule multiqc:
    input:
        expand(
            f"{OUTDIR}/{PROJECT}/qc/fastqc/{{sample}}_R1_fastqc.html",
            sample=SAMPLES
        )
    output:
        f"{OUTDIR}/{PROJECT}/qc/multiqc_report.html"
    log:
        f"{OUTDIR}/{PROJECT}/logs/multiqc.log"
    shell:
        """
        mkdir -p {OUTDIR}/{PROJECT}/qc
        mkdir -p {OUTDIR}/{PROJECT}/logs

        multiqc {OUTDIR}/{PROJECT}/qc/fastqc \
            -o {OUTDIR}/{PROJECT}/qc \
            > {log} 2>&1
        """