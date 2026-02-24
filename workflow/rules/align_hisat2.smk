########################################
# HISAT2 genome alignment
########################################

GENOME = config["reference"]["genome_fasta"]
THREADS = config["alignment"]["threads"]
LAYOUT = config["sequencing"]["layout"]

INDEX_PREFIX = f"{OUTDIR}/{PROJECT}/genome/index/hisat2_index"

########################################
# Build HISAT2 index (once)
########################################

rule hisat2_index:
    input:
        GENOME
    output:
        expand(INDEX_PREFIX + ".{i}.ht2", i=range(1,9))
    log:
        f"{OUTDIR}/{PROJECT}/genome/logs/hisat2_index.log"
    threads: THREADS
    shell:
        """
        mkdir -p {OUTDIR}/{PROJECT}/genome/index
        mkdir -p {OUTDIR}/{PROJECT}/genome/logs

        hisat2-build {input} {INDEX_PREFIX} \
            > {log} 2>&1
        """


########################################
# Align reads
########################################

rule hisat2_align:
    input:
        index = rules.hisat2_index.output,
        r1 = f"{OUTDIR}/{PROJECT}/trim/{{sample}}_R1.trimmed.fastq.gz",
        r2 = f"{OUTDIR}/{PROJECT}/trim/{{sample}}_R2.trimmed.fastq.gz" if LAYOUT == "paired" else None
    output:
        bam = f"{OUTDIR}/{PROJECT}/genome/bam/{{sample}}.sorted.bam"
    log:
        f"{OUTDIR}/{PROJECT}/genome/logs/{{sample}}.hisat2.log"
    threads: THREADS
    shell:
        """
        mkdir -p {OUTDIR}/{PROJECT}/genome/bam
        mkdir -p {OUTDIR}/{PROJECT}/genome/logs

        if [ "{LAYOUT}" = "paired" ]; then
            hisat2 \
                -x {INDEX_PREFIX} \
                -1 {input.r1} \
                -2 {input.r2} \
                -p {threads} 2>> {log} \
            | samtools view -bS - \
            | samtools sort -@ {threads} -o {output.bam}
        else
            hisat2 \
                -x {INDEX_PREFIX} \
                -U {input.r1} \
                -p {threads} 2>> {log} \
            | samtools view -bS - \
            | samtools sort -@ {threads} -o {output.bam}
        fi

        samtools index {output.bam}
        """