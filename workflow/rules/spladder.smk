########################################
# SplAdder - alternative splicing
########################################

import os

ANNOTATION = config["reference"]["annotation_gtf"]

BAM_FILES = expand(
    f"{OUTDIR}/{PROJECT}/genome/bam/{{sample}}.sorted.bam",
    sample=SAMPLES
)

BAI_FILES = expand(
    f"{OUTDIR}/{PROJECT}/genome/bam/{{sample}}.sorted.bam.bai",
    sample=SAMPLES
)

SPLADDER_OUT = f"{OUTDIR}/{PROJECT}/spladder"

########################################
# Index BAMs
########################################

rule index_bam:
    input:
        f"{OUTDIR}/{PROJECT}/genome/bam/{{sample}}.sorted.bam"
    output:
        f"{OUTDIR}/{PROJECT}/genome/bam/{{sample}}.sorted.bam.bai"
    shell:
        """
        samtools index {input}
        """

########################################
# SplAdder build
########################################

rule spladder_build:
    input:
        bam = BAM_FILES,
        gtf = ANNOTATION
    output:
         directory(SPLADDER_OUT)
    threads: 8
    params:
        bam_list = lambda wildcards, input: ",".join(input.bam)
    shell:
        """
        mkdir -p {output}

        spladder build \
            --bams {params.bam_list} \
            --annotation {input.gtf} \
            --outdir {output} \
            --event-types exon_skip,mutex_exons,alt_3prime,alt_5prime,mult_exon_skip \
            --parallel {threads}
        """