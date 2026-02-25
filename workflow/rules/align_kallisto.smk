########################################
# Kallisto transcriptome quantification
########################################

TRANSCRIPTOME = config["reference"]["transcriptome_fasta"]
THREADS = config["alignment"]["threads"]
LAYOUT = config["sequencing"]["layout"]

INDEX_PATH = f"{OUTDIR}/{PROJECT}/trans/index/kallisto_index.idx"

########################################
# Build kallisto index 
########################################

rule kallisto_index:
    input:
        TRANSCRIPTOME
    output:
        INDEX_PATH
    log:
        f"{OUTDIR}/{PROJECT}/trans/logs/kallisto_index.log"
    threads: 1
    shell:
        """
        mkdir -p {OUTDIR}/{PROJECT}/trans/index
        mkdir -p {OUTDIR}/{PROJECT}/trans/logs

        kallisto index \
        -k 15 \
            -i {output} \
            {input} \
            > {log} 2>&1
        """


########################################
# Quantify samples
########################################

if LAYOUT == "paired":

    rule kallisto_quant:
        input:
            index = rules.kallisto_index.output,
            r1 = f"{OUTDIR}/{PROJECT}/trim/{{sample}}_R1.trimmed.fastq.gz",
            r2 = f"{OUTDIR}/{PROJECT}/trim/{{sample}}_R2.trimmed.fastq.gz"
        output:
            f"{OUTDIR}/{PROJECT}/trans/kallisto/{{sample}}/abundance.tsv"
        log:
            f"{OUTDIR}/{PROJECT}/trans/logs/{{sample}}.kallisto.log"
        threads: THREADS
        shell:
            """
            mkdir -p {OUTDIR}/{PROJECT}/trans/kallisto/{wildcards.sample}
            mkdir -p {OUTDIR}/{PROJECT}/trans/logs

            kallisto quant \
                -i {input.index} \
                -o {OUTDIR}/{PROJECT}/trans/kallisto/{wildcards.sample} \
                -t {threads} \
                {input.r1} {input.r2} \
                > {log} 2>&1
            """

else:

    rule kallisto_quant:
        input:
            index = rules.kallisto_index.output,
            r1 = f"{OUTDIR}/{PROJECT}/trim/{{sample}}_R1.trimmed.fastq.gz"
        output:
            f"{OUTDIR}/{PROJECT}/trans/kallisto/{{sample}}/abundance.tsv"
        log:
            f"{OUTDIR}/{PROJECT}/trans/logs/{{sample}}.kallisto.log"
        threads: THREADS
        shell:
            """
            mkdir -p {OUTDIR}/{PROJECT}/trans/kallisto/{wildcards.sample}
            mkdir -p {OUTDIR}/{PROJECT}/trans/logs

            kallisto quant \
                -i {input.index} \
                -o {OUTDIR}/{PROJECT}/trans/kallisto/{wildcards.sample} \
                -t {threads} \
                {input.r1} \
                > {log} 2>&1
            """