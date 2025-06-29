indexes = list(range(1, 9))

rule hi_alignment_end:
    input:
        GTF = expand(final_path_genome + "/countFile/{sample}/{sample}.gtf",sample=samples),
        report = final_path_genome + "/report_align_count.html",
        gene_mat = final_path_genome + "/gene_count_matrix.csv",
        trans_mat = final_path_genome + "/transcript_count_matrix.csv"

rule indexGenome:
    input:
        genome = config["GENOME"]
    output:
        indexes = expand(final_path_genome + "/indexes/index.{index}.ht2", index = indexes),

    params:
        index = index_path_genome + "/indexes/index"
    priority: 10
    shell:
        "hisat2-build -p {config[NCORE]} {input.genome} {params.index}"

rule alignment_hisat:
    input:
        index = expand(index_path_genome + "/indexes/index.{index}.ht2", index = indexes),
        # Input from trimmed uncompressed files
        read1 = f"{final_path}/trimmed_uncompressed/{{sample}}_val_1.fastq",
        read2 = f"{final_path}/trimmed_uncompressed/{{sample}}_val_2.fastq",
        single = f"{final_path}/trimmed_uncompressed/{{sample}}_trimmed.fastq"
    output:
        sam = temp(final_path_genome + "/samFile/{sample}.sam"),
        bam = temp(final_path_genome + "/bamFile/{sample}.bam"),
        unaligned = temp(final_path_genome + "/unaligned/{sample}.sam"),
        unaligned_bam = final_path + "/unaligned/{sample}.bam",
        log = final_path_genome + "/countFile/{sample}/{sample}_log.txt"
    params:
        index = index_path_genome + "/indexes/index",
        is_paired = lambda wildcards: config.get("end", "pair") == "pair"
    benchmark:
        final_path_genome + "/benchmarks/{sample}.hisat2.benchmark.txt"
    shell:
        """
        if [ "{params.is_paired}" == "True" ]; then
            hisat2 -q -p {config[NCORE]} -x {params.index} -1 {input.read1} -2 {input.read2} -S {output.sam} --un {output.unaligned}
            samtools view -@ {config[NCORE]} -b -S {output.sam} > {output.bam}
            samtools view -@ {config[NCORE]} -b -S {output.unaligned} > {output.unaligned_bam}
        else
            hisat2 -q -p {config[NCORE]} -x {params.index} --summary-file {output.log} -U {input.single} -S {output.sam} --un {output.unaligned}
            samtools view -@ {config[NCORE]} -b -S {output.sam} > {output.bam}
        fi
        """

rule sortBAM:
    input:
        bam = final_path_genome + "/bamFile/{sample}.bam"
    output:
        sort = final_path_genome + "/bamFileSort/{sample}.sort.bam"
    shell:
        "samtools sort -@ {config[NCORE]} {input.bam} -o {output.sort}"

rule alignmentQC:
    input:
        sort = final_path_genome + "/bamFileSort/{sample}.sort.bam"
    output:
        bamqc = directory(final_path_genome + "/alignmentQC/{sample}_BAMqc")
    shell:
        "qualimap bamqc -bam {input.sort} --java-mem-size=4G -outdir {output.bamqc}"

rule featureCount:
    input:
        sort = final_path_genome + "/bamFileSort/{sample}.sort.bam",
        annotation = config["ANNOTATION"]
    output:
        GENE = final_path_genome + "/countFile/{sample}/{sample}_gene.out",
        GTF = final_path_genome + "/countFile/{sample}/{sample}.gtf"
    run:
        shell("stringtie {input.sort} -G {input.annotation} -o {output.GTF} -A {output.GENE} -p 8 -e -B")

rule prepDE:
    input:
        GTF = expand(final_path_genome + "/countFile/{sample}/{sample}.gtf",sample=samples)
    params:
        GTF = final_path_genome + "/countFile"
    output:
        gene_mat = final_path_genome + "/gene_count_matrix.csv",
        trans_mat = final_path_genome + "/transcript_count_matrix.csv"
    shell:
        "python ./scripts/prepDE.py -i {params.GTF} -g {output.gene_mat} -t {output.trans_mat} "

rule Report_hisatt:
    input:
        bamqc = expand(final_path_genome + "/alignmentQC/{sample}_BAMqc", sample = samples),
        count_summary = expand(final_path_genome + "/countFile/{sample}/{sample}_log.txt", sample = samples)
    output:
        report = final_path_genome + "/report_align_count.html"
    shell:
        "multiqc -f {input.bamqc} {input.count_summary} --filename {output.report}"
