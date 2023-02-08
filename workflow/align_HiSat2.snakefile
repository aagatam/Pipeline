import pandas as pd
configfile: "configs/config.yaml"

samples = pd.read_csv(config["METAFILE"], sep = ',', header = 0)['Group']

indexes = list(range(1, 9))
end = config["END"]
compression = config['COMPRESSION_TYPE']
compressed = config["COMPRESSED"]
index_path = config["FINALOUTPUT"] + "/" + config["PROJECT"] + "/genome"
final_path = config["FINALOUTPUT"] + "/" + config["PROJECT"] + "/genome"
input_path = config["INPUTPATH"]


rule end:
    input:
        GTF = expand(final_path + "/countFile/{sample}/{sample}.gtf",sample=samples),
        report = final_path + "/report_align_count.html",
        gene_mat = final_path + "/Hisat_results/gene_count_matrix.csv",
        trans_mat = final_path + "/Hisat_results/transcript_count_matrix.csv"

rule indexGenome:
    input:
        genome = config["GENOME"]
    output:
        indexes = expand(final_path + "/indexes/index.{index}.ht2", index = indexes),

    params:
        index = index_path + "/indexes/index"
    priority: 10
    shell:
        "hisat2-build -p {config[NCORE]} {input.genome} {params.index}"
if compressed == "yes":
    rule uncompress:
        input:
            read = input_path + "/{sample}.fastq.dsrc"
        output:
            uncompress = temp(final_path + "/uncompressed/{sample}.out.fastq"),
        shell:
            "dsrc d -t{config[NCORE]} -s {input.read} >>{output.uncompress} "
else:
    rule uncompress:
        input:
            read = input_path + "/{sample}.fastq"
        output:
            uncompress =  temp(final_path + "/uncompressed/{sample}.out.fastq")
        shell:
            "ln -s {input.read} {output.uncompress}"


if end == "pair":
    rule alignment:
        input:
            index = expand(index_path + "/indexes/index.{index}.ht2", index = indexes),
            forward = input_path + "{sample}.out.fastq",
            reverse = input_path + "{sample}.out.fastq"
        output:
            sam = temp(final_path + "/samFile/{sample}.sam"),
            bam = temp(final_path + "/bamFile/{sample}.bam"),
            unaligned = temp(final_path + "/unaligned/{sample}.sam"),
            unaligned_bam = final_path + "/unaligned/{sample}.bam"
        params:
            index = index_path + "/indexes/index"
        benchmark:
            final_path + "/benchmarks/{sample}.hisat2.benchmark.txt"
        run:
            shell("hisat2 -q -p {config[NCORE]}  -x {params.index} -1 {input.forward} -2 {input.reverse} -S {output.sam} --un {output.unaligned}")
            shell("samtools view -@ {config[NCORE]} -b -S {output.sam} > {output.bam}")
            shell("samtools view -@ {config[NCORE]} -b -S {output.unaligned} > {output.unaligned_bam}")
else:
    rule alignment:
        input:
            index = expand(index_path + "/indexes/index.{index}.ht2", index = indexes),
            uncompress = temp(final_path + "/uncompressed/{sample}.out.fastq")
        output:
            sam = temp(final_path + "/samFile/{sample}.sam"),
            bam = temp(final_path + "/bamFile/{sample}.bam"),
            unaligned = temp(final_path + "/unaligned/{sample}.sam"),
            log = final_path + "/countFile/{sample}/{sample}_log.txt"
        params:
            index = index_path + "/indexes/index"
        benchmark:
            final_path + "/benchmarks/{sample}.hisat2.benchmark.txt"
        run:
            shell("hisat2 -q -p {config[NCORE]} -x {params.index} --summary-file {output.log} -U {input.uncompress} -S {output.sam} --un {output.unaligned}")
            shell("samtools view -@ {config[NCORE]} -b -S {output.sam} > {output.bam}")


rule sortBAM:
    input:
        bam = final_path + "/bamFile/{sample}.bam"
    output:
        sort = final_path + "/bamFileSort/{sample}.sort.bam"
    shell:
        "samtools sort -@ {config[NCORE]} {input.bam} -o {output.sort}"

rule alignmentQC:
    input:
        sort = final_path + "/bamFileSort/{sample}.sort.bam"
    output:
        bamqc = directory(final_path + "/alignmentQC/{sample}_BAMqc")
    shell:
        "qualimap bamqc -bam {input.sort} --java-mem-size=4G -outdir {output.bamqc}"

rule featureCount:
    input:
        sort = final_path + "/bamFileSort/{sample}.sort.bam",
        annotation = config["ANNOTATION"]
    output:
        GENE = final_path + "/countFile/{sample}/{sample}_gene.out",
        GTF = final_path + "/countFile/{sample}/{sample}.gtf"
    run:
        shell("nice stringtie {input.sort} -G {input.annotation} -o {output.GTF} -A {output.GENE} -p 8 -e -B")

rule prepDE:
    input:
        GTF = expand(final_path + "/countFile/{sample}/{sample}.gtf",sample=samples)
    params:
        GTF = final_path + "/countFile"
    output:
        gene_mat = final_path + "/Hisat_results/gene_count_matrix.csv",
        trans_mat = final_path + "/Hisat_results/transcript_count_matrix.csv"
    shell:
        "python ../scripts/prepDE.pr -i {params.GTF} -g {output.gene_mat} -t {output.trans_mat} "

rule summaryReport:
    input:
        bamqc = expand(final_path + "/alignmentQC/{sample}_BAMqc", sample = samples),
        count_summary = expand(final_path + "/countFile/{sample}/{sample}_log.txt", sample = samples)
    output:
        report = final_path + "/report_align_count.html"
    shell:
        "multiqc -f {input.bamqc} {input.count_summary} --filename {output.report}"
