import pandas as pd

configfile: "configs/config.yaml"
samples = pd.read_csv(config["METAFILE"], sep = ',', header = 0)['Group']
indexes = list(range(1, 9))
end = config["END"]
compression = config['COMPRESSION_TYPE']
trimmed=config['TRIMMED']
index_path = config["FINALOUTPUT"] + "/" + config["PROJECT"] + "/genome"
final_path = config["FINALOUTPUT"] + "/" + config["PROJECT"] + "/genome"
input_path = config["INPUTPATH"]
intermediate_path = config["FINALOUTPUT"] + "/" + config["PROJECT"] + "/trim"


rule end:
    input:
        GTF = expand(final_path + "/countFile/{sample}/{sample}.gtf",sample=samples),
        report = final_path + "/report_align_count.html",
        gene_mat = final_path + "/gene_count_matrix.csv",
        trans_mat = final_path + "/transcript_count_matrix.csv"

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

if end == "pair":
    if compression == "dsrc":
        rule uncompress:
            input:
                forward = input_path + "{sample}_R1.fastq.dsrc",
                reverse = input_path + "/{sample}_R2.fastq.dsrc"
            output:
                uncompress1 = temp(final_path + "/uncompressed/{sample}_R1.out.fastq"),
                uncompress2 = temp(final_path + "/uncompressed/{sample}_R2.out.fastq")
            run:
                shell("dsrc d -t{config[NCORE]} -s {input.forward} >>{output.uncompress1} ")
                shell("dsrc d -t{config[NCORE]} -s {input.reverse} >>{output.uncompress2} ")
    elif compression == 'gz':
        rule uncompress:
            input:
                forward = input_path + "{sample}_R1.fastq.gz",
                reverse = input_path + "/{sample}_R2.fastq.gz"
            output:
                uncompress1 = temp(final_path + "/uncompressed/{sample}_R1.out.fastq"),
                uncompress2 = temp(final_path + "/uncompressed/{sample}_R2.out.fastq")
            run:
                shell("pigz -d -k -c -p{config[NCORE]} {input.forward} > {output.uncompress1}")
                shell("pigz -d -k -c -p{config[NCORE]} {input.reverse} > {output.uncompress1}")
    elif trimmed == 'yes':
        rule uncompress:
            input:
                read_trim_forward = intermediate_path + "/{sample}_val_1.fq.gz",
                read_trim_reverse = intermediate_path + "/{sample}_val_2.fq.gz"
            output:
                uncompress1 = temp(final_path + "/uncompressed/{sample}_R1.out.fastq"),
                uncompress2 = temp(final_path + "/uncompressed/{sample}_R2.out.fastq")
            run:
                shell("pigz -d -k -c -p{config[NCORE]} {input.read_trim_forward} > {output.uncompress1}")
                shell("pigz -d -k -c -p{config[NCORE]} {input.read_trim_reverse} > {output.uncompress1}")

    else:
        rule uncompress:
            input:
                forward = input_path + "/{sample}_R1.fastq",
                reverse = input_path + "/{sample}_R2.fastq"
            output:
                uncompress1 =  temp(final_path + "/uncompressed/{sample}_R1.out.fastq"),
                uncompress2 =  temp(final_path + "/uncompressed/{sample}_R2.out.fastq")
            run:
                shell("ln -s {input.forward} {output.uncompress1}")
                shell("ln -s {input.reverse} {output.uncompress2}")
else:
    if compression == "dsrc":
        rule uncompress:
            input:
                read = input_path + "/{sample}.fastq.dsrc"
            output:
                uncompress = temp(final_path + "/uncompressed/{sample}.out.fastq"),
            shell:
                "dsrc d -t{config[NCORE]} -s {input.read} >>{output.uncompress} "
    elif compression == 'gz':
        rule uncompress:
            input:
                read = input_path + "/{sample}.fastq.gz"
            output:
                uncompress =  temp(final_path + "/uncompressed/{sample}.out.fastq")
            shell:
                "pigz -d -k -c -p{config[NCORE]} {input.read} > {output.uncompress}"
    elif trimmed == 'yes':
        rule uncompress:
            input:
                read_trim = intermediate_path + "/{sample}_trimmed.fq"
            output:
                uncompress =  temp(final_path + "/uncompressed/{sample}.out.fastq")
            shell:
                "pigz -d -k -c -p{config[NCORE]} {input.read_trim} > {output.uncompress}"
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
            forward = temp(final_path + "/uncompressed/{sample}_R1.out.fastq"),
            reverse = temp(final_path + "/uncompressed/{sample}_R2.out.fastq"),
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
        gene_mat = final_path + "/gene_count_matrix.csv",
        trans_mat = final_path + "/transcript_count_matrix.csv"
    shell:
        "python ./scripts/prepDE.py -i {params.GTF} -g {output.gene_mat} -t {output.trans_mat} "

rule summaryReport:
    input:
        bamqc = expand(final_path + "/alignmentQC/{sample}_BAMqc", sample = samples),
        count_summary = expand(final_path + "/countFile/{sample}/{sample}_log.txt", sample = samples)
    output:
        report = final_path + "/report_align_count.html"
    shell:
        "multiqc -f {input.bamqc} {input.count_summary} --filename {output.report}"
