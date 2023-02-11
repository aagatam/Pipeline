import pandas as pd
configfile: "configs/config.yaml"

samples = pd.read_csv(config["METAFILE"], sep = ',', header = 0)['Group']
end = config["END"]
trimmed=config['TRIMMED']
index_path = config["FINALOUTPUT"] + "/" + config["PROJECT"] + "/trans"
compression = config['COMPRESSION_TYPE']
final_path = config["FINALOUTPUT"] + "/" + config["PROJECT"] + "/trans"
input_path = config["INPUTPATH"]
intermediate_path = config["FINALOUTPUT"] + "/" + config["PROJECT"] + "/trim"


rule end:
    input:
        report = final_path + "/report_align_count.html",
        index = index_path + "/indexes/kallisto_index.idx"

rule indexTrans:
    input:
        trans = config["TRANS"]
    output:
        index = index_path + "/indexes/kallisto_index.idx"
    priority: 10
    shell:
        "kallisto index -i {output.index} {input.trans}"

if end == "pair":
    if compression == "dsrc":
        rule uncompress:
            input:
                forward = input_path + "{sample}_R1.fastq.dsrc",
                reverse = input_path + "/{sample}_R2.fastq.dsrc"
            output:
                uncompress1 = temp(final_path + "/uncompressed/{sample}_R1.out.fastq"),
                uncompress1 = temp(final_path + "/uncompressed/{sample}_R2.out.fastq")
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
                uncompress1 = temp(final_path + "/uncompressed/{sample}_R2.out.fastq")
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
                uncompress1 = temp(final_path + "/uncompressed/{sample}_R2.out.fastq")
            run:
                shell("pigz -d -k -c -p{config[NCORE]}  {input.read_trim_forward} > {output.uncompress1}")
                shell("pigz -d -k -c -p{config[NCORE]}  {input.read_trim_reverse} > {output.uncompress1}")
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
                read_trim = intermediate_path + "/{sample}._trimmed.fq.gz"
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
            index = index_path + "/indexes/kallisto_index.idx",
            uncompress1 = temp(final_path + "/uncompressed/{sample}_R1.out.fastq"),
            uncompress2 = temp(final_path + "/uncompressed/{sample}_R2.out.fastq")
        output:
            out = directory(final_path + "/kallisto/{sample}"),
            log = final_path + "/kallisto/{sample}_log.txt"
        benchmark:
            final_path + "/benchmarks/{sample}.kallisto.benchmark.txt"
        shell:
            "kallisto quant --bias --single-overhang --fusion -i {input.index} -o {output.out}  -t {config[NCORE]} {input.uncompress1} {input.uncompress2} &>{output.log}"

else:
    rule alignment:
        input:
            index = index_path + "/indexes/kallisto_index.idx",
            uncompress = temp(final_path + "/uncompressed/{sample}.out.fastq")
        output:
            out = directory(final_path + "/kallisto/{sample}"),
            log = final_path + "/kallisto/{sample}_log.txt"
        benchmark:
            final_path + "/benchmarks/{sample}.kallisto.benchmark.txt"
        shell:
            "kallisto quant --single --bias --single-overhang --fusion -i {input.index} -o {output.out}  -t {config[NCORE]} -l 130 -s 30 {input.uncompress} &>{output.log}"

rule summaryReport:
    input:
        count_summary = expand(final_path + "/kallisto/{sample}_log.txt", sample = samples)
    output:
        report = final_path + "/report_align_count.html"
    shell:
        "multiqc -f {input.count_summary} --filename {output.report}"
