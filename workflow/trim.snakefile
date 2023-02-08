import pandas as pd

configfile: "configs/config.yaml"
samples = pd.read_csv(config["METAFILE"], sep = ',', header = 0)['Group']
end = config["END"]
input_path = config["INPUTPATH"]
compression = config['COMPRESSION_TYPE']
intermediate_path = config["FINALOUTPUT"] + "/" + config["PROJECT"] + "/trim"
final_path = config["FINALOUTPUT"] + "/" + config["PROJECT"] + "/genome"

def trimFiles(wildcards):
    if (end == "pair"):
        forward_trim = expand(intermediate_path + "/{sample}_R1.fq.dsrc", sample = samples)
        return forward_trim
    else:
        read_trim = expand(intermediate_path + "/{sample}.out_trimmed.fq.dsrc", sample = samples)
        return read_trim

rule all:
    input:
        report = final_path + "/trim/fastqc_after_trimming/report_quality_control_after_trimming.html",
        # compressed = expand(final_path + "/trim/compressed/{sample}.fastq.dsrc",sample=samples)

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
                shell("gunzip -c {input.forward} > {output.uncompress1}")
                shell("gunzip -c {input.reverse} > {output.uncompress1}")
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
                "gunzip -c {input.read} > {output.uncompress}"
    else:
        rule uncompress:
            input:
                read = input_path + "/{sample}.fastq"
            output:
                uncompress =  temp(final_path + "/uncompressed/{sample}.out.fastq")
            shell:
                "ln -s {input.read} {output.uncompress}"

if end == "pair":
    rule trim:
        input:
            forward = temp(final_path + "/uncompressed/{sample}_R1.out.fastq"),
            reverse = temp(final_path + "/uncompressed/{sample}_R2.out.fastq"),
        output:
            read_trim_forward = intermediate_path + "/{sample}_R1.fq.gz",
            read_trim_reverse = intermediate_path + "/{sample}_R2.fq.gz"
        params:
            outputpath = intermediate_path
        shell:
            "trim_galore --fastqc -j 4 --gzip --paired --basename {wildcards.sample} -o {params.outputpath} {input.forward} {input.reverse}"

else:
    rule trim:
        input:
            uncompress = temp(final_path + "/uncompressed/{sample}.out.fastq")
        output:
            read_trim = temp(intermediate_path + "/{sample}.out_trimmed.fq")
        params:
            outputpath = intermediate_path
        shell:
            "trim_galore --fastqc -j 32 --gzip --basename {wildcards.sample} -o {params.outputpath} {input.uncompress}"

rule summaryReport:
    input:
        trimFiles
    output:
        report = final_path + "/trim/fastqc_after_trimming/report_quality_control_after_trimming.html"
    params:
        path = intermediate_path
    shell:
        "multiqc {params.path} --filename {output.report}"
