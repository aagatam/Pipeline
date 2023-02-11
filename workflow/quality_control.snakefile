import pandas as pd
configfile: "configs/config.yaml"
samples = pd.read_csv(config["METAFILE"], sep = ',', header = 0)['Group']
end = config["END"]
compressed = config["COMPRESSED"]
final_path = config["FINALOUTPUT"] + "/" + config["PROJECT"] + "/genome"
input_path = config["INPUTPATH"]

rule end:
    input:
        report = final_path + "/fastqc/report_quality_control.html"

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
                read_trim = intermediate_path + "/{sample}_trimmed.fq.gz"
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

    rule qualityControl:
        input:
            forward = input_path + "/{sample}_R1.fastq.gz",
            reverse = input_path + "/{sample}_R2.fastq.gz"
        output:
            fastqc_forward = final_path + "/fastqc/{sample}_R1_fastqc.html",
            fastqc_reverse = final_path + "/fastqc/{sample}_R2_fastqc.html"
        params:
            outputpath = final_path + "/fastqc"
        shell:
            "fastqc -t $(({config[NCORE]}+0)) -o {params.outputpath} {input.forward} && "
            "fastqc -t $(({config[NCORE]}+0)) -o {params.outputpath} {input.reverse}"

    rule summaryReport:
        input:
            fastqc_forward = expand(final_path + "/fastqc/{sample}_R1_fastqc.html", sample = samples),
            fastqc_reverse = expand(final_path + "/fastqc/{sample}_R2_fastqc.html", sample = samples)
        output:
            report = final_path + "/fastqc/report_quality_control.html"
        params:
            path = final_path + "/fastqc"
        shell:
            "multiqc {params.path} --filename {output.report}"

else:

    rule qualityControl:
        input:
            uncompress = temp(final_path + "/uncompressed/{sample}.out.fastq"),
        output:
            fastqc = final_path + "/fastqc/{sample}.out_fastqc.html"
        params:
            outputpath = final_path + "/fastqc"
        run:
            shell("fastqc -t {config[NCORE]} -o {params.outputpath} {input.uncompress}")

    rule summaryReport:
        input:
            fastqc = expand(final_path + "/fastqc/{sample}.out_fastqc.html", sample = samples)
        output:
            report = final_path + "/fastqc/report_quality_control.html"
        params:
            path = final_path + "/fastqc"
        shell:
            "multiqc {params.path} --filename {output.report}"
