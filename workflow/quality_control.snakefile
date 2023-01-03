import pandas as pd
configfile: "configs/config.yaml"
samples = pd.read_csv(config["METAFILE"], sep = ';', header = 0)['Sample']
end = config["END"]
compressed = config["COMPRESSED"]
final_path = config["FINALOUTPUT"] + "/" + config["PROJECT"] + "/genome"
input_path = config["INPUTPATH"]

rule end:
    input:
        report = final_path + "/fastqc/report_quality_control.html"

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
            shell("fastqc -t 8 -o {params.outputpath} {input.uncompress}")

    rule summaryReport:
        input:
            fastqc = expand(final_path + "/fastqc/{sample}.out_fastqc.html", sample = samples)
        output:
            report = final_path + "/fastqc/report_quality_control.html"
        params:
            path = final_path + "/fastqc"
        shell:
            "multiqc {params.path} --filename {output.report}"
