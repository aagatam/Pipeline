import pandas as pd

configfile: "configs/config.yaml"

samples = pd.read_csv(config["METAFILE"], sep=',', header=0)['Group']
end = config["END"]
input_path = config["INPUTPATH"]
final_path = f"{config['FINALOUTPUT']}/{config['PROJECT']}/genome"

rule end:
    input:
        report = final_path + "/fastqc/report_quality_control.html"

if end == "pair":

    rule qualityControl:
        input:
            forward = final_path + "/uncompressed/{sample}_R1.out.fastq",
            reverse = final_path + "/uncompressed/{sample}_R2.out.fastq"
        output:
            fastqc_forward = final_path + "/fastqc/{sample}_R1_fastqc.html",
            fastqc_reverse = final_path + "/fastqc/{sample}_R2_fastqc.html"
        params:
            outputpath = final_path + "/fastqc"
        shell:
            "fastqc -t {config[NCORE]} -o {params.outputpath} {input.forward} && "
            "fastqc -t {config[NCORE]} -o {params.outputpath} {input.reverse}"

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
            uncompress = final_path + "/uncompressed/{sample}.out.fastq"
        output:
            fastqc = final_path + "/fastqc/{sample}.out_fastqc.html"
        params:
            outputpath = final_path + "/fastqc"
        shell:
            "fastqc -t {config[NCORE]} -o {params.outputpath} {input.uncompress}"

    rule summaryReport:
        input:
            fastqc = expand(final_path + "/fastqc/{sample}.out_fastqc.html", sample = samples)
        output:
            report = final_path + "/fastqc/report_quality_control.html"
        params:
            path = final_path + "/fastqc"
        shell:
            "multiqc {params.path} --filename {output.report}"
