samples = pd.read_csv(config["METAFILE"], sep=',', header=0)['Group']
input_path = config["INPUTPATH"]
final_path = f"{config['FINALOUTPUT']}/{config['PROJECT']}/genome"

rule end:
    input:
        report = final_path + "/fastqc/report_quality_control.html"

# ===== SINGLE RULE FOR QUALITY CONTROL (FASTQC) =====
rule qualityControl:
    input:
        # Dynamic input: Paired-end or single-end?
        files = (
            [
                final_path + "/uncompressed/{sample}_R1.out.fastq",
                final_path + "/uncompressed/{sample}_R2.out.fastq"
            ]
            if config.get("end", "pair") == "pair"
            else final_path + "/uncompressed/{sample}.out.fastq"
        )
    output:
        # Dynamic output: Paired-end or single-end?
        html_files = (
            [
                final_path + "/fastqc/{sample}_R1_fastqc.html",
                final_path + "/fastqc/{sample}_R2_fastqc.html"
            ]
            if config.get("end", "pair") == "pair"
            else final_path + "/fastqc/{sample}.out_fastqc.html"
        )
    params:
        outputpath = final_path + "/fastqc"
    shell:
        # Dynamic shell command: Run FastQC on 1 or 2 files
        "fastqc -t {config[NCORE]} -o {params.outputpath} {input.files}"

# ===== SINGLE RULE FOR SUMMARY REPORT (MULTIQC) =====
rule summaryReport:
    input:
        # Dynamic input: Collect all FastQC outputs
        fastqc_files = (
            # Paired-end: R1 and R2 for all samples
            expand(final_path + "/fastqc/{sample}_R1_fastqc.html", sample=samples) +
            expand(final_path + "/fastqc/{sample}_R2_fastqc.html", sample=samples)
            if config.get("end", "pair") == "pair"
            # Single-end: One file per sample
            else expand(final_path + "/fastqc/{sample}.out_fastqc.html", sample=samples)
        )
    output:
        report = final_path + "/fastqc/report_quality_control.html"
    params:
        path = final_path + "/fastqc"
    shell:
        "multiqc {params.path} --filename {output.report}"