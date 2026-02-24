# workflow/rules/quality_control.snakefile

rule qc_end:
    input:
        report = final_path + "/fastqc/report_quality_control.html"

# ===== SINGLE RULE FOR QUALITY CONTROL (FASTQC) =====
rule qualityControl:
    input:
        # Now takes input from trimmed files instead of uncompressed
        files = (
            [
                f"{intermediate_path}/{{sample}}_val_1.fq.gz",
                f"{intermediate_path}/{{sample}}_val_2.fq.gz"
            ]
            if config.get("end", "pair") == "pair"
            else f"{intermediate_path}/{{sample}}_trimmed.fq.gz"
        )
    output:
        # Updated output file names to reflect trimmed input
        html_files = (
            [
                final_path + "/fastqc/{sample}_val_1_fastqc.html",
                final_path + "/fastqc/{sample}_val_2_fastqc.html"
            ]
            if config.get("end", "pair") == "pair"
            else final_path + "/fastqc/{sample}_trimmed_fastqc.html"
        )
    params:
        outputpath = final_path + "/fastqc"
    shell:
        # Dynamic shell command: Run FastQC on 1 or 2 files
        "fastqc -t {config[NCORE]} -o {params.outputpath} {input.files}"

# ===== SINGLE RULE FOR SUMMARY REPORT (MULTIQC) =====
rule summaryReport:
    input:
        # Dynamic input: Collect all FastQC outputs from trimmed files
        fastqc_files = (
            # Paired-end: val_1 and val_2 for all samples
            expand(final_path + "/fastqc/{sample}_val_1_fastqc.html", sample=samples) +
            expand(final_path + "/fastqc/{sample}_val_2_fastqc.html", sample=samples)
            if config.get("end", "pair") == "pair"
            # Single-end: One trimmed file per sample
            else expand(final_path + "/fastqc/{sample}_trimmed_fastqc.html", sample=samples)
        )
    output:
        report = final_path + "/fastqc/report_quality_control.html"
    params:
        path = final_path + "/fastqc"
    shell:
        "multiqc {params.path} --filename {output.report}"