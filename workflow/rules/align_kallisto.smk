rule al_kallisto_end:
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

# Single rule for alignment that handles both paired and single-end
rule alignment:
    input:
        index = index_path + "/indexes/kallisto_index.idx",
        # Input from trimmed uncompressed files
        read1 = f"{final_path}/trimmed_uncompressed/{{sample}}_val_1.fastq",
        read2 = f"{final_path}/trimmed_uncompressed/{{sample}}_val_2.fastq",
        single = f"{final_path}/trimmed_uncompressed/{{sample}}_trimmed.fastq"
    output:
        out = directory(final_path + "/kallisto/{sample}"),
        log = final_path + "/kallisto/{sample}_log.txt"
    params:
        is_paired = lambda wildcards: config.get("end", "pair") == "pair"
    benchmark:
        final_path + "/benchmarks/{sample}.kallisto.benchmark.txt"
    shell:
        """
        if [ "{params.is_paired}" == "True" ]; then
            kallisto quant --bias --single-overhang --fusion -i {input.index} -o {output.out} -t {config[NCORE]} {input.read1} {input.read2} &>{output.log}
        else
            kallisto quant --single --bias --single-overhang --fusion -i {input.index} -o {output.out} -t {config[NCORE]} -l 130 -s 30 {input.single} &>{output.log}
        fi
        """

rule multiqc_Report:
    input:
        count_summary = expand(final_path + "/kallisto/{sample}_log.txt", sample = samples)
    output:
        report = final_path + "/report_align_count.html"
    shell:
        "multiqc -f {input.count_summary} --filename {output.report}"