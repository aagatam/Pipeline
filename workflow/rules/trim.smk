# workflow/rules/trim.snakefile

rule trim:
    input:
        # Explicit input files based on paired/single end
        read1 = f"{final_path}/uncompressed/{{sample}}_R1.out.fastq",
        read2 = f"{final_path}/uncompressed/{{sample}}_R2.out.fastq",
        single = f"{final_path}/uncompressed/{{sample}}.out.fastq"
    output:
        # Static output - all possible outputs, but only relevant ones will be created
        trimmed_1 = f"{intermediate_path}/{{sample}}_val_1.fq.gz",
        trimmed_2 = f"{intermediate_path}/{{sample}}_val_2.fq.gz", 
        trimmed_single = f"{intermediate_path}/{{sample}}_trimmed.fq.gz"
    params:
        outputpath = intermediate_path,
        is_paired = lambda wildcards: config.get("end", "pair") == "pair"
    conda: "../configs/trim_env.yaml"
    shell:
        """
        if [ "{params.is_paired}" == "True" ]; then
            trim_galore --fastqc -j {config[NCORE]} --gzip --paired \
                --basename {wildcards.sample} -o {params.outputpath} \
                {input.read1} {input.read2}
        else
            trim_galore --fastqc -j {config[NCORE]} --gzip \
                --basename {wildcards.sample} -o {params.outputpath} \
                {input.single}
        fi
        """