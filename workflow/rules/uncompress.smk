# First decompression: Raw files → uncompressed files (for trimming)
rule uncompress:
    input:
        # Pass wildcards AND config to the lambda
        lambda wildcards, config=config: get_input_files(wildcards, config, input_path, config["INTERMEDIATEPATH"])
    output:
        uncompressed1 = temp(f"{final_path}/uncompressed/{{sample}}_R1.out.fastq"),
        uncompressed2 = temp(f"{final_path}/uncompressed/{{sample}}_R2.out.fastq"),
        uncompressed_single = temp(f"{final_path}/uncompressed/{{sample}}.out.fastq")
    params:
        is_paired = lambda wildcards, config=config: config.get("end", "pair") == "pair"
    run:
        input_files = input[0]
        if params.is_paired(wildcards, config):
            decompress_command(
                {"forward": input_files[0], "reverse": input_files[1]},
                {"uncompress1": output.uncompressed1, "uncompress2": output.uncompressed2},
                config
            )
        else:
            decompress_command(
                {"read": input_files[0]},
                {"uncompress": output.uncompressed_single},
                config
            )

# Second decompression: Trimmed compressed files → uncompressed trimmed files (for QC and alignment)
rule uncompress_trimmed:
    input:
        # Use the static outputs from trim rule
        trimmed_1 = f"{intermediate_path}/{{sample}}_val_1.fq.gz",
        trimmed_2 = f"{intermediate_path}/{{sample}}_val_2.fq.gz",
        trimmed_single = f"{intermediate_path}/{{sample}}_trimmed.fq.gz"
    output:
        uncompressed1 = temp(f"{final_path}/trimmed_uncompressed/{{sample}}_val_1.fastq"),
        uncompressed2 = temp(f"{final_path}/trimmed_uncompressed/{{sample}}_val_2.fastq"),
        uncompressed_single = temp(f"{final_path}/trimmed_uncompressed/{{sample}}_trimmed.fastq")
    params:
        is_paired = lambda wildcards, config=config: config.get("end", "pair") == "pair"
    run:
        if params.is_paired(wildcards, config):
            decompress_trimmed_command(
                {"forward": input.trimmed_1, "reverse": input.trimmed_2},
                {"uncompress1": output.uncompressed1, "uncompress2": output.uncompressed2},
                config
            )
        else:
            decompress_trimmed_command(
                {"read": input.trimmed_single},
                {"uncompress": output.uncompressed_single},
                config
            )