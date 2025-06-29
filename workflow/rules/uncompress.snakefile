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