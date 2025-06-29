rule uncompress:
    input:
        lambda wildcards: get_input_path(wildcards)
    output:
        uncompress1 = temp(f"{final_path}/uncompressed/{{sample}}_R1.out.fastq") if end == "pair" else temp(f"{final_path}/uncompressed/{{sample}}.out.fastq"),
        uncompress2 = temp(f"{final_path}/uncompressed/{{sample}}_R2.out.fastq") if end == "pair" and not single_end else None
    run:
        decompress_command(input, output)
