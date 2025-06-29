def get_input_path(wildcards):
    sample = wildcards.sample
    if end == "pair":
        if compression == "dsrc":
            return {
                "forward": f"{input_path}/{sample}_R1.fastq.dsrc",
                "reverse": f"{input_path}/{sample}_R2.fastq.dsrc"
            }
        elif compression == "gz":
            return {
                "forward": f"{input_path}/{sample}_R1.fastq.gz",
                "reverse": f"{input_path}/{sample}_R2.fastq.gz"
            }
        elif trimmed == "yes":
            return {
                "read_trim_forward": f"{intermediate_path}/{sample}_val_1.fq.gz",
                "read_trim_reverse": f"{intermediate_path}/{sample}_val_2.fq.gz"
            }
        else:
            return {
                "forward": f"{input_path}/{sample}_R1.fastq",
                "reverse": f"{input_path}/{sample}_R2.fastq"
            }
    else:
        if compression == "dsrc":
            return {"read": f"{input_path}/{sample}.fastq.dsrc"}
        elif compression == "gz":
            return {"read": f"{input_path}/{sample}.fastq.gz"}
        elif trimmed == "yes":
            return {"read_trim": f"{intermediate_path}/{sample}_trimmed.fq"}
        else:
            return {"read": f"{input_path}/{sample}.fastq"}


def decompress_command(input, output):
    if end == "pair":
        if compression == "dsrc":
            shell(f"dsrc d -t{{config[NCORE]}} -s {input.forward} >> {output.uncompress1}")
            shell(f"dsrc d -t{{config[NCORE]}} -s {input.reverse} >> {output.uncompress2}")
        elif compression == "gz":
            shell(f"pigz -d -k -c -p{{config[NCORE]}} {input.forward} > {output.uncompress1}")
            shell(f"pigz -d -k -c -p{{config[NCORE]}} {input.reverse} > {output.uncompress2}")
        elif trimmed == "yes":
            shell(f"pigz -d -k -c -p{{config[NCORE]}} {input.read_trim_forward} > {output.uncompress1}")
            shell(f"pigz -d -k -c -p{{config[NCORE]}} {input.read_trim_reverse} > {output.uncompress2}")
        else:
            shell(f"ln -s {input.forward} {output.uncompress1}")
            shell(f"ln -s {input.reverse} {output.uncompress2}")
    else:
        if compression == "dsrc":
            shell(f"dsrc d -t{{config[NCORE]}} -s {input.read} >> {output.uncompress}")
        elif compression == "gz":
            shell(f"pigz -d -k -c -p{{config[NCORE]}} {input.read} > {output.uncompress}")
        elif trimmed == "yes":
            shell(f"pigz -d -k -c -p{{config[NCORE]}} {input.read_trim} > {output.uncompress}")
        else:
            shell(f"ln -s {input.read} {output.uncompress}")
