def get_input_files(wildcards, config, input_path, intermediate_path):
    """Returns a list of input file paths for raw samples."""
    sample = wildcards.sample
    if config.get("end", "pair") == "pair":
        if config.get("compression") == "dsrc":
            return [
                f"{input_path}/{sample}_R1.fastq.dsrc",
                f"{input_path}/{sample}_R2.fastq.dsrc"
            ]
        elif config.get("compression") == "gz":
            return [
                f"{input_path}/{sample}_R1.fastq.gz",
                f"{input_path}/{sample}_R2.fastq.gz"
            ]
        else:
            return [
                f"{input_path}/{sample}_R1.fastq",
                f"{input_path}/{sample}_R2.fastq"
            ]
    else:
        if config.get("compression") == "dsrc":
            return [f"{input_path}/{sample}.fastq.dsrc"]
        elif config.get("compression") == "gz":
            return [f"{input_path}/{sample}.fastq.gz"]
        else:
            return [f"{input_path}/{sample}.fastq"]

def get_trimmed_input_files(wildcards, config, intermediate_path):
    """Returns a list of trimmed compressed file paths."""
    sample = wildcards.sample
    if config.get("end", "pair") == "pair":
        return [
            f"{intermediate_path}/{sample}_val_1.fq.gz",
            f"{intermediate_path}/{sample}_val_2.fq.gz"
        ]
    else:
        return [f"{intermediate_path}/{sample}_trimmed.fq.gz"]

def decompress_command(input_dict, output_dict, config):
    """Executes decompression commands."""
    if config.get("end", "pair") == "pair":
        if config.get("compression") == "dsrc":
            shell(f"dsrc d -t{config['NCORE']} -s {input_dict['forward']} > {output_dict['uncompress1']}")
            shell(f"dsrc d -t{config['NCORE']} -s {input_dict['reverse']} > {output_dict['uncompress2']}")
        elif config.get("compression") == "gz":
            shell(f"pigz -d -k -c -p{config['NCORE']} {input_dict['forward']} > {output_dict['uncompress1']}")
            shell(f"pigz -d -k -c -p{config['NCORE']} {input_dict['reverse']} > {output_dict['uncompress2']}")
        else:
            shell(f"ln -s {input_dict['forward']} {output_dict['uncompress1']}")
            shell(f"ln -s {input_dict['reverse']} {output_dict['uncompress2']}")
    else:
        if config.get("compression") == "dsrc":
            shell(f"dsrc d -t{config['NCORE']} -s {input_dict['read']} > {output_dict['uncompress']}")
        elif config.get("compression") == "gz":
            shell(f"pigz -d -k -c -p{config['NCORE']} {input_dict['read']} > {output_dict['uncompress']}")
        else:
            shell(f"ln -s {input_dict['read']} {output_dict['uncompress']}")

def decompress_trimmed_command(input_dict, output_dict, config):
    """Executes decompression commands for trimmed files (always .gz)."""
    if config.get("end", "pair") == "pair":
        shell(f"pigz -d -k -c -p{config['NCORE']} {input_dict['forward']} > {output_dict['uncompress1']}")
        shell(f"pigz -d -k -c -p{config['NCORE']} {input_dict['reverse']} > {output_dict['uncompress2']}")
    else:
        shell(f"pigz -d -k -c -p{config['NCORE']} {input_dict['read']} > {output_dict['uncompress']}")