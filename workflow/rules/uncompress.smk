########################################
# Uncompress input reads
########################################

COMPRESSION = config["input"]["compression"]
LAYOUT = config["sequencing"]["layout"]

def get_R1(wildcards):
    return samples_df.loc[
        samples_df["sample"] == wildcards.sample, "R1"
    ].values[0]

def get_R2(wildcards):
    return samples_df.loc[
        samples_df["sample"] == wildcards.sample, "R2"
    ].values[0]


########################################
# Paired-end
########################################

if LAYOUT == "paired":

    rule uncompress:
        input:
            r1 = get_R1,
            r2 = get_R2
        output:
            r1 = f"{OUTDIR}/{PROJECT}/reads/{{sample}}_R1.fastq",
            r2 = f"{OUTDIR}/{PROJECT}/reads/{{sample}}_R2.fastq"
        log:
            f"{OUTDIR}/{PROJECT}/logs/uncompress/{{sample}}.log"
        threads: 1
        shell:
            """
            mkdir -p {OUTDIR}/{PROJECT}/reads
            mkdir -p {OUTDIR}/{PROJECT}/logs/uncompress

            if [ "{COMPRESSION}" = "gz" ]; then
                gunzip -c {input.r1} > {output.r1}
                gunzip -c {input.r2} > {output.r2}
            elif [ "{COMPRESSION}" = "dsrc" ]; then
                dsrc d -i {input.r1} -o {output.r1}
                dsrc d -i {input.r2} -o {output.r2}
            else
                cp {input.r1} {output.r1}
                cp {input.r2} {output.r2}
            fi
            """

########################################
# Single-end
########################################

else:

    rule uncompress:
        input:
            r1 = get_R1
        output:
            r1 = f"{OUTDIR}/{PROJECT}/reads/{{sample}}_R1.fastq"
        log:
            f"{OUTDIR}/{PROJECT}/logs/uncompress/{{sample}}.log"
        threads: 1
        shell:
            """
            mkdir -p {OUTDIR}/{PROJECT}/reads
            mkdir -p {OUTDIR}/{PROJECT}/logs/uncompress

            if [ "{COMPRESSION}" = "gz" ]; then
                gunzip -c {input.r1} > {output.r1}
            elif [ "{COMPRESSION}" = "dsrc" ]; then
                dsrc d -i {input.r1} -o {output.r1}
            else
                cp {input.r1} {output.r1}
            fi
            """