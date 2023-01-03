import pandas as pd
configfile: "configs/config.yaml"

samples = pd.read_csv(config["METAFILE"], sep = ';', header = 0)['Sample']
end = config["END"]
index_path = config["FINALOUTPUT"] + "/" + config["PROJECT"] + "/trans"
compressed = config["COMPRESSED"]
final_path = config["FINALOUTPUT"] + "/" + config["PROJECT"] + "/trans"
input_path = config["INPUTPATH"]

rule end:
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
if compressed == "yes":
    rule uncompress:
        input:
            read = input_path + "/{sample}.fastq.dsrc"
        output:
            uncompress = temp(final_path + "/uncompressed/{sample}.out.fastq"),
        shell:
            "dsrc d -t{config[NCORE]} -s {input.read} >>{output.uncompress} "
else:
    rule uncompress:
        input:
            read = input_path + "/{sample}.fastq"
        output:
            uncompress =  temp(final_path + "/uncompressed/{sample}.out.fastq")
        shell:
            "ln -s {input.read} {output.uncompress}"

rule alignment:
    input:
        index = index_path + "/indexes/kallisto_index.idx",
        uncompress = temp(final_path + "/uncompressed/{sample}.out.fastq")
    output:
        out = directory(final_path + "/kallisto/{sample}"),
        log = final_path + "/kallisto/{sample}_log.txt"
    benchmark:
        final_path + "/benchmarks/{sample}.kallisto.benchmark.txt"
    shell:
        "kallisto quant --single --bias --single-overhang --fusion -i {input.index} -o {output.out}  -t 32 -l 130 -s 30 {input.uncompress} &>{output.log}"

rule summaryReport:
    input:
        count_summary = expand(final_path + "/kallisto/{sample}_log.txt", sample = samples)
    output:
        report = final_path + "/report_align_count.html"
    shell:
        "multiqc -f {input.count_summary} --filename {output.report}"
