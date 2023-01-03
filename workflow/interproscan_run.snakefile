import pandas as pd
configfile: "configs/config.yaml"

#samples = pd.read_csv(config["METAFILE"], sep = ';', header = 0)['Sample']
events = ['alt_5prime','exon_skip','mult_exon_skip','mutex_exons']
bisbee_out = config["FINALOUTPUT"] + "/" + config["PROJECT"] + "/genome/bisbee"
interpro_out = config["FINALOUTPUT"] + "/" + config["PROJECT"] + "/genome/InterProScan"
interpro_path = config["IP_PATH"]
bb_path = config["BB_PATH"]


rule end:
    input:
        expand(interpro_out + "/{event}.changeSeq.fasta.tsv",event=events),

rule removeDuplicates:
    input:
        FASTAIN = bisbee_out + "/{event}.changeSeq.fasta"
    output:
        FASTAOUT = bisbee_out + "/{event}.changeSeq.nodup.fasta"
    params:
        BB = bb_path + "/prot",
        name= bisbee_out + "/{event}"
    shell:
        "python {params.BB}/remove_duplicates.py {input.FASTAIN} {params.name}"

rule wrapFasta:
    input:
        FASTAIN = bisbee_out + "/{event}.changeSeq.nodup.fasta"
    output:
        FASTAOUT = bisbee_out + "/{event}.wrapped.fasta"
    shell:
        """ sed -e 's/\\(^>.*$\\)/#\\1#/' {input.FASTAIN} | tr -d "\\r" | tr -d "\\n" | sed -e 's/$/#/' | tr "#" "\\n" | sed -e '/^$/d' > {output.FASTAOUT} """ #wrap FASTA

rule grepFasta:
    input:
        FASTAIN = bisbee_out + "/{event}.wrapped.fasta"
    output:
        FASTAFILT = bisbee_out + "/{event}.grepped.filtered.fasta"
    params:
        GREP = bisbee_out + "/{event}_newAndOld_bisbee_filtered_to_grep.txt"
    run:
        shell("cat {input.FASTAIN} | pcregrep -M -A 1 --buffer-size=40000 -f {params.GREP} > {output.FASTAFILT}") #grep only those common
        shell("sed -i 's/--//' {output.FASTAFILT}") #sed --
        shell("sed -i 's/*//' {output.FASTAFILT}") #sed *

rule runInterPro:
    input:
        FASTAIN = bisbee_out + "/{event}.grepped.filtered.fasta"
    output:
        TSVOUT = interpro_out + "/{event}.changeSeq.fasta.tsv"
    params:
        IP = interpro_path + "/interproscan.sh",
        NAMEOUT = interpro_out + "/{event}.changeSeq.fasta"
    shell:
        "{params.IP} -i {input.FASTAIN} -b {params.NAMEOUT} -f tsv -goterms"
