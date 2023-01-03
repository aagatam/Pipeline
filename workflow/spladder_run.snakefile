import pandas as pd
configfile: "configs/config.yaml"

samples = pd.read_csv(config["METAFILE"], sep = ',', header = 0)['Group']
events = ['alt_3prime','alt_5prime','exon_skip','intron_retention','mult_exon_skip','mutex_exons']
final_path = config["FINALOUTPUT"] + "/" + config["PROJECT"] + "/genome"
splad_out = config["FINALOUTPUT"] + "/" + config["PROJECT"] + "/genome/spladder"
bisbee_out = config["FINALOUTPUT"] + "/" + config["PROJECT"] + "/genome/bisbee"
bb_path = config["BB_PATH"]

rule end:
    input:
        expand(splad_out+"/merge_graphs_{event}_C3.confirmed.txt.gz", event=events),
        expand(final_path + "/bamFileSort/{sample}.sort.bam.bai", sample=samples),
        'scripts/Spladder.pdf'

rule indexBams:
    input:
        BAM = final_path + "/bamFileSort/{sample}.sort.bam"
    output:
        BAI = final_path + "/bamFileSort/{sample}.sort.bam.bai"
    shell:
        "samtools index {input.BAM} > {output.BAI}"

rule spladderBuild:
    input:
        BAM = expand(final_path+"/bamFileSort/{sample}.sort.bam", sample=samples)
    output:
        expand(splad_out+"/merge_graphs_{event}_C3.confirmed.txt.gz", event=events)
    params:
        GTF = config["ANNOTATION"],
        files=lambda wildcards, input: ','.join(input),
        BAI = expand(final_path + "/bamFileSort/{sample}.sort.bam.bai", sample=samples),
        out = splad_out
    shell:
        "spladder build --bams {params.files} --annotation {params.GTF} --outdir {params.out}" # no paralel for dumb mac
        #"spladder build --parallel {config[NCORE]} --bams {params.files} --annotation {params.GTF} --outdir {params.out}"

rule markdown:
    input:
        "scripts/Spladder.Rmd"
    output:
        'scripts/Spladder.pdf'
    shell:
        "Rscript -e \"rmarkdown::render('scripts/Spladder.Rmd')\""
