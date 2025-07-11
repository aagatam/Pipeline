import pandas as pd
configfile: "configs/config.yaml"
events = ['alt_3prime','alt_5prime','exon_skip','mult_exon_skip','mutex_exons']#'intron_retention',
splad_out = config["FINALOUTPUT"] + "/" + config["PROJECT"] + "/genome/spladder"
bisbee_out = config["FINALOUTPUT"] + "/" + config["PROJECT"] + "/genome/bisbee"
bb_path = config["BB_PATH"]

rule end:
    input:
        expand(bisbee_out + "/{event}.effects.csv",event=events),
        expand(bisbee_out + "/{event}.peptides.csv",event=events),
        'scripts/Spladder.pdf'


rule bisbeePrep:
    input:
        HDF5 = splad_out + "/merge_graphs_{event}_C3.counts.hdf5"
    output:
        CSV = bisbee_out + "/{event}.bisbeeCounts.csv"
    params:
        BB = bb_path+ "/utils",
        EV = "{event}"
    shell:
        "python {params.BB}/prep.py {input.HDF5} {params.EV} {output.CSV} 2"

rule bisbeeProt:
    input:
        CSV = bisbee_out + "/{event}.bisbeeCounts.csv"
    output:
        CSV_EFF = bisbee_out + "/{event}.effects.csv",
        CSV_PROT = bisbee_out + "/{event}.peptides.csv"
    params:
        BB = bb_path + "/prot",
        genome = config["GENOME"],
        name= bisbee_out + "/{event}",
        release = config['RELEASE'],
        EV = "{event}"
    shell:
        "python {params.BB}/build.py {input.CSV} {params.EV} 9 {params.name} {params.release} {params.genome} "

rule markdown:
    input:
        CSV_PROT = bisbee_out + "/{event}.peptides.csv",
        "scripts/Spladder.Rmd"
    output:
        'scripts/Spladder.pdf'
    run:
        shell("touch {input.CSV_PROT}")
        shell("Rscript -e \"rmarkdown::render('scripts/analysis/Spladder.Rmd')\"")
