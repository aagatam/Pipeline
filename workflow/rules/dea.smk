########################################
# Differential expression analysis
########################################

THREADS = config["alignment"]["threads"]

if config["alignment"]["method"] == "kallisto":

    rule differential_expression:
        input:
            expand(
                f"{OUTDIR}/{PROJECT}/trans/kallisto/{{sample}}/abundance.tsv",
                sample=SAMPLES
            )
        output:
            f"{OUTDIR}/{PROJECT}/dea/Expression_RNAseq.pdf"
        params:
            samples = config["input"]["samples"]
        conda:
            "../envs/dea_env.yaml"
        shell:
            """
            mkdir -p {OUTDIR}/{PROJECT}/dea

            Rscript -e "rmarkdown::render(
                'scripts/analysis/Expression_RNAseq.Rmd',
                knit_root_dir = getwd(),
                params = list(
                    quant_method = 'kallisto',
                    config_path = 'config',
                    samples_path = '{params.samples}',
                    utils_path = 'scripts/utils.R'
                ),
                output_file = 'Expression_RNAseq.pdf',
                output_dir = '{OUTDIR}/{PROJECT}/dea'
            )"
            """

elif config["alignment"]["method"] == "hisat2":

    rule differential_expression:
        input:
            f"{OUTDIR}/{PROJECT}/genome/gene_count_matrix.csv"
        output:
            f"{OUTDIR}/{PROJECT}/dea/Expression_RNAseq.pdf"
        params:
            samples = config["input"]["samples"]
        conda:
            "../envs/dea_env.yaml"
        shell:
            """
            mkdir -p {OUTDIR}/{PROJECT}/dea

            Rscript -e "rmarkdown::render(
                'scripts/analysis/Expression_RNAseq.Rmd',
                knit_root_dir = getwd(),
                params = list(
                    quant_method = 'hisat2',
                    config_path = 'config',
                    samples_path = '{params.samples}',
                    utils_path = 'scripts/utils.R'
                ),
                output_file = 'Expression_RNAseq.pdf',
                output_dir = '{OUTDIR}/{PROJECT}/dea'
            )"
            """