# High- throughput data analysis pipeline

## How does this work?
Snakeamke wrapper scripts (locates in the `workflow` folder) enable for automatic RNA-seq data analysis in terms of quality control, assembly, quantification, gene ontology, differential gene expression and alternative splicing and it's effects on protein level. Additional 'RMarkdown' script enables final visualization for AS.
Additional Rmarkdown script allowas for Illumina microarrays analysis.

## Requirements
 - `conda` for building the environment

## 1. Installation
### 1. Clone the repository
 - `git clone www.github.com/aagatam/Pipeline.git`
### 2. Create a conda environment
 - `conda env create -f environment.yml`
### 3. Activate environment
 - `conda activate Pipeline`

## 2. Config files configuration
 You need to adjust two files `configs/config.yaml` and `configs/Description.csv` to match your data. All the fields in first config are explained within the original file. Currently only analysis for samples with equal repetitions is available.

## 3. Running the pipeline.
## Quality control, assebmle and quantification.
### 1. Quality contol on raw reads
 - `snakemake -p -s workflow/quality_control.snakefile`
### 2. Alignment to genome
 - `snakemake -p -s workflow/align_HiSat2.snakefile`
### 3. Alignment to transcriptome
 - `snakemake -p -s workflow/align_kallisto.snakefile`

## Alternative splicing and effects analysis
### To perform ASEs analysis alignment to genome must be done!

### 4. Alternative splicing discovery with Spladder and analysis with R
 - `snakemake -p -s workflow/spladder_run.snakefile`
### 5. Alternative events analysis with Bisbee
Install desired species release, for example:\
 - `pyensembl install --release 104 --species musmusculus`
 - `snakemake -p -s workflow/bisbee_run.snakefile`
### 6. Protein domains affected by ASEs with InterProScan
 - `snakemake -p -s workflow/interproscan_run.snakefile`
### 7. Final visualization with RMarkdown:
 - `R -e "rmarkdown::render('scripts/Plots.Rmd',params=list(event_type='event_type', event='event_no'),output_file='Out_name.pdf')"` \
 Here `event_no` is the event you want to visualize (for example mutex_exons_168) and `event_type` is one of: alt_3_prime, alt_5_prime, exon_skip, mult_exon_skip, mutex_exons (in this case mutex_exons).

## Differential gene expression and Gene Ontology analysis for RNA-seq

## Differential gene expression and Gene Ontology for Illumina microarrays
