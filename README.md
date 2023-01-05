# High- throughput data analysis pipeline

## How does this work?
Snakeamke wrapper scripts (located in the `workflow` folder) enable for automatic RNA-seq data analysis in terms of quality control, assembly, quantification, gene ontology, differential gene expression and alternative splicing and it's effects on protein level. Additional 'RMarkdown' script enables final visualization for AS.
Additional Rmarkdown script allowas for Illumina microarrays analysis.

## Requirements
 - `conda` for building the environment
 - [`interproscan`](https://interproscan-docs.readthedocs.io/en/latest/UserDocs.html) for finding protein domains

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
 - `snakemake --cores all -p -s workflow/quality_control.snakefile`
### 2. Alignment to genome
 - `snakemake --cores all -p -s workflow/align_HiSat2.snakefile`
### 3. Alignment to transcriptome
 - `snakemake --cores all -p -s workflow/align_kallisto.snakefile`

## Alternative splicing and effects analysis
### To perform ASEs analysis alignment to genome must be done!

### 4. Alternative splicing discovery with Spladder and analysis with R
This part will index all BAM files, run Spladder and analyse its output files with Rmarkdown script producing report and also bunch of csv files with results and files needed for further analysis with InterProscan and visualization.
 - `snakemake --cores all -p -s workflow/spladder_run.snakefile`
 - Outputs:
  - In `FINALOUTPUT`/`PROJECT`/genome/spladder:
    - To_plots.csv <- file with all common events from new+new and new+old group
    - To_plots.Rmd <- for further visualizations with Plots.Rmd
  - In `FINALOUTPUT`/`PROJECT`/genome/spladder/Script_output:
    - .txt files with common genes for three groups and all events,
    - .csv with GO terms detected for each group,
    - .txt files with common GO terms between events within three groups,
    - .pdf files showing how GO terms change for first 10 terms from old+old group when adding events first from new+old group, then new+new and also top 10 terms from each group and how significant are they in others.
  - In `FINALOUTPUT`/`PROJECT`/genome/bisbee/Filetered:
    - bisbee results filtered with respect to valid events.
  - In `FINALOUTPUT`/`PROJECT`/genome/bisbee:
    - files  _to_grep.txt used for filtering fasta files for further InterProScan analysis.

### 5. Alternative events analysis with Bisbee
Install desired species release, for example:
 - `pyensembl install --release 104 --species musmusculus`\
 THEN
 - `snakemake --cores all -p -s workflow/bisbee_run.snakefile`
### 6. Protein domains affected by ASEs with InterProScan
 - `snakemake --cores all -p -s workflow/interproscan_run.snakefile`
### 7. Final visualization with RMarkdown:
 - `R -e "rmarkdown::render('scripts/Plots.Rmd',params=list(event_type='event_type', event='event_no'),output_file='Out_name.pdf')"` \
 Here `event_no` is the event you want to visualize (for example mutex_exons_168) and `event_type` is one of: alt_3_prime, alt_5_prime, exon_skip, mult_exon_skip, mutex_exons (in this case mutex_exons).

## Differential gene expression and Gene Ontology analysis for RNA-seq
 - `R -e "rmarkdown::render('scripts/Expression_HiSat.Rmd')"`
 - Output: `scripts/Expression_HiSat.pdf`

 OR

 - `R -e "rmarkdown::render('scripts/Expression_Kallisto.Rmd')"`
 - Output: `scripts/Expression_Kallisto.pdf`

## Differential gene expression and Gene Ontology for Illumina microarrays
 - `R -e "rmarkdown::render('scripts/Expression_Illumina_ microarrays.Rmd')"`
