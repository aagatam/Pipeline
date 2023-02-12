# High- throughput data analysis pipeline
# Contents
- [Requirements](#requirements)
- [Installation](#installation)
- [Configuration](#config)
- [Running the pipeline](#run)
  - [Trimming](#trim)
  - [QC raw reads](#qc)
  - [Alignment to genome](#genome)
  - [Alignment to transcriptome](#trans)
  - [Alternative splicing discovery](#asdis)
  - [Alternative splicing analysis](#asan)
  - [Associated protein domains analysis](#prot)
  - [Final plots](#plots)
  - [Differential expression analysis](#dif)
  - [Microarray differential expression](#micro)

## How does this work?
Snakeamke wrapper scripts (located in the `workflow` folder) enable for automatic RNA-seq data analysis in terms of quality control, assembly, quantification, gene ontology, differential gene expression and alternative splicing and it's effects on protein level. Additional 'RMarkdown' script enables final visualization for AS.
Additional Rmarkdown script allowas for Illumina microarrays analysis. Pipeline accepts either .fastq, fastq.gz or .fastq.dsrc files.

## 1. Requirements <a name="requirements"></a>
 - `conda` for building the environment
 - [`interproscan`](https://interproscan-docs.readthedocs.io/en/latest/UserDocs.html) and `pcregrep` for finding protein domains
 - [`dsrc`](https://github.com/refresh-bio/DSRC) if your files are dsrc compressed  

## 2. Installation <a name="installation"></a>
### 1. Clone the repository
 - `git clone www.github.com/aagatam/Pipeline.git`
### 2. Create a conda environment
 - `conda env create -f environment.yml`
### 3. Activate environment
 - `conda activate Pipeline`

## 3. Config files configuration <a name="config"></a>
 You need to adjust two files `configs/config.yaml` and `configs/Description.csv` to match your data for RNA-seq data analysis. For microarray analysis adjust  `configs/config_Illumina.yaml` and `configs/Description_Illumina.csv` All the fields in yaml configs are explained within the original file. Currently only analysis for samples with equal repetitions is available.

## 4. Running the pipeline <a name="run"></a>
## Trimming, quality control, assembly and quantification.

### 1. Trimming  <a name="trim"></a>
  - `snakemake --cores all --use-conda --conda-frontend conda -p -j 1 -s workflow/trim.snakefile`
#### **Outputs**
  - In `FINALOUTPUT`/`PROJECT`/trim:
    - trimmed fq.gz files,
    - quality report in `fastqc_after_trimming` folder

### 2. Quality contol on raw reads <a name="qc"></a>
 - `snakemake --cores all -p -s workflow/quality_control.snakefile`
#### **Outputs**
  - In `FINALOUTPUT`/`PROJECT`/fastqc:
    - summary report_quality_control.html,
    - report for each sample
### 3. Alignment to genome  <a name="genome"></a>
 - `snakemake --cores all -p -s workflow/align_HiSat2.snakefile`
#### **Outputs**
  - In `FINALOUTPUT`/`PROJECT`/genome:
    - multiqc summary report  report_align_count.html,
    - results prepared for analysis with R in `Hisat_results` folder,
    - Strintie results in `countFile` folder,
    - index in `indexes` folder,
    - sorted BAM files in `bamFIleSort`,
    - qualimap alignment QC in `alignmentQC` folder

### 4. Alignment to transcriptome <a name="trans"></a>
 - `snakemake --cores all -p -s workflow/align_kallisto.snakefile`
#### **Outputs**
  - In `FINALOUTPUT`/`PROJECT`/trans:
    - multiqc summary report report_align_count.html,
    - benchmarks folder,
    - kallisto index in `indexes` folder,
    - log files in `kallisto` folder,
    - folders with Kallisto results for all samples in `kallisto` folder

#### Alternative splicing and effects analysis
#### To perform ASEs analysis alignment to genome must be done!

### 5. Alternative splicing discovery with Spladder  <a name="asdis"></a>
This part will index all BAM files and run Spladder.
 - `snakemake --cores all -p -s workflow/spladder_run.snakefile`
#### **Outputs**
 - .bai index files for all BAM files,
 - Spladder output files in `FINALOUTPUT`/`PROJECT`/genome/spladder

### 6. Alternative events analysis with Bisbee and analysis with R <a name="asan"></a>
This part will run Bisbbe and analyse its output files with Rmarkdown script producing report and also bunch of csv files with results and files needed for further analysis with InterProscan and visualization. Depending on dataset size, this step might take a few hours, especially during first run when necessary libraries are downloaded. \
Install desired species release, for example:
 - `pyensembl install --release 104 --species mouse`\
 THEN
 - `snakemake --cores all -p -s workflow/bisbee_run.snakefile`
#### **Outputs**
   - In `FINALOUTPUT`/`PROJECT`/genome/bisbee:
    - csv files with bisbee results,
    - fasta files with transcripts including novel events.
  - In `FINALOUTPUT`/`PROJECT`/genome/spladder:
    - To_plots.csv <- file with all common events from new+new and new+old group
    - To_plots.RData <- for further visualizations with Plots.Rmd
  - In `FINALOUTPUT`/`PROJECT`/genome/spladder/Script_output:
    - .txt files with common genes for three groups and all events,
    - .csv with GO terms detected for each group,
    - .txt files with common GO terms between events within three groups,
    - .pdf files showing how GO terms change for first 10 terms from old+old group when adding events first from new+old group, then new+new and also top 10 terms from each group and how significant are they in others.
  - In `FINALOUTPUT`/`PROJECT`/genome/bisbee/Filetered:
    - bisbee results filtered with respect to valid events.
  - In `FINALOUTPUT`/`PROJECT`/genome/bisbee:
    - files  _to_grep.txt used for filtering fasta files for further InterProScan analysis.
  - Spladder.pdf - report in `scripts` folder

### 7. Protein domains affected by ASEs analysis with InterProScan <a name="prot"></a>
 - `snakemake --cores all -p -s workflow/interproscan_run.snakefile`
#### **Outputs**
   - In `FINALOUTPUT`/`PROJECT`/genome/bisbee:
    - .grepped.filtered.fasta - new fasta files with only intresting events, prepared for InterProScan analysis,
  - In `FINALOUTPUT`/`PROJECT`/genome/InterProScan:
   - .tsv files with InterProScan results

### 8. Final visualization with RMarkdown: <a name="plots"></a>
 - `R -e "rmarkdown::render('scripts/Plots.Rmd',params=list(event_type='event_type', event='event_no'),output_file='Out_name.pdf')"` \
 Here `event_no` is the event you want to visualize (for example mutex_exons_168) and `event_type` is one of: alt_3_prime, alt_5_prime, exon_skip, mult_exon_skip, mutex_exons (in this case mutex_exons). List of most interesting events is provided in `To_plots.csv` file from 6th step.
#### **Outputs**
   - Plots.pdf file with two plots for a given event. One for the whole transcript, and a close-up on the second one.

### 9. Differential gene expression and Gene Ontology analysis for RNA-seq <a name="dif"></a>
 - `R -e "rmarkdown::render('scripts/Expression_HiSat.Rmd')"`

 OR

 - `R -e "rmarkdown::render('scripts/Expression_Kallisto.Rmd')"`
#### **Outputs for both scripts**
   - pdf report in `scripts` folder,
   - results for limma, edgeR and DeSeq2 DEG (all and below given p-value) in `FINALOUTPUT`/`PROJECT`/genome/Hisat_results or `FINALOUTPUT`/`PROJECT`/trans/kallisto,
   - results for GO terms analysis in folders like above.

### 10. Differential gene expression and Gene Ontology for Illumina HumanHT-12 v4 microarrays <a name="micro"></a>
 - Change `config_Illumina.yaml` and `Description_Illumina.csv`,
 - `R -e "rmarkdown::render('scripts/Expression_Illumina_ microarrays.Rmd')"`
#### **Outputs**
   - Expression_Illumina_ microarrays.pdf report in `scripts` folder
