# High- Throughput data analysis pipeline

## How does this work?
Snakeamke wrapper scripts (locates in the `workflow` folder) enable for automatic RNA-seq data analysis in terms of quality control, assembly, quantification, gene ontology, differential gene expression and alternative splicing and it's effects on protein level. Additional 'RMarkdown' script enables final visualization for AS.
Additional Rmarkdown script allowas for Illumina microarrays analysis.

## Requirements
 - `conda` for building the environment
 - Python
 - interproscan

## 1. Installation
### 1. Clone the repository
 - `git clone www.github.com/aagatam/Pipeline.git`
### 2. Create a conda environment
 - `conda env create -f environment.yml`
