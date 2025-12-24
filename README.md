# Measuring the Age of Individual Breast Cancers Using an Entropy-Based Molecular Clock
# EpiClockInvasiveBRCA

[![DOI](https://zenodo.org/badge/758688116.svg)](https://doi.org/10.5281/zenodo.16813782)

## 1. Introduction

The code in this repository can be used to generate all figures and results found in Monyak et al. (2025). Preprocessing of some of the data sources is performed by Python scripts, while all generation of figures and results must be done in Jupyter notebooks and R Markdown files.

Software requirements:
- Python 3
- R

## 2. Setup

Fork and clone this repository locally as normal.

### Python

Use a bash shell to run all scripts and Jupyter notebooks. To see what shell is running, use ```echo $SHELL```. Run the following line to append the path to the parent directory of the repository clone to the Python path:

```
repo_parent_dir=/PATH/TO/REPO/PARENT/DIR
echo "export PYTHONPATH=$PYTHONPATH:$repo_parent_dir" >> ~/.bash_profile
```

**Note**: the local clone of EpiClockInvasiveBRCA should be located directly in *repo_parent_dir*.

### R

In your R environment, preferably Rstudio, run the following line and copy the path outputted:

```
file.path(Sys.getenv("R_HOME"), 'etc', 'Rprofile.site')
```

Append the following line to the file at the path outputted (create the file if necessary):

```
repo_dir <- '/PATH/TO/REPO/PARENT/DIR/EpiClockInvasiveBRCA'
```

replacing the code above appropriately with the path to the local repository clone.

### Path variables

Open src/consts.json in a text editor and insert appropriate paths for the following attributes:
- **repo_dir** — Path to the repository (same as the R variable "repo_dir" in the previous step)
- **official_indir** — Path to a directory in an external file location (preferably Box) that can hold terabytes of data
- **TCGA_datadir** — Path to a directory that will hold the TCGA data (preferably a subdirectory of official_indir)
- **Lund_datadir** — Path to a directory that will hold the Lund cohort data (preferably a subdirectory of official_indir)

## 2. Supplementary Data Retrieval

In the external file location (preferably Box), create directories with the name of each cohort, to which one should download the relevant data from Gene Expression Omnibus (GEO):
- **Aurora**: GSE212370
  - [Methylation](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE212nnn/GSE212370/suppl/GSE212370%5F20220828%5FMethylationEPIC%5FGEO%5Fdata%5Fmatrix146.tsv.gz)
  - Clinical - Download supplementary data directory [43018_2022_491_MOESM2_ESM](https://static-content.springer.com/esm/art%3A10.1038%2Fs43018-022-00491-x/MediaObjects/43018_2022_491_MOESM2_ESM.rar) from https://doi.org/10.1038/s43018-022-00491-x
- **Desmedt**: GSE39451
  - [Series Matrix](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE39nnn/GSE39451/matrix/GSE39451-GPL13534_series_matrix.txt.gz)
- **Germany**: GSE69914
  - [Series Matrix](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE69nnn/GSE69914/matrix/GSE69914_series_matrix.txt.gz)
- **Lund**: GSE25307
  - [Methylation](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE25nnn/GSE25307/matrix/GSE25307_series_matrix.txt.gz)
  - [Clinical](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE25nnn/GSE25307/suppl/GSE25307%5Fsample%5Fannotations.txt.gz)
- **Luo**: GSE106360
  - [Series Matrix](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE106nnn/GSE106360/matrix/GSE106360_series_matrix.txt.gz)
- **Reyngold**: GSE58999
  - [Series Matrix](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE58nnn/GSE58999/matrix/GSE58999_series_matrix.txt.gz)

## 3. Pipeline

### 1. Simulation

To run all simulations, do:
```
sh runAllSimulations.sh
```

### 2. TCGA Retrieval

To retrieve the TCGA data and generate the HTML output, set the header parameters accordingly in Data_Prep.Rmd, and in bash, do:
```
Rscript -e "rmarkdown::render('Data_Prep.Rmd', output_format = 'html_document', output_file = paste0('Data_Prep ', Sys.time(), '.html'))"
```

It could take up to a few hours to run, though it will likely take less than 1 hour. This script should be run on a machine of at least 16 GB of memory.

### 3. Select fCpGs

Run all cells in the Jupyter notebook "Select_fCpGs-Revision.ipynb"

### 4. Process Supplementary Data

To process all supplementary data, do:
```
sh processAllData.sh
```

### 5. Subtyping

To calculate PAM50 subtype for the TCGA tumors, do:

```
Rscript "subtype.R"
```

### 6. Beta Mixture Model

To run the beta mixture model decomposition analysis on the TCGA and Lund data and generate the HTML output, set the header parameters accordingly in Fit_BetaMixture.Rmd, and in bash, do:
```
Rscript -e "rmarkdown::render('Fit_BetaMixture.Rmd', output_format = 'html_document', output_file = paste0('Fit_BetaMixture ', Sys.time(), '.html'))"
```

It should take no more than 10 minutes to render.

### 7. Analysis

**In this order:**

Run all cells in:
1. c_beta Analysis.ipynb
2. Estimate ages.ipynb
3. Multi-sample.ipynb


To perform beta value adjustment and subsequent analysis, and generate the HTML output, set the header parameters accordingly in beta_adjustment.Rmd, and in bash, do:

```
Rscript -e "rmarkdown::render('beta_adjustment.Rmd', output_format = 'html_document', output_file = paste0('beta_adjustment ', Sys.time(), '.html'))"
```

To generate the GSEA-related figures and generate the HTML output, set the header parameters accordingly in GSEA_Figure.Rmd, and in bash, do:

```
Rscript -e "rmarkdown::render('GSEA_Figure.Rmd', output_format = 'html_document', output_file = paste0('GSEA_Figure ', Sys.time(), '.html'))"
```

