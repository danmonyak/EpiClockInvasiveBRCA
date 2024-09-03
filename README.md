# Measuring the Age of Individual Breast Cancers Using an Entropy-Based Molecular Clock

## 1. Introduction

The code in this repository can be used to generate all figures and results found in Monyak et al., 2024. Preprocessing of some of the data sources is performed by Python scripts, while all generation of figures and results must be done in Jupyter notebooks and R Markdown files.

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

## 3. Pipeline

### 1.

### 2.

### 3.

### 4.

### 5.

### 6.

### 7.

### 8.
