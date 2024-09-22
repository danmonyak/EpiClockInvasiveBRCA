# Use genefu to get PAM50 calls from RNA data

# Load libraries

library(tidyverse)
library(DESeq2)
library(genefu)

consts <- read_json(file.path(repo_dir, 'src', 'consts.json'))

# Load pam50 params -------------------------------------------------------

#- See https://rdrr.io/github/bhklab/genefu/man/pam50.html
#- List of parameters defining the PAM50 classifier for identification of breast cancer molecular subtypes (**Parker et al 2009**).
#- Three versions of the model are provided, each of ones differs by the gene expressions standardization method since it has an important impact on the subtype classification. I selected "pam50.robust".
#- `pam50.robust`: Use of the official centroids with robust scaling of the gene expressions (see rescale()) The model 'pam50.robust“ has been shown to reach the best concordance with the traditional clinical parameters (ER IHC, HER2 IHC/FISH and histological grade). However the use of this model is recommended only when the dataset is representative of a global population of breast cancer patients (no sampling bias, the 5 subtypes should be present) 
#- **Parker et al (2009)** "Supervised Risk Predictor of Breast Cancer Based on Intrinsic Subtypes", Journal of Clinical Oncology, 27(8):1160–1167
#- 5 Classes: "Basal"  "Her2"   "LumA"   "LumB"   "Normal"
#- 50 of 50 probes have a unique EntrezGeneID
data(pam50.robust)
pam50.genes <- pam50.robust$centroids.map$probe

# Import TPM data for TCGA breast samples
TCGA.tpm <- read.table(file.path(consts['TCGA_datadir'], 'TCGA.rnaseq_tpm.tsv'))

# > pam50.genes[!(pam50.genes %in% TCGA.tpm$gene_name)]
# [1] "CDCA1" "KNTC2" "ORC6L"

# CDCA1: ENSG00000143228.13
TCGA.tpm$gene_name[TCGA.tpm$gene_id == 'ENSG00000143228.13'] = 'CDCA1'
# KNTC2: ENSG00000080986.13
TCGA.tpm$gene_name[TCGA.tpm$gene_id == 'ENSG00000080986.13'] = 'KNTC2'
# ORC6L: ENSG00000091651.9
TCGA.tpm$gene_name[TCGA.tpm$gene_id == 'ENSG00000091651.9'] = 'ORC6L'

# > TCGA.tpm$gene_name[TCGA.tpm$gene_id == 'ENSG00000143228.13']
# [1] "NUF2"
# > TCGA.tpm$gene_name[TCGA.tpm$gene_id == 'ENSG00000080986.13']
# [1] "NDC80"
# > TCGA.tpm$gene_name[TCGA.tpm$gene_id == 'ENSG00000091651.9']
# [1] "ORC6"

all(table(TCGA.tpm$gene_name[TCGA.tpm$gene_name %in% pam50.genes]) == 1)
TCGA.tpm <- TCGA.tpm[!duplicated(TCGA.tpm$gene_name), ]
rownames(TCGA.tpm) <- TCGA.tpm$gene_name
TCGA.tpm <- TCGA.tpm[, !(names(TCGA.tpm) %in% c('gene_id', 'gene_name'))]

# # Load TCGA gene expression data
# gene.expr <- read.table(file.path(consts['TCGA_datadir'], 'gene_expr.txt'))
# 
# split.rownames = strsplit(rownames(gene.expr), '[|]')
# gene.name = c()
# ncbi.gene.id = c()
# for (sp.row in split.rownames) {
#   gene.name <- append(gene.name, sp.row[1])
#   ncbi.gene.id <- append(ncbi.gene.id, sp.row[2])
# }
# 
# gene.name[ncbi.gene.id == '83540'] = 'CDCA1'
# gene.name[ncbi.gene.id == '10403'] = 'KNTC2'
# gene.expr$gene.name <- gene.name
# 
# gene.expr <- gene.expr[isUnique(gene.expr$gene.name), ]
# rownames(gene.expr) <- gene.expr$gene.name
# gene.expr <- gene.expr[, names(gene.expr) != 'gene.name']

# Prep the count matrix ---------------------------------------------------

# expr.data <- gene.expr
# expr.data <- TCGA.tpm
expr.data <- log2(1 + TCGA.tpm)

# Check that your dataset rownames must match the probes
sum(pam50.genes %in% rownames(expr.data)) # 50 genes matches

# Remove empty features
ridx <- rowSums(expr.data) > 0
expr.data.c <- expr.data[ridx,]
expr <- expr.data.c
dim(expr)

# Subset to pam50 features
expr.sub <- expr[row.names(expr) %in% pam50.genes, ]

# Reformat expression matrix
ddata <- t(expr.sub)
dim(ddata)

# Run subtyping ---------------------------------------------------

set.seed(123456)
preds.pam50.robust <- molecular.subtyping(
  sbt.model = "pam50",
  data = ddata,
  annot = NULL,
  do.mapping = FALSE,
  verbose = TRUE
)
preds.pam50.robust

# Reformat results
probs <- preds.pam50.robust$subtype.proba
colnames(probs) <- paste0("subtype.pam50_proba_", colnames(probs))
df.pam50.fb <- data.frame(libName = names(preds.pam50.robust$subtype),
                          subtype.pam50 = preds.pam50.robust$subtype) %>%
  bind_cols(probs)

# write.table(df.pam50.fb, 'TCGA.genefu.pam50.subtypes.txt', quote=F, sep='\t', row.names=F)
# write.table(df.pam50.fb, 'TCGA.genefu.pam50.subtypes_from.expr.txt', quote=F, sep='\t', row.names=F)
write.table(df.pam50.fb, 'TCGA.genefu.pam50.subtypes_log2.norm.txt', quote=F, sep='\t', row.names=F)
