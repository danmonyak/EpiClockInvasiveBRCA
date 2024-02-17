# install-packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# if (!requireNamespace("devtools"))
#   install.packages('devtools')
# 
# if (!require("dbplyr", quietly = TRUE))
#   devtools::install_version("dbplyr", version = "2.3.4")

if (!require("TCGAbiolinks", quietly = TRUE))
  BiocManager::install("TCGAbiolinks", force = TRUE)

# settings
# This package streamlines the data pull from GDC
library(TCGAbiolinks)


ductal_only <- FALSE
female_only <- TRUE
project <- "TCGA-BRCA"

# directory where all files will be saved
save_path <- "data_objects/"
dir.create(save_path, showWarnings = FALSE)

# Do not add file type ending -- just name
methyl_filename <- "cohort1.methyl"
# Do not add file type ending -- just name
gene_filename <- "cohort1.rnaseq"


# query_clinical_data
query_clinical <- GDCquery_clinic(project = project, type = "clinical")

# we are only interest in primary tumor_subset tumors
if (ductal_only) {
  ductal <- !is.na(query_clinical$primary_diagnosis) & 
    query_clinical$primary_diagnosis == "Infiltrating duct carcinoma, NOS"
} else {
  ductal <- rep(TRUE, length(query_clinical$primary_diagnosis))
}

if (female_only) {
  female <- !is.na(query_clinical$gender) & query_clinical$gender == "female"
} else {
  female <- rep(TRUE, length(query_clinical$gender))
}
tumor_subset <- ductal & female

# query_gene_data
query_gene <- GDCquery(project = project, 
                       data.category = 'Transcriptome Profiling', 
                       data.type = 'Gene Expression Quantification', 
                       workflow.type = 'STAR - Counts',
                       experimental.strategy = "RNA-Seq",
                       sample.type = c("Primary Tumor"))
gene_ids <- substr(getResults(query_gene, cols = "cases"), 1, 12)

#query_methyl_data
query_methyl <- GDCquery(project = project, 
                         data.category = 'DNA Methylation', 
                         platform = c("Illumina Human Methylation 450"),
                         data.type = "Methylation Beta Value",
                         sample.type = c("Primary Tumor"))
methyl_ids <- substr(getResults(query_methyl, cols = "cases"), 1, 12)


# must have both gene data and methylation data
common.patients <- intersect(
  substr(getResults(query_methyl, cols = "cases"), 1, 12),
  substr(getResults(query_gene, cols = "cases"), 1, 12)
)

# limit_to_tumor_subset
# must be tumor_subset
common.patients <- intersect(
  common.patients,
  query_clinical$submitter_id[tumor_subset]
)

# redefine queries with only the common samples
query_gene <- GDCquery(project = project, 
                       data.category = 'Transcriptome Profiling', 
                       data.type = 'Gene Expression Quantification', 
                       workflow.type = 'STAR - Counts',
                       experimental.strategy = "RNA-Seq",
                       sample.type = c("Primary Tumor"),
                       barcode = common.patients)

query_methyl <- GDCquery(project = project, 
                         data.category = 'DNA Methylation', 
                         platform = c("Illumina Human Methylation 450"),
                         sample.type = c("Primary Tumor"),
                         data.type = "Methylation Beta Value",
                         barcode = common.patients)

# save the clinical data
clinical_final <- query_clinical[query_clinical$bcr_patient_barcode %in% common.patients, ]
write.table(clinical_final, paste0(save_path, 'sample_annotations.txt'), quote=FALSE)

# download and save the gene expression data
GDCdownload(query_gene)
data.seq <- GDCprepare(query_gene, save = TRUE, 
                       save.filename = paste0(save_path, gene_filename, ".RData"))

# use the tpm data for methylation analysis
tpm <- data.seq@assays@data$tpm_unstrand
colnames(tpm) <- rownames(data.seq@colData)
rm_pary <- grepl("PAR_Y", data.seq@rowRanges$gene_id)
rownames(tpm) <- substr(data.seq@rowRanges$gene_id, 1, 15)
tpm <- tpm[!rm_pary, ]

saveRDS(tpm, file = paste0(save_path, gene_filename, "_tpm.rds"))

# download and save the DNA methylation data
GDCdownload(query_methyl)
data.methyl <- GDCprepare(query_methyl, save = TRUE, 
                          save.filename = paste0(save_path, methyl_filename, ".RData"))
methyl <- data.methyl@assays@data[[1L]]
saveRDS(methyl, file = paste0(save_path, methyl_filename, ".rds"))
