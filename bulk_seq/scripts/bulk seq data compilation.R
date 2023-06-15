# install packages, skip if all already installed
install.packages("tidyverse")
install.packages("here")
install.packages("readr")
install.packages("readxl")

bioPkgs <- c("TCGAbiolinks", "SummarizedExperiment", "spqn")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(bioPkgs)

library(tidyverse)
library(here)
library(readr)
library(readxl)
library(TCGAbiolinks)
library(SummarizedExperiment)

# seq data tsv download/unzipping/pre-processing
setwd(here::here("tnbc_bulk_rna_seq_data"))
clinical <- read_excel("organized_tnbc_clinical_data.xlsx")

set.seed(1234)
patients <- sample(clinical$"BCR patient barcode", 80)
# might have messed up here, need to make sure that the patients are the same for each session
# if not set seed to null, restart, change the query

query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts",
                  barcode = patients)
GDCdownload(query = query)
testdata <- GDCprepare(query = query) # all seq/clinical data as a SummarizedExperiment (SE)

ref <- read_tsv("gencode.gene.info.v22.tsv")  # gene lengths for RPKM calculation
ref <- ref[ref$gene_type %in% c("protein_coding", "miRNA"), ] %>% 
  distinct(gene_name, .keep_all = TRUE) %>%
  select(gene_name, gene_type, exon_length)

assays(testdata) <- assays(testdata)[-(2:6)]
counts <- as.data.frame(assay(testdata)) 
counts <- tibble::rownames_to_column(counts, "gene_name")
counts$gene_name <- rowData(testdata)$gene_name
counts <- distinct(counts, gene_name, .keep_all = TRUE) 
counts <- subset(counts, counts$gene_name %in% ref$gene_name)

ref <- ref[(ref$gene_name %in% counts$gene_name), ]
ref <- ref[order(match(ref$gene_name, counts$gene_name)), ]

patInfo <- colData(testdata)
colnames(counts)[-1] <- patInfo$patient
counts0 <- as.matrix(counts[ , colnames(counts)[-1]])
rownames(patInfo) <- colnames(counts0)
rownames(counts0) <- counts$gene_name

# create a new SE with only the valid genes
fullTnbc <- SummarizedExperiment(assays = SimpleList(exprs = counts0),
                                 rowData = ref,
                                 colData = patInfo)

save(fullTnbc, file = "fullTnbc.rda")
