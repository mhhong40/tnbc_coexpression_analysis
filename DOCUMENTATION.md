# TNBC-RNA-seq-analysis
'23 spring -> present research with Dr. Gabor Balazsi @ Stony Brook University

[DOCUMENTATION]

## Handling bulk RNA-seq data from the GDC Repository:
### Preliminaries:
  - Retrieve TCGA-BRCA clinical data from [the NIH GDC Repository](https://portal.gdc.cancer.gov/repository?filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22content%22%3A%7B%22field%22%3A%22cases.project.project_id%22%2C%22value%22%3A%5B%22TCGA-BRCA%22%5D%7D%2C%22op%22%3A%22in%22%7D%2C%7B%22content%22%3A%7B%22field%22%3A%22files.cases.primary_site%22%2C%22value%22%3A%5B%22breast%22%5D%7D%2C%22op%22%3A%22in%22%7D%2C%7B%22content%22%3A%7B%22field%22%3A%22files.data_type%22%2C%22value%22%3A%5B%22Gene%20Expression%20Quantification%22%5D%7D%2C%22op%22%3A%22in%22%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.experimental_strategy%22%2C%22value%22%3A%5B%22RNA-Seq%22%5D%7D%7D%5D%7D&searchTableTab=files)
  - Transfer clinical .txt file data into workable Excel spreadsheet
  - Sort by IHC results to select for just TNBC (total: n = 115 patients)
  
  - Define genes of interest:
      - BACH1 (pro-metastatic)
      - ZEB1 (pro)
      - SNAI1 (pro)
      - TWIST as TWIST1 (pro)
      - OCT4 as POU5F1 (pro)
      - LIN28 as LIN28A (pro)
      - RKIP as PEBP1 (anti)
    
### Bulk seq data compilation:
  - Randomly select 80 patients from the clinical for seq data; form/download/prepare a GDC query (via TCGABiolinks functionalities)
      - (Query data is provided in a SummarizedExperiment object; however, it can't directly be used for analysis for reasons below)
 
  - Here, we closely follow the methodology taken by [Wang, Hicks & Hansen (2020)](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009954), with the aim of first determining the presence of the mean-correlation bias in the bulk seq data, and then applying spatial quantile normalization (SpQN) to remove it if necessary
  - The authors choose log2(RPKM) as the basic normalization method for the read counts, so we need gene (exon) length data. However, length data is not contained within any attribute of the SummarizedExperiment object 
  - So, obtain the data from the reference genome annotation file, and read it into R as well
  - Subset both the annotation file and seq data such that they:
      - only contain protein-coding and miRNA-coding genes whose information are also found in both files
      - contain no duplicate gene names
  - Construct a new counts matrix and row data from these subset conditions, then coerce them into a new SummarizedExperiment object and preserve as a .rda file

### Bulk seq data subsetting:
TO DO: determine number of principal components to be removed, then randomly select 3993 other genes for the corr matrix for a total of 4000 genes (when added to the genes of interest)

### Bulk seq data corr(elation) analysis:
TO DO: plot the diagonal submatrices of the corr matrix to determine the presence of the mean-correlation bias, apply SpQN if necessary. then extract correlations of interest and visualize using heatmap. anything else?

## Single-cell RNA-seq data from the GEO: 
[TBD]



