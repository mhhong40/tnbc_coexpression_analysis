# TNBC-RNA-seq-analysis
'23 spring -> present research with Dr. Gabor Balazsi @ Stony Brook University

[DOCUMENTATION]

## Handling bulk RNA-seq data from the GDC Repository:
  - Retrieve TCGA-BRCA clinical data from [the NIH GDC Repository](https://portal.gdc.cancer.gov/repository?filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22content%22%3A%7B%22field%22%3A%22cases.project.project_id%22%2C%22value%22%3A%5B%22TCGA-BRCA%22%5D%7D%2C%22op%22%3A%22in%22%7D%2C%7B%22content%22%3A%7B%22field%22%3A%22files.cases.primary_site%22%2C%22value%22%3A%5B%22breast%22%5D%7D%2C%22op%22%3A%22in%22%7D%2C%7B%22content%22%3A%7B%22field%22%3A%22files.data_type%22%2C%22value%22%3A%5B%22Gene%20Expression%20Quantification%22%5D%7D%2C%22op%22%3A%22in%22%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.experimental_strategy%22%2C%22value%22%3A%5B%22RNA-Seq%22%5D%7D%7D%5D%7D&searchTableTab=files)
  - Transfer clinical .txt file data into workable Excel spreadsheet
  - Sort by IHC results to select for just TNBC (total: n = 115 patients)
  
  - List genes of interest:
      - BACH1 (pro-metastatic)
      - ZEB1 (pro)
      - SNAI1 (pro)
      - TWIST as TWIST1 (pro)
      - OCT4 as POU5F1 (pro)
      - LIN28 as LIN28A (pro)
      - RKIP as PEBP1 (anti)
      
  - Construct expression matrix in R (start with first 10 patients in spreadsheet (bolded names))
      - Download STAR Augmented Gene Counts (.tsv) files for each patient. Some manual file preprocessing:
          - Rename each file to match patient ID (for later convenience)
          - Remove first line (i.e., "# gene-model: GENCODE v36") of each tsv file so R processes it with the correct number of columns
      - Final result: gene names as rows, patients as columns, unstranded read counts* as matrix values

*[The GDC Docs](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/) recommends that users normalize raw read count values if a subset of genes is investigated (as in our case)
  - GDC treats reads as unstranded for all analyses for uniformity
  - Since we are not performing a task where transcript directionality is important, unstranded reads should be satisfactory for us as well
  - However, if it is appropriate, perhaps we may just use the normalized values (I am getting mixed results from searching how to handle gene subset RNA-seq data)?

## Single-cell RNA-seq data from the GEO: 
[TBD]



