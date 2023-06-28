# TNBC RNA seq analysis
'23 spring -> present research with Dr. Gabor Balazsi @ Stony Brook University

## Bulk RNA-seq data from the GDC Repository:
`bulk seq data compilation.R` : contains code for the GDC bulk RNA-seq data query (87 samples) as well as encapsulation of raw counts, patient data, and gene lengths from the [GRCh38 reference genome](https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files) within a SummarizedExperiment object.

`bulk seq data subsetting 87.R` : here we closely follow the sample data preprocessing steps taken by [Wang et al. (2020)](https://s3.jcloud.sjtu.edu.cn/899a892efef34b1b944a19981040f55b-oss01/bioconductor/3.13/bioc/vignettes/spqn/inst/doc/spqn.html). includes log2(RPKM) calculations, conversion to z-scores, removal of principal components (PCs), and subsetting of data to our genes of interest and several others for a total of 4000 genes.

`bulk seq spqn 87.R` : includes the calculation of the three correlation matrices (Pearson r) and ridge plots comparing correlation distributions across expression levels, pre and post-SpQN. also includes subsets of the pre-SpQN and post-SpQN correlation matrices containing only our genes of interest, and heatmap visualizations of the latter.

`bulk seq corr stats 87.R` : for evaluation of the Pearson r as the correlation coefficient given the dataset with 0 PCs removed. includes Q-Q plots, Shapiro-Wilk tests, and scatter plots.

`bulk seq unexpressed 87.R` : exploratory analysis of the unexpressed genes of interest, LIN28 and OCT4 (as POU5F1).

## Single-cell RNA-seq data from the GEO: 
[TBD]



