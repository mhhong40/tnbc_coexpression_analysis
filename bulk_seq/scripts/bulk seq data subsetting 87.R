# skip if packages already installed
BiocManager::install(c("WGCNA", "sva")) # * had issues installing GO.db dependency for WGCNA?

library(WGCNA)  
library(SummarizedExperiment)
library(sva)
library(matrixStats)

# assume setwd(here::here("tnbc_bulk_seq_data")) has already been called in one of the other scripts

# expression matrix data 
load("fullTnbc.rda")

# manually calculate log2(RPKM)s
cSums <- colSums(assay(fullTnbc))
log_rpkm <- sweep(log2(assay(fullTnbc) + 0.5), 2, FUN = "-", STATS = log2(cSums / 10^6)) 
log_rpkm <- log_rpkm - log2(rowData(fullTnbc)$exon_length / 1000)  

# keep only genes are above a certain expression lvl threshold
expressed <- which(rowMedians(log_rpkm) > 0)

# possible TO-DO: which percentage of genes are protein-coding vs miRNA-coding?
tnbc.0pcs <- fullTnbc[expressed, ]
logrpkm.0pcs <- log_rpkm[expressed, ]
aveLogrpkm <- rowMeans(logrpkm.0pcs)

# mean-centering, then variance scaling
# this converts each log2(RPKM) to its z-score according to the 
# expression profile distribution for the respective gene
logrpkm.0pcs <- (logrpkm.0pcs - aveLogrpkm) / matrixStats::rowSds(logrpkm.0pcs) 

assays(tnbc.0pcs) <- SimpleList(logrpkm = logrpkm.0pcs)
rowData(tnbc.0pcs)$aveLogrpkm <- aveLogrpkm
save(tnbc.0pcs, file = "tnbc.0pcs.rda", compress = "xz")


# removing principal components (PCs) is recommended to reduce (possible) batch effect
# estimate the # of surrogate variables using sva package
# (open q: exactly how many PCs to remove - may revisit later?)
tnbc.Xpcs <- tnbc.0pcs
tnbc.4pcs <- tnbc.0pcs # also try removal of 4 PCs (a la spqn authors)

dat <- t(SummarizedExperiment::assay(tnbc.Xpcs))
mod <- matrix(1, nrow = dim(dat)[1], ncol = 1)  # need to double check
colnames(mod) <- "Intercept"
numPC <- num.sv(t(dat), mod, method = "be") # yields 18

assay(tnbc.Xpcs) <- removePrincipalComponents(t(scale(t(logrpkm.0pcs))), n = numPC) 
save(tnbc.Xpcs, file = "tnbc.Xpcs.rda", compress = "xz")

assay(tnbc.4pcs) <- removePrincipalComponents(t(scale(t(logrpkm.0pcs))), n = 4) 
save(tnbc.4pcs, file = "tnbc.4pcs.rda", compress = "xz")


# besides our genes in the neighborhood of BACH1...
# randomly select 3995 genes for a total of 4000 (a la spqn authors)
# we can expand the sample size if appropriate
genesInterest <- c("BACH1", "ZEB1", "SNAI1", "LIN28A", "PEBP1", "POU5F1", "TWIST1")

set.seed(39872)
others <- which(!rownames(tnbc.0pcs) %in% genesInterest)

# no PCs removed, temporarily using these
tnbc.4k0 <- tnbc.0pcs[sample(others, size = 3995, replace = FALSE), ]
tnbc.4k0 <- append(tnbc.4k0, tnbc.0pcs[rownames(tnbc.0pcs) %in% genesInterest, ]) 
save(tnbc.4k0, file = "tnbc.4k0.rda", compress = "xz")

# these should have PCs removed:
## numPCs (18) as estimated by sva
tnbc.4kX <- tnbc.Xpcs[rownames(tnbc.4k0), ]
save(tnbc.4kX, file = "tnbc.4kX.rda", compress = "xz")

## 4 PCs
tnbc.4k4 <- tnbc.4pcs[rownames(tnbc.4k0), ]
save(tnbc.4k4, file = "tnbc.4k4.rda", compress = "xz")
