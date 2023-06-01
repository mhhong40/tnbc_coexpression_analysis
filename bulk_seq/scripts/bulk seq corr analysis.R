library(tidyverse)
library(here)
library(SummarizedExperiment)
library(spqn)

setwd(here::here("tnbc_bulk_rna_seq_data"))

# same for PC-removed data, obtain soon
# load("tnbc.Xpcs.rda")
# cor_matX <- cor(t(assay(tnbc.Xpcs)))load("tnbc.0pcs.rda")
# aveLogrpkmX <- rowData(tnbc.Xpcs)$aveLogrpkm

cor_mat0 <- cor(t(assay(tnbc.0pcs)))  # correlation matrix (Pearson r)
aveLogrpkm0 <- rowData(tnbc.0pcs)$aveLogrpkm

###############################################################################

# check how correlation distributions differ across different expression levels
plot_signal_condition_exp(cor_mat0, aveLogrpkm0, signal=0)

plot_signal_condition_exp(cor_mat0, aveLogrpkm0, signal=0.001)  # haven't tested this out yet.  

# apply spqn - in any case the mean-correlation relationship is present
cor_m0_spqn <- normalize_correlation(cor_mat0, ave_exp = aveLogrpkm, ngrp = x, size_grp = y, ref_grp = )

plot_signal_condition_exp(cor_m0_spqn, ave_logrpkm, signal=0)
plot_signal_condition_exp(cor_m0_spqn, ave_logrpkm, signal=0.001)

# TO DO: form correlation submatrices only with genes of interest
# and a heatmap for visuals
# anything else?
