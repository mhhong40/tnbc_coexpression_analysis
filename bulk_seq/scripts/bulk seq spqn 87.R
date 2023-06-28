library(tidyverse)
library(here)
library(reshape)
library(SummarizedExperiment)
library(spqn)

# setwd(here::here("tnbc_bulk_rna_seq_data"))
load("tnbc.4k0.rda")
load("tnbc.4k4.rda")
load("tnbc.4kX.rda")

cor_mat0 <- cor(t(assay(tnbc.4k0)))  # correlation matrices (Pearson r)
aveLogrpkm0 <- rowData(tnbc.4k0)$aveLogrpkm

# we'll obtain the results for correlations with different # PCs removed
# for further exploratory analysis
cor_mat4 <- cor(t(assay(tnbc.4k4)))
aveLogrpkm4 <- rowData(tnbc.4k4)$aveLogrpkm

cor_matX <- cor(t(assay(tnbc.4kX)))
aveLogrpkmX <- rowData(tnbc.4kX)$aveLogrpkm

# check how correlation distributions differ across different expression levels
plot_signal_condition_exp(cor_mat0, aveLogrpkm0, signal=0)
plot_signal_condition_exp(cor_mat0, aveLogrpkm0, signal=0.001)  

IQRlist <- get_IQR_condition_exp(cor_mat0, aveLogrpkm0) # 2d boxplot
plot_IQR_condition_exp(IQRlist)

# apply spqn since mean-correlation relationship is present
# condense code?
cor_m0_spqn <- normalize_correlation(cor_mat0, ave_exp = aveLogrpkm0, ngrp = 10, size_grp = 15, ref_grp = 9)
cor_m4_spqn <- normalize_correlation(cor_mat4, ave_exp = aveLogrpkm0, ngrp = 10, size_grp = 15, ref_grp = 9)
cor_mX_spqn <- normalize_correlation(cor_matX, ave_exp = aveLogrpkm0, ngrp = 10, size_grp = 15, ref_grp = 9)

plot_signal_condition_exp(cor_m0_spqn, aveLogrpkm0, signal=0.001)

plot_signal_condition_exp(cor_m4_spqn, aveLogrpkm4, signal=0.001)

plot_signal_condition_exp(cor_mX_spqn, aveLogrpkmX, signal=0.001)

# obtain correlations for genes of interest
genesInterest <- c("BACH1", "ZEB1", "SNAI1", "LIN28A", "PEBP1", "POU5F1", "TWIST1")
isExpr <- genesInterest[genesInterest %in% rownames(tnbc.4kX)]  

# no spqn correlation subsets
# 0 PCs removed
noSpQN0 <- cor_mat0[rownames(cor_mat0) %in% isExpr, colnames(cor_mat0) %in% isExpr]
rownames(noSpQN0) <- isExpr
colnames(noSpQN0) <- isExpr
save(noSpQN0, file = "noSpQN 87.rda")

# 4 PCs removed
noSpQN4 <- cor_mat4[rownames(cor_mat0) %in% isExpr, colnames(cor_mat0) %in% isExpr]
rownames(noSpQN4) <- isExpr
colnames(noSpQN4) <- isExpr
save(noSpQN4, file = "noSpQN4 87.rda")

# 18 PCs removed
noSpQNX <- cor_matX[rownames(cor_mat0) %in% isExpr, colnames(cor_mat0) %in% isExpr]
rownames(noSpQNX) <- isExpr
colnames(noSpQNX) <- isExpr
save(noSpQNX, file = "noSpQNX 87.rda")

# subsets with pcs removed
# 0 pcs
b0Subset <- cor_m0_spqn[rownames(cor_mat0) %in% isExpr, colnames(cor_mat0) %in% isExpr]
rownames(b0Subset) <- isExpr
colnames(b0Subset) <- isExpr
save(b0Subset, file = "b0Subset.rda")

# 4 pcs
b4Subset <- cor_m4_spqn[rownames(cor_mat0) %in% isExpr, colnames(cor_mat0) %in% isExpr]
rownames(b4Subset) <- isExpr
colnames(b4Subset) <- isExpr
save(b4Subset, file = "b4Subset.rda")

# 18 pcs
bXSubset <- cor_mX_spqn[rownames(cor_mat0) %in% isExpr, colnames(cor_mat0) %in% isExpr]
rownames(bXSubset) <- isExpr
colnames(bXSubset) <- isExpr
save(bXSubset, file = "bXSubset.rda")


# heatmap visualizations of Pearson r correlation matrices
b0SLong <- melt(b0Subset, as.is = TRUE)
hmap0 <- ggplot(b0SLong, aes(X1, X2)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue")
print(hmap0)

b4SLong <- melt(b4Subset)
hmap4 <- ggplot(b4SLong, aes(X1, X2)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue")
print(hmap4)

bXSLong <- melt(bXSubset)
hmapX <- ggplot(bXSLong, aes(X1, X2)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue")
print(hmapX)
