library(matrixStats)
library(SummarizedExperiment)

load("fullTnbc.rda")

# log2(RPKM) values
cSums <- colSums(assay(fullTnbc))
log_rpkm <- sweep(log2(assay(fullTnbc) + 0.5), 2, FUN = "-", STATS = log2(cSums / 10^6)) 
log_rpkm <- log_rpkm - log2(rowData(fullTnbc)$exon_length / 1000)  

unExp <- c("LIN28A", "POU5F1")
unexpOnly <- as.data.frame(t(log_rpkm[rownames(log_rpkm) %in% unExp, ]))

print(paste(sum(unexpOnly$LIN28A > 0), "samples express LIN28")) # yields 1 
print(paste(sum(unexpOnly$POU5F1 > 0), "samples express OCT4"))  # yields 29