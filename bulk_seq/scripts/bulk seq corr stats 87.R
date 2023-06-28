library(tidyverse)
library(here)
library(SummarizedExperiment)

load("tnbc.4k0.rda")

genesInterest <- c("BACH1", "ZEB1", "SNAI1", "LIN28A", "PEBP1", "POU5F1", "TWIST1")
isExpr <- genesInterest[genesInterest %in% rownames(tnbc.4k0)]

subsetExpr0 <- assay(tnbc.4k0)[rownames(assay(tnbc.4k0)) %in% isExpr, ]
subsetExpr0 <- subsetExpr0[order(match(rownames(subsetExpr0), isExpr)), , drop = FALSE]

# Q-Q plots
for (i in 1:5) {
  qqnorm(subsetExpr0[i, ], 
         main = base::paste(isExpr[i], "Expression Profile"))
  qqline(subsetExpr0[i, ])
}

# Shapiro-Wilk tests
for (i in 1:5) {
  print(shapiro.test(subsetExpr0[i, ]))
}

# scatter plots and corrrelation p-values for each gene pair
for (i in 1:4) {
  for (j in 5:(i+1)) {
    print(cor.test(subsetExpr0[i, ], subsetExpr0[j, ]))
    xData <- as.data.frame(subsetExpr0[i, ])
    yData <- as.data.frame(subsetExpr0[j, ])
    scatterData <- base::cbind(xData, yData)
    print(ggplot(data = scatterData, aes(x = subsetExpr0[i, ], y = subsetExpr0[j, ])) + 
      geom_point(color = 'blue') +
      labs(x = isExpr[i], y = isExpr[j]))
  }
}
