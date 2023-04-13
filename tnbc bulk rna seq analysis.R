# install packages
install.packages("readr")
install.packages("tidyverse")
install.packages("tidyr")
install.packages("dplyr")
install.packages("data.table")

# load libraries
library(readr)
library(tidyverse)
library(tidyr)
library(dplyr)
library(data.table)

# function to read in a seq data tsv and add the file name (patient ID) as a column
customized_read_tsv <- function(file){
    	read_tsv(file) %>%
        	mutate(fileName = substr(file, start = 3, stop = 14)) # extract just patient ID
}

# temporary expression matrix w/ all genes
temp <- list.files(pattern = ".tsv", full.names = TRUE) %>% # list all the files
	lapply(customized_read_tsv) %>% # read them all in w/ custom function
  reduce(bind_rows) %>% # stack them all on top of each other
  select(gene_name, fileName, unstranded) # list gene name, patient ID, unstranded read count

# create final w/ only genes of interest
# also: unique identifier row for each patient ID so unstranded values appear properly
genes_interest <- c('BACH1', 'ZEB1', 'SNAI1', 'TWIST1', 'POU5F1', 'LIN28A', 'PEBP1')

bulk_exp_matrix <- subset(temp, gene_name %in% genes_interest) %>%
  group_by(fileName) %>%
  mutate(row = row_number()) %>%	
  pivot_wider(names_from = fileName, values_from = unstranded) %>%
  select(-row)

# TODO: normalize unstranded read counts appropriately 