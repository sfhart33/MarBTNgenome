# usage: 
# dndscv_script.R filename location

library(dndscv)
# library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
location <- args[2]

mutations_list <- read.table(paste(location,filename,sep="/"))

dndsout <- dndscv(mutations_list,
                  refdb="/ssd3/Mar_genome_analysis/dnds/dndscv/Mar.3.4.6.p1_snap02_refCDS.rda",
                  cv=NULL,
                  #outp = 1,
                  max_coding_muts_per_sample = 1000000000,
                  max_muts_per_gene_per_sample = 1000000000)

saveRDS(dndsout, file = paste(filename,".rds", sep = ""))