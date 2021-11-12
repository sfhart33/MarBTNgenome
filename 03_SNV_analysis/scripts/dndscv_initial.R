# module load R
# R
library(devtools)
install_github("im3sanger/dndscv")
library(dndscv)
library(tidyverse)

buildref(cdsfile="/ssd3/Mar_genome_analysis/dnds/2020-09-11_Mar_genome_snap02.all.dNdScv.input",
         genomefile="/ssd3/Mar_genome_analysis/genomes/Mar.3.4.6.p1_Q30Q30A.fasta",
         outfile = "/ssd3/Mar_genome_analysis/dnds/dndscv/Mar.3.4.6.p1_snap02_refCDS.rda")