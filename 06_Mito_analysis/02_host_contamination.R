# module load R/3.6.0

library(tidyverse)
library(gridExtra)
setwd("/ssd3/Mar_genome_analysis/bwa_mapping/mito/all_samples/somatypus")
setwd("/ssd3/Mar_genome_analysis/bwa_mapping/mito/all_samples_old/somatypus")
snvs <- read.delim("Somatypus_SNVs_final.counts", header = TRUE)
head(snvs)

snvs <- filter(snvs, pos < 12100 | pos > 12900) # exclude dloop region

# plotting functions
mt_scatter <- function(data, sample, name){
	ggplot(data, aes(get(paste0("C",sample,"_f")),get(paste0("T",sample,"_f"))))+
		geom_point()+
		ggtitle(name)+
  		theme_classic() +
  		theme(axis.text=element_text(size=12,face="bold"),
        		axis.title=element_text(size=16,face="bold"),
        		text=element_text(size=18,face="bold")) +
		xlab("Hemolymph SNV freq") +
		ylab("Tissue SNV freq")
}
mt_hist_lower <- function(data, sample, name){
	filter(data, get(paste0("T",sample,"_f")) > 0.1) %>%
	ggplot(aes(get(paste0("C",sample,"_f"))))+
		geom_histogram(binwidth = 0.005)+
		xlim(-0.01,0.2) +
		ggtitle(name) +
  		theme_classic() +
  		theme(axis.text=element_text(size=12,face="bold"),
        		axis.title=element_text(size=16,face="bold"),
        		text=element_text(size=18,face="bold")) +
		xlab("Hemolymph SNV freq")
}
mt_hist_upper <- function(data, sample, name){
	filter(data, get(paste0("T",sample,"_f")) < 0.9) %>%
	ggplot(aes(get(paste0("C",sample,"_f"))))+
		geom_histogram(binwidth = 0.005)+
		xlim(0.75,1.01) +
		ggtitle(name) +
  		theme_classic() +
  		theme(axis.text=element_text(size=12,face="bold"),
        		axis.title=element_text(size=16,face="bold"),
        		text=element_text(size=18,face="bold")) +
		xlab("Hemolymph SNV freq")
}

p <- list()
p[[1]] <- mt_scatter(snvs, "pei1","PEI-DN03")
p[[2]] <- mt_scatter(snvs, "pei2","PEI-DN07")
p[[3]] <- mt_scatter(snvs, "pei3","PEI-DN08")
p[[4]] <- mt_scatter(snvs, "usa1","FFM-22A10")
p[[5]] <- mt_scatter(snvs, "usa2","FFM-22F10")
p[[6]] <- mt_scatter(snvs, "usa3","MELC-A10")
p[[7]] <- mt_scatter(snvs, "usa4","MELC-A11")
p[[8]] <- mt_scatter(snvs, "usa5","NYTC-C9")

q <- list()
q[[1]] <- mt_hist_lower(snvs, "pei1","PEI-DN03")
q[[2]] <- mt_hist_lower(snvs, "pei2","PEI-DN07")
q[[3]] <- mt_hist_lower(snvs, "pei3","PEI-DN08")
q[[4]] <- mt_hist_lower(snvs, "usa1","FFM-22A10")
q[[5]] <- mt_hist_lower(snvs, "usa2","FFM-22F10")
q[[6]] <- mt_hist_lower(snvs, "usa3","MELC-A10")
q[[7]] <- mt_hist_lower(snvs, "usa4","MELC-A11")
q[[8]] <- mt_hist_lower(snvs, "usa5","NYTC-C9")

s <- list()
s[[1]] <- mt_hist_upper(snvs, "pei1","PEI-DN03")
s[[2]] <- mt_hist_upper(snvs, "pei2","PEI-DN07")
s[[3]] <- mt_hist_upper(snvs, "pei3","PEI-DN08")
s[[4]] <- mt_hist_upper(snvs, "usa1","FFM-22A10")
s[[5]] <- mt_hist_upper(snvs, "usa2","FFM-22F10")
s[[6]] <- mt_hist_upper(snvs, "usa3","MELC-A10")
s[[7]] <- mt_hist_upper(snvs, "usa4","MELC-A11")
s[[8]] <- mt_hist_upper(snvs, "usa5","NYTC-C9")

pdf("mito_snp_freq.pdf", height=10, width=10)
	do.call(grid.arrange,p)
	do.call(grid.arrange,q)
	do.call(grid.arrange,s)
dev.off()
