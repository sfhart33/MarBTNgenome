library(tidyverse)
setwd("/ssd3/Mar_genome_analysis/dnds/dndscv/outputs")
dnds_files <- list.files(pattern = ".vcf.dnds.rds")

start <- readRDS("USA_1234_noP_noH_nopeiLOH.vcf.dnds.rds")
dnds_summary <- data.frame(start$globaldnds[1,], stringsAsFactors = FALSE)%>%
	mutate(name = "delete_me")


for(dnds_file in dnds_files){
	dnds_name <- str_split(dnds_file,".vcf")[[1]][1]
	print(paste("Loading",dnds_name), quote = FALSE)
	dnds_output <- readRDS(dnds_file)
	dnds_value <- data.frame(dnds_output$globaldnds[1,], stringsAsFactors = FALSE) %>%
		mutate(name = dnds_name)
	dnds_summary <- rbind(dnds_summary,dnds_value)
	print("     done", quote = FALSE)	
}
dnds_summary <- dnds_summary[-1,]
rm(dnds_output)
print(dnds_summary)

saveRDS(dnds_summary, "dnds_global_summary.rds")