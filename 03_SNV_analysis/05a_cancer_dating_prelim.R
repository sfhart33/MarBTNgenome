# Original file here: C:\Users\shart\Metzger Lab Dropbox\Sam_data\Mya_genome\SNPs\somatypus_output_dating.r

library(tidyverse)
library(ape)
library(sigfit)
library(lsa)
data(cosmic_signatures_v2)
data(cosmic_signatures_v3)

setwd("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run")


# Load data if running on all SNVs
    # snvs.nr <- read.table("Somatypus_SNVs_final_NR.txt", header=T, check.names=F)
    # snvs.nv <- read.table("Somatypus_SNVs_final_NV.txt", header=T, check.names=F)
    # snvs.metadata <- read.table("Somatypus_SNVs_final_Metadata.txt", header=T)
    # setwd("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/pairwise")  
# Load data if running on SNVs outside possible LOH regions
    snvs.nr <- read.table("Somatypus_SNVs_final_noLOH2_NR.txt", header=T, check.names=F)
    snvs.nv <- read.table("Somatypus_SNVs_final_noLOH2_NV.txt", header=T, check.names=F)
    snvs.metadata <- read.table("Somatypus_SNVs_final_noLOH2_Metadata.txt", header=T)
    setwd("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/pairwise/noLOH")

# Define sample names and combine into one table
    samples <- c("H_ref","H_usa","H_pei",
    "C_pei1","C_pei2","C_pei4","C_pei3", # C_pei4 has high host contamination
    "C_usa1","C_usa2","C_usa6","C_usa3","C_usa7","C_usa4","C_usa5", # C_usa6 & 7 have high host contamination
    "C_pei2T","C_pei4T","C_pei3T","C_usa6T","C_usa3T","C_usa7T","C_usa4T","C_usa5T")
    samplesNR <- paste0(samples, "_NR")
    samplesNV <- paste0(samples, "_NV")
    colnames(snvs.nr) <- samplesNR
    colnames(snvs.nv) <- samplesNV
    colnames(indels.nr) <- samplesNR
    colnames(indels.nv) <- samplesNV
    snvs.metadata2 <- snvs.metadata %>%
        select(CHROM,POS,REF,ALT)
    snvs <- cbind(snvs.metadata2,snvs.nv,snvs.nr)

# Calculate average read depth across the genome
    for(sample in samples){
        depth <- mean(snvs[,paste0(sample,"_NR")])
        assign(paste0(sample,"_avg"), depth)
    }
    for(sample in samples){
        print(sample)
        print(get(paste0(sample,"_avg")))
    }

#Set thresholds
    min_cutoff=3
    stringent_depth=6 # was 8
    weak_depth=16

# Not found in any healthy, found in at least one cancer above stringent threshold
    realSNVs <- filter(snvs, !(H_ref_NV > min_cutoff |
                        H_usa_NV > min_cutoff |
                        H_pei_NV > min_cutoff) &
                        (C_pei1_NV > C_pei1_avg/stringent_depth |
                        C_pei2_NV > C_pei2_avg/stringent_depth |
                        C_pei3_NV > C_pei3_avg/stringent_depth |
                        C_usa1_NV > C_usa1_avg/stringent_depth |
                        C_usa2_NV > C_usa2_avg/stringent_depth |
                        C_usa3_NV > C_usa3_avg/stringent_depth |
                        C_usa4_NV > C_usa4_avg/stringent_depth |
                        C_usa5_NV > C_usa5_avg/stringent_depth))

################# Old notes in original file


# Counts since divergence to get time of divergence estimate
    noPEI <- filter(realSNVs,
                        !(C_pei1_NV > C_pei1_avg/weak_depth |
                        C_pei2_NV > C_pei2_avg/weak_depth |
                        C_pei3_NV > C_pei3_avg/weak_depth))
    usa1 <- filter(noPEI,
                        (C_usa1_NV > C_usa1_avg/stringent_depth))
    usa2 <- filter(noPEI,
                        (C_usa2_NV > C_usa2_avg/stringent_depth))
    usa3 <- filter(noPEI,
                        (C_usa3_NV > C_usa3_avg/stringent_depth))
    usa4 <- filter(noPEI,
                        (C_usa4_NV > C_usa4_avg/stringent_depth))
    usa5 <- filter(noPEI,
                        (C_usa5_NV > C_usa5_avg/stringent_depth))    
    noUSA <- filter(realSNVs,
                        !(C_usa1_NV > C_usa1_avg/weak_depth |
                        C_usa2_NV > C_usa2_avg/weak_depth |
                        C_usa3_NV > C_usa3_avg/weak_depth |
                        C_usa4_NV > C_usa4_avg/weak_depth |
                        C_usa5_NV > C_usa5_avg/weak_depth))
    pei1 <- filter(noUSA,
                        (C_pei1_NV > C_pei1_avg/stringent_depth))
    pei2 <- filter(noUSA,
                        (C_pei2_NV > C_pei2_avg/stringent_depth))
    pei3 <- filter(noUSA,
                        (C_pei3_NV > C_pei3_avg/stringent_depth))
# # counts since sublineage divergence
#     notallUSA <- filter(noPEI,
#                         !(C_usa1_NV > C_usa1_avg/weak_depth) |
#                         !(C_usa2_NV > C_usa2_avg/weak_depth) |
#                         !(C_usa3_NV > C_usa3_avg/weak_depth) |
#                         !(C_usa4_NV > C_usa4_avg/weak_depth) |
#                         !(C_usa5_NV > C_usa5_avg/weak_depth))
#     usa1b <- filter(notallUSA,
#                         (C_usa1_NV > C_usa1_avg/stringent_depth))
#     usa2b <- filter(notallUSA,
#                         (C_usa2_NV > C_usa2_avg/stringent_depth))
#     usa3b <- filter(notallUSA,
#                         (C_usa3_NV > C_usa3_avg/stringent_depth))
#     usa4b <- filter(notallUSA,
#                         (C_usa4_NV > C_usa4_avg/stringent_depth))
#     usa5b <- filter(notallUSA,
#                         (C_usa5_NV > C_usa5_avg/stringent_depth))    
#     notallPEI <- filter(noUSA,
#                         !(C_pei1_NV > C_pei1_avg/weak_depth) |
#                         !(C_pei2_NV > C_pei2_avg/weak_depth) |
#                         !(C_pei3_NV > C_pei3_avg/weak_depth))
#     pei1b <- filter(notallPEI,
#                         (C_pei1_NV > C_pei1_avg/stringent_depth))
#     pei2b <- filter(notallPEI,
#                         (C_pei2_NV > C_pei2_avg/stringent_depth))
#     pei3b <- filter(notallPEI,
#                         (C_pei3_NV > C_pei3_avg/stringent_depth))

# Counts check
    div_counts <- numeric(0)
    #div_counts <- append(div_counts,1:3)
    for(i in c("pei1","pei2","pei3","usa1","usa2","usa3","usa4","usa5")){
        count1 <- get(i) %>% nrow()
        #count2 <- get(paste0(i,"b")) %>% nrow()
        div_counts <- append(div_counts,count1)
        #div_counts <- append(div_counts,count2)
    }
    div_counts

# Histograms check
    plotalt2 <- function(snpfile, sampleID){
        plot1 <- ggplot(get(snpfile), aes(get(paste0(sampleID,"_NV")))) +
            geom_histogram(binwidth=1)+
            xlim(0,100)+
            xlab("Allele count")+
            theme_classic() +
                theme(axis.text=element_text(size=12,face="bold"),
                axis.title=element_text(size=16,face="bold"),
                text=element_text(size=18,face="bold")) +
            ggtitle(paste0(snpfile," ", sampleID, " SNVs count: ", nrow(snpfile)))
        return(plot1)
    }
    pdf("post_divergence_histograms.pdf")
    for(i in c("usa1","usa2","usa3","usa4","usa5","pei1","pei2","pei3")){
        plotalt2(i,paste0("C_",i)) %>% print()
        plotalt2(paste0(i,"b"),paste0("C_",i)) %>% print()
    }
    dev.off()

# Build helmsman input
    helmsmanfile_div <- function(sample){
        name1 <- paste0(sample,"_div1")
        get(sample) %>% select(CHROM, POS, REF, ALT) %>%
            mutate(ID = name1) %>%
            write.table(file = paste0("helmsman.divergence.txt"), append = TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep = '\t')
        # name2 <- paste0(sample,"_div2")
        # get(paste0(sample,"b")) %>% select(CHROM, POS, REF, ALT) %>%
        #     mutate(ID = name2) %>%
        #     write.table(file = paste0("helmsman.divergence.txt"), append = TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep = '\t')
    }
    if(file.exists("helmsman.divergence.txt")) {
        file.remove("helmsman.divergence.txt")
    }
    for(i in c("usa1","usa2","usa3","usa4","usa5","pei1","pei2","pei3")){
        helmsmanfile_div(i)
    }
