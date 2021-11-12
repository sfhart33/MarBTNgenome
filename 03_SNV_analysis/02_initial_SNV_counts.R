# From original file: C:\Users\shart\Metzger Lab Dropbox\Sam_data\Mya_genome\SNPs\somatypus_output_initial_counts.r

library(tidyverse)


# Load full genome data (same as when loading for signature extraction)
    setwd("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run")
    snvs.nr <- read.table("Somatypus_SNVs_final_NR.txt", header=T, check.names=F)
    snvs.nv <- read.table("Somatypus_SNVs_final_NV.txt", header=T, check.names=F)
    snvs.metadata <- read.table("Somatypus_SNVs_final_Metadata.txt", header=T)
    samples <- c("H_ref","H_usa","H_pei",
    "C_pei1","C_pei2","C_pei4","C_pei3", # C_pei4 has high host contamination
    "C_usa1","C_usa2","C_usa6","C_usa3","C_usa7","C_usa4","C_usa5", # C_usa6 & 7 have high host contamination
    "C_pei2T","C_pei4T","C_pei3T","C_usa6T","C_usa3T","C_usa7T","C_usa4T","C_usa5T")
    samplesNR <- paste0(samples, "_NR")
    samplesNV <- paste0(samples, "_NV")
    colnames(snvs.nr) <- samplesNR
    colnames(snvs.nv) <- samplesNV
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

setwd("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/bins")
    min_cutoff=3
    stringent_depth=6 # was 8
    weak_depth=16


################

#BINNING - same as for signature extraction
    # Is SNV found in any cancer
        anyC <- filter(snvs, C_pei1_NV > C_pei1_avg/stringent_depth |
                            C_pei2_NV > C_pei2_avg/stringent_depth |
                            C_pei3_NV > C_pei3_avg/stringent_depth |
                            C_usa1_NV > C_usa1_avg/stringent_depth |
                            C_usa2_NV > C_usa2_avg/stringent_depth |
                            C_usa3_NV > C_usa3_avg/stringent_depth |
                            C_usa4_NV > C_usa4_avg/stringent_depth |
                            C_usa5_NV > C_usa5_avg/stringent_depth)
        noC <- filter(snvs, !(C_pei1_NV > C_pei1_avg/stringent_depth |
                            C_pei2_NV > C_pei2_avg/stringent_depth |
                            C_pei3_NV > C_pei3_avg/stringent_depth |
                            C_usa1_NV > C_usa1_avg/stringent_depth |
                            C_usa2_NV > C_usa2_avg/stringent_depth |
                            C_usa3_NV > C_usa3_avg/stringent_depth |
                            C_usa4_NV > C_usa4_avg/stringent_depth |
                            C_usa5_NV > C_usa5_avg/stringent_depth))

    # Is SNV found in any healthy clam    
        anyHanyC <- filter(anyC, H_ref_NV > min_cutoff |
                            H_usa_NV > min_cutoff |
                            H_pei_NV > min_cutoff)
        noHanyC <- filter(anyC, !(H_ref_NV > min_cutoff |
                            H_usa_NV > min_cutoff |
                            H_pei_NV > min_cutoff))
        anyHnoC <- filter(noC, H_ref_NV > min_cutoff |
                            H_usa_NV > min_cutoff |
                            H_pei_NV > min_cutoff)
        noHnoC <- filter(noC, !(H_ref_NV > min_cutoff |
                            H_usa_NV > min_cutoff |
                            H_pei_NV > min_cutoff))
    # which healthy clam is it found in
        anyH <- rbind(anyHnoC, anyHanyC)  %>%
            arrange(CHROM,POS)
        H_R <- filter(noC, (H_ref_NV > min_cutoff) &
                            !(H_usa_NV > min_cutoff) &
                            !(H_pei_NV > min_cutoff))
        H_RU <- filter(noC, (H_ref_NV > min_cutoff) &
                            (H_usa_NV > min_cutoff) &
                            !(H_pei_NV > min_cutoff))
        H_RP <- filter(noC, (H_ref_NV > min_cutoff) &
                            !(H_usa_NV > min_cutoff) &
                            (H_pei_NV > min_cutoff))
        H_RUP <- filter(noC, (H_ref_NV > min_cutoff) &
                            (H_usa_NV > min_cutoff) &
                            (H_pei_NV > min_cutoff))
        H_U <- filter(noC, !(H_ref_NV > min_cutoff) &
                            (H_usa_NV > min_cutoff) &
                            !(H_pei_NV > min_cutoff))
        H_UP <- filter(noC, !(H_ref_NV > min_cutoff) &
                            (H_usa_NV > min_cutoff) &
                            (H_pei_NV > min_cutoff))      
        H_P <- filter(noC, !(H_ref_NV > min_cutoff) &
                            !(H_usa_NV > min_cutoff) &
                            (H_pei_NV > min_cutoff))                                            
    # Likely founder germline
        allCanyH <- filter(anyHanyC, C_pei1_NV > C_pei1_avg/weak_depth &
                            C_pei2_NV > C_pei2_avg/weak_depth &
                            C_pei3_NV > C_pei3_avg/weak_depth &
                            C_usa1_NV > C_usa1_avg/weak_depth &
                            C_usa2_NV > C_usa2_avg/weak_depth &
                            C_usa3_NV > C_usa3_avg/weak_depth &
                            C_usa4_NV > C_usa4_avg/weak_depth &
                            C_usa5_NV > C_usa5_avg/weak_depth)
    # Likely founder germline in LOH region or mismap
        subsetCanyH <- filter(anyHanyC, !(C_pei1_NV > C_pei1_avg/weak_depth &
                            C_pei2_NV > C_pei2_avg/weak_depth &
                            C_pei3_NV > C_pei3_avg/weak_depth &
                            C_usa1_NV > C_usa1_avg/weak_depth &
                            C_usa2_NV > C_usa2_avg/weak_depth &
                            C_usa3_NV > C_usa3_avg/weak_depth &
                            C_usa4_NV > C_usa4_avg/weak_depth &
                            C_usa5_NV > C_usa5_avg/weak_depth))
    # Founder germline OR early somatic mutation
        allCnoH <- filter(noHanyC, C_pei1_NV > C_pei1_avg/weak_depth &
                            C_pei2_NV > C_pei2_avg/weak_depth &
                            C_pei3_NV > C_pei3_avg/weak_depth &
                            C_usa1_NV > C_usa1_avg/weak_depth &
                            C_usa2_NV > C_usa2_avg/weak_depth &
                            C_usa3_NV > C_usa3_avg/weak_depth &
                            C_usa4_NV > C_usa4_avg/weak_depth &
                            C_usa5_NV > C_usa5_avg/weak_depth)
    # Somtic mutations OR mismaps, maybe some founder germline that dropped out of some samples
        subsetCnoH <- filter(noHanyC, !(C_pei1_NV > C_pei1_avg/weak_depth &
                            C_pei2_NV > C_pei2_avg/weak_depth &
                            C_pei3_NV > C_pei3_avg/weak_depth &
                            C_usa1_NV > C_usa1_avg/weak_depth &
                            C_usa2_NV > C_usa2_avg/weak_depth &
                            C_usa3_NV > C_usa3_avg/weak_depth &
                            C_usa4_NV > C_usa4_avg/weak_depth &
                            C_usa5_NV > C_usa5_avg/weak_depth))
    # Sublineage specific
        anyUSAnoPEInoH <- filter(subsetCnoH, !(C_pei1_NV > C_pei1_avg/weak_depth |
                            C_pei2_NV > C_pei2_avg/weak_depth |
                            C_pei3_NV > C_pei3_avg/weak_depth) &
                            (C_usa1_NV > C_usa1_avg/weak_depth |
                            C_usa2_NV > C_usa2_avg/weak_depth |
                            C_usa3_NV > C_usa3_avg/weak_depth |
                            C_usa4_NV > C_usa4_avg/weak_depth |
                            C_usa5_NV > C_usa5_avg/weak_depth))
        anyPEInoUSAnoH <- filter(subsetCnoH, (C_pei1_NV > C_pei1_avg/weak_depth |
                            C_pei2_NV > C_pei2_avg/weak_depth |
                            C_pei3_NV > C_pei3_avg/weak_depth) &
                            !(C_usa1_NV > C_usa1_avg/weak_depth |
                            C_usa2_NV > C_usa2_avg/weak_depth |
                            C_usa3_NV > C_usa3_avg/weak_depth |
                            C_usa4_NV > C_usa4_avg/weak_depth |
                            C_usa5_NV > C_usa5_avg/weak_depth))        
        allUSAnoPEInoH <- filter(subsetCnoH, !(C_pei1_NV > C_pei1_avg/weak_depth |
                            C_pei2_NV > C_pei2_avg/weak_depth |
                            C_pei3_NV > C_pei3_avg/weak_depth) &
                            (C_usa1_NV > C_usa1_avg/weak_depth &
                            C_usa2_NV > C_usa2_avg/weak_depth &
                            C_usa3_NV > C_usa3_avg/weak_depth &
                            C_usa4_NV > C_usa4_avg/weak_depth &
                            C_usa5_NV > C_usa5_avg/weak_depth))
        allPEInoUSAnoH <- filter(subsetCnoH, (C_pei1_NV > C_pei1_avg/weak_depth &
                            C_pei2_NV > C_pei2_avg/weak_depth &
                            C_pei3_NV > C_pei3_avg/weak_depth) &
                            !(C_usa1_NV > C_usa1_avg/weak_depth |
                            C_usa2_NV > C_usa2_avg/weak_depth |
                            C_usa3_NV > C_usa3_avg/weak_depth |
                            C_usa4_NV > C_usa4_avg/weak_depth |
                            C_usa5_NV > C_usa5_avg/weak_depth))

nrow(anyC)
nrow(noC)
nrow(anyHanyC)
nrow(noHanyC)
nrow(anyHnoC)
nrow(noHnoC)
nrow(anyH)
nrow(allCanyH)
nrow(allCnoH)# From original file: C:\Users\shart\Metzger Lab Dropbox\Sam_data\Mya_genome\SNPs\somatypus_output_initial_counts.r

library(tidyverse)


# Load full genome data (same as when loading for signature extraction)
    setwd("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run")
    snvs.nr <- read.table("Somatypus_SNVs_final_NR.txt", header=T, check.names=F)
    snvs.nv <- read.table("Somatypus_SNVs_final_NV.txt", header=T, check.names=F)
    snvs.metadata <- read.table("Somatypus_SNVs_final_Metadata.txt", header=T)
    samples <- c("H_ref","H_usa","H_pei",
    "C_pei1","C_pei2","C_pei4","C_pei3", # C_pei4 has high host contamination
    "C_usa1","C_usa2","C_usa6","C_usa3","C_usa7","C_usa4","C_usa5", # C_usa6 & 7 have high host contamination
    "C_pei2T","C_pei4T","C_pei3T","C_usa6T","C_usa3T","C_usa7T","C_usa4T","C_usa5T")
    samplesNR <- paste0(samples, "_NR")
    samplesNV <- paste0(samples, "_NV")
    colnames(snvs.nr) <- samplesNR
    colnames(snvs.nv) <- samplesNV
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

setwd("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/bins")
    min_cutoff=3
    stringent_depth=6 # was 8
    weak_depth=16


################

#BINNING - same as for signature extraction
    # Is SNV found in any cancer
        anyC <- filter(snvs, C_pei1_NV > C_pei1_avg/stringent_depth |
                            C_pei2_NV > C_pei2_avg/stringent_depth |
                            C_pei3_NV > C_pei3_avg/stringent_depth |
                            C_usa1_NV > C_usa1_avg/stringent_depth |
                            C_usa2_NV > C_usa2_avg/stringent_depth |
                            C_usa3_NV > C_usa3_avg/stringent_depth |
                            C_usa4_NV > C_usa4_avg/stringent_depth |
                            C_usa5_NV > C_usa5_avg/stringent_depth)
        noC <- filter(snvs, !(C_pei1_NV > C_pei1_avg/stringent_depth |
                            C_pei2_NV > C_pei2_avg/stringent_depth |
                            C_pei3_NV > C_pei3_avg/stringent_depth |
                            C_usa1_NV > C_usa1_avg/stringent_depth |
                            C_usa2_NV > C_usa2_avg/stringent_depth |
                            C_usa3_NV > C_usa3_avg/stringent_depth |
                            C_usa4_NV > C_usa4_avg/stringent_depth |
                            C_usa5_NV > C_usa5_avg/stringent_depth))

    # Is SNV found in any healthy clam    
        anyHanyC <- filter(anyC, H_ref_NV > min_cutoff |
                            H_usa_NV > min_cutoff |
                            H_pei_NV > min_cutoff)
        noHanyC <- filter(anyC, !(H_ref_NV > min_cutoff |
                            H_usa_NV > min_cutoff |
                            H_pei_NV > min_cutoff))
        anyHnoC <- filter(noC, H_ref_NV > min_cutoff |
                            H_usa_NV > min_cutoff |
                            H_pei_NV > min_cutoff)
        noHnoC <- filter(noC, !(H_ref_NV > min_cutoff |
                            H_usa_NV > min_cutoff |
                            H_pei_NV > min_cutoff))
    # which healthy clam is it found in
        anyH <- rbind(anyHnoC, anyHanyC)  %>%
            arrange(CHROM,POS)
        H_R <- filter(noC, (H_ref_NV > min_cutoff) &
                            !(H_usa_NV > min_cutoff) &
                            !(H_pei_NV > min_cutoff))
        H_RU <- filter(noC, (H_ref_NV > min_cutoff) &
                            (H_usa_NV > min_cutoff) &
                            !(H_pei_NV > min_cutoff))
        H_RP <- filter(noC, (H_ref_NV > min_cutoff) &
                            !(H_usa_NV > min_cutoff) &
                            (H_pei_NV > min_cutoff))
        H_RUP <- filter(noC, (H_ref_NV > min_cutoff) &
                            (H_usa_NV > min_cutoff) &
                            (H_pei_NV > min_cutoff))
        H_U <- filter(noC, !(H_ref_NV > min_cutoff) &
                            (H_usa_NV > min_cutoff) &
                            !(H_pei_NV > min_cutoff))
        H_UP <- filter(noC, !(H_ref_NV > min_cutoff) &
                            (H_usa_NV > min_cutoff) &
                            (H_pei_NV > min_cutoff))      
        H_P <- filter(noC, !(H_ref_NV > min_cutoff) &
                            !(H_usa_NV > min_cutoff) &
                            (H_pei_NV > min_cutoff))                                            
    # Likely founder germline
        allCanyH <- filter(anyHanyC, C_pei1_NV > C_pei1_avg/weak_depth &
                            C_pei2_NV > C_pei2_avg/weak_depth &
                            C_pei3_NV > C_pei3_avg/weak_depth &
                            C_usa1_NV > C_usa1_avg/weak_depth &
                            C_usa2_NV > C_usa2_avg/weak_depth &
                            C_usa3_NV > C_usa3_avg/weak_depth &
                            C_usa4_NV > C_usa4_avg/weak_depth &
                            C_usa5_NV > C_usa5_avg/weak_depth)
    # Likely founder germline in LOH region or mismap
        subsetCanyH <- filter(anyHanyC, !(C_pei1_NV > C_pei1_avg/weak_depth &
                            C_pei2_NV > C_pei2_avg/weak_depth &
                            C_pei3_NV > C_pei3_avg/weak_depth &
                            C_usa1_NV > C_usa1_avg/weak_depth &
                            C_usa2_NV > C_usa2_avg/weak_depth &
                            C_usa3_NV > C_usa3_avg/weak_depth &
                            C_usa4_NV > C_usa4_avg/weak_depth &
                            C_usa5_NV > C_usa5_avg/weak_depth))
    # Founder germline OR early somatic mutation
        allCnoH <- filter(noHanyC, C_pei1_NV > C_pei1_avg/weak_depth &
                            C_pei2_NV > C_pei2_avg/weak_depth &
                            C_pei3_NV > C_pei3_avg/weak_depth &
                            C_usa1_NV > C_usa1_avg/weak_depth &
                            C_usa2_NV > C_usa2_avg/weak_depth &
                            C_usa3_NV > C_usa3_avg/weak_depth &
                            C_usa4_NV > C_usa4_avg/weak_depth &
                            C_usa5_NV > C_usa5_avg/weak_depth)
    # Somtic mutations OR mismaps, maybe some founder germline that dropped out of some samples
        subsetCnoH <- filter(noHanyC, !(C_pei1_NV > C_pei1_avg/weak_depth &
                            C_pei2_NV > C_pei2_avg/weak_depth &
                            C_pei3_NV > C_pei3_avg/weak_depth &
                            C_usa1_NV > C_usa1_avg/weak_depth &
                            C_usa2_NV > C_usa2_avg/weak_depth &
                            C_usa3_NV > C_usa3_avg/weak_depth &
                            C_usa4_NV > C_usa4_avg/weak_depth &
                            C_usa5_NV > C_usa5_avg/weak_depth))
    # Sublineage specific
        anyUSAnoPEInoH <- filter(subsetCnoH, !(C_pei1_NV > C_pei1_avg/weak_depth |
                            C_pei2_NV > C_pei2_avg/weak_depth |
                            C_pei3_NV > C_pei3_avg/weak_depth) &
                            (C_usa1_NV > C_usa1_avg/weak_depth |
                            C_usa2_NV > C_usa2_avg/weak_depth |
                            C_usa3_NV > C_usa3_avg/weak_depth |
                            C_usa4_NV > C_usa4_avg/weak_depth |
                            C_usa5_NV > C_usa5_avg/weak_depth))
        anyPEInoUSAnoH <- filter(subsetCnoH, (C_pei1_NV > C_pei1_avg/weak_depth |
                            C_pei2_NV > C_pei2_avg/weak_depth |
                            C_pei3_NV > C_pei3_avg/weak_depth) &
                            !(C_usa1_NV > C_usa1_avg/weak_depth |
                            C_usa2_NV > C_usa2_avg/weak_depth |
                            C_usa3_NV > C_usa3_avg/weak_depth |
                            C_usa4_NV > C_usa4_avg/weak_depth |
                            C_usa5_NV > C_usa5_avg/weak_depth))        
        allUSAnoPEInoH <- filter(subsetCnoH, !(C_pei1_NV > C_pei1_avg/weak_depth |
                            C_pei2_NV > C_pei2_avg/weak_depth |
                            C_pei3_NV > C_pei3_avg/weak_depth) &
                            (C_usa1_NV > C_usa1_avg/weak_depth &
                            C_usa2_NV > C_usa2_avg/weak_depth &
                            C_usa3_NV > C_usa3_avg/weak_depth &
                            C_usa4_NV > C_usa4_avg/weak_depth &
                            C_usa5_NV > C_usa5_avg/weak_depth))
        allPEInoUSAnoH <- filter(subsetCnoH, (C_pei1_NV > C_pei1_avg/weak_depth &
                            C_pei2_NV > C_pei2_avg/weak_depth &
                            C_pei3_NV > C_pei3_avg/weak_depth) &
                            !(C_usa1_NV > C_usa1_avg/weak_depth |
                            C_usa2_NV > C_usa2_avg/weak_depth |
                            C_usa3_NV > C_usa3_avg/weak_depth |
                            C_usa4_NV > C_usa4_avg/weak_depth |
                            C_usa5_NV > C_usa5_avg/weak_depth))

nrow(anyC)
nrow(noC)
nrow(anyHanyC)
nrow(noHanyC)
nrow(anyHnoC)
nrow(noHnoC)
nrow(anyH)
nrow(allCanyH)
nrow(allCnoH)
nrow(anyUSAnoPEInoH)
nrow(anyPEInoUSAnoH)
nrow(allUSAnoPEInoH)
nrow(allPEInoUSAnoH)


# for(sample in c("C_pei1","C_pei2","C_pei3","C_usa1","C_usa2","C_usa3","C_usa4","C_usa5","H_ref","H_usa","H_pei")){
#     # print(sample)
#     count <- filter(snvs, get(paste0(sample, "_NV")) > get(paste0(sample, "_avg"))/stringent_depth) %>% nrow()
#     print(paste(sample, count))
# }

# for(sample in c("C_pei1","C_pei2","C_pei3","C_usa1","C_usa2","C_usa3","C_usa4","C_usa5","H_ref","H_usa","H_pei")){
#     count <- filter(snvs, get(paste0(sample, "_NV")) > get(paste0(sample, "_avg"))/weak_depth)) %>% nrow()
#     print(paste(sample, count))
# }




# for(sample in c("C_pei1","C_pei2","C_pei3","C_usa1","C_usa2","C_usa3","C_usa4","C_usa5","H_ref","H_usa","H_pei")){
#     # print(sample)
#     count <- filter(snvs, get(paste0(sample, "_NV")) > get(paste0(sample, "_avg"))/stringent_depth) %>% nrow()
#     print(paste(sample, count))
# }

# for(sample in c("C_pei1","C_pei2","C_pei3","C_usa1","C_usa2","C_usa3","C_usa4","C_usa5","H_ref","H_usa","H_pei")){
#     count <- filter(snvs, get(paste0(sample, "_NV")) > get(paste0(sample, "_avg"))/weak_depth)) %>% nrow()
#     print(paste(sample, count))
# }

