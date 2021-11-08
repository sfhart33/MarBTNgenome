# original file here: C:\Users\shart\Metzger Lab Dropbox\Sam_data\Mya_genome\SNPs\somatypus_output_sigS_upstream.r

library(tidyverse)

setwd("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run")

# remove last helmsman file   
    file.remove("helmsman.txt")
    #file.remove("dndscv.txt")

# run full binning loop on all
for(region in c("genome","gene","exon","CDS","five_prime_UTR","three_prime_UTR")){
#for(region in c("genome")){ If you just want to run for genome
    print(paste("Starting",region))
# Load full genome data
    if(region=="genome"){
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
    }
# Load gene regions data
    if(region!="genome"){
        snvs.nr <- read.table(paste0("Somatypus_SNVs_final.",region,"_NR.txt"), header=T, check.names=F)
        snvs.nv <- read.table(paste0("Somatypus_SNVs_final.",region,"_NV.txt"), header=T, check.names=F)
        snvs.metadata <- read.table(paste0("Somatypus_SNVs_final.",region,"_Metadata.txt"), header=T)
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
    }
    print("     data loaded")
#Set thresholds
    min_cutoff=3
    stringent_depth=6 # was 8
    weak_depth=16
#BINNING
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
        anyH <- rbind(anyHnoC, anyHanyC) %>% nrow()
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
    # Found in both sublineages, but not all samples: could be missed allC, mismaps, or LOH for just a few samples
        subsetCPandUnoH <- filter(subsetCnoH, (C_pei1_NV > C_pei1_avg/weak_depth |
                            C_pei2_NV > C_pei2_avg/weak_depth |
                            C_pei3_NV > C_pei3_avg/weak_depth) &
                            (C_usa1_NV > C_usa1_avg/weak_depth |
                            C_usa2_NV > C_usa2_avg/weak_depth |
                            C_usa3_NV > C_usa3_avg/weak_depth |
                            C_usa4_NV > C_usa4_avg/weak_depth |
                            C_usa5_NV > C_usa5_avg/weak_depth))
    # USA sublineage
        allUSAnoPEInoH <- filter(anyUSAnoPEInoH, 
                            (C_usa1_NV > C_usa1_avg/weak_depth &
                            C_usa2_NV > C_usa2_avg/weak_depth &
                            C_usa3_NV > C_usa3_avg/weak_depth &
                            C_usa4_NV > C_usa4_avg/weak_depth &
                            C_usa5_NV > C_usa5_avg/weak_depth))
        subsetUSAnoPEInoH <- filter(anyUSAnoPEInoH, 
                            !(C_usa1_NV > C_usa1_avg/weak_depth &
                            C_usa2_NV > C_usa2_avg/weak_depth &
                            C_usa3_NV > C_usa3_avg/weak_depth &
                            C_usa4_NV > C_usa4_avg/weak_depth &
                            C_usa5_NV > C_usa5_avg/weak_depth))
        multipleUSAnoPEInoH <- filter(subsetUSAnoPEInoH, 
                            !((C_usa1_NV > C_usa1_avg/weak_depth & !(C_usa2_NV > C_usa2_avg/weak_depth | C_usa3_NV > C_usa3_avg/weak_depth | C_usa4_NV > C_usa4_avg/weak_depth | C_usa5_NV > C_usa5_avg/weak_depth)) |
                            (C_usa2_NV > C_usa2_avg/weak_depth & !(C_usa1_NV > C_usa1_avg/weak_depth | C_usa3_NV > C_usa3_avg/weak_depth | C_usa4_NV > C_usa4_avg/weak_depth | C_usa5_NV > C_usa5_avg/weak_depth)) |
                            (C_usa3_NV > C_usa3_avg/weak_depth & !(C_usa2_NV > C_usa2_avg/weak_depth | C_usa1_NV > C_usa1_avg/weak_depth | C_usa4_NV > C_usa4_avg/weak_depth | C_usa5_NV > C_usa5_avg/weak_depth)) |
                            (C_usa4_NV > C_usa4_avg/weak_depth & !(C_usa2_NV > C_usa2_avg/weak_depth | C_usa3_NV > C_usa3_avg/weak_depth | C_usa1_NV > C_usa1_avg/weak_depth | C_usa5_NV > C_usa5_avg/weak_depth)) |
                            (C_usa5_NV > C_usa5_avg/weak_depth & !(C_usa2_NV > C_usa2_avg/weak_depth | C_usa3_NV > C_usa3_avg/weak_depth | C_usa4_NV > C_usa4_avg/weak_depth | C_usa1_NV > C_usa1_avg/weak_depth))))
        uniqUSAnoPEInoH <- filter(subsetUSAnoPEInoH, 
                            (C_usa1_NV > C_usa1_avg/weak_depth & !(C_usa2_NV > C_usa2_avg/weak_depth | C_usa3_NV > C_usa3_avg/weak_depth | C_usa4_NV > C_usa4_avg/weak_depth | C_usa5_NV > C_usa5_avg/weak_depth)) |
                            (C_usa2_NV > C_usa2_avg/weak_depth & !(C_usa1_NV > C_usa1_avg/weak_depth | C_usa3_NV > C_usa3_avg/weak_depth | C_usa4_NV > C_usa4_avg/weak_depth | C_usa5_NV > C_usa5_avg/weak_depth)) |
                            (C_usa3_NV > C_usa3_avg/weak_depth & !(C_usa2_NV > C_usa2_avg/weak_depth | C_usa1_NV > C_usa1_avg/weak_depth | C_usa4_NV > C_usa4_avg/weak_depth | C_usa5_NV > C_usa5_avg/weak_depth)) |
                            (C_usa4_NV > C_usa4_avg/weak_depth & !(C_usa2_NV > C_usa2_avg/weak_depth | C_usa3_NV > C_usa3_avg/weak_depth | C_usa1_NV > C_usa1_avg/weak_depth | C_usa5_NV > C_usa5_avg/weak_depth)) |
                            (C_usa5_NV > C_usa5_avg/weak_depth & !(C_usa2_NV > C_usa2_avg/weak_depth | C_usa3_NV > C_usa3_avg/weak_depth | C_usa4_NV > C_usa4_avg/weak_depth | C_usa1_NV > C_usa1_avg/weak_depth)))
    # Unique to one sample
        uniqU1noPEInoH <- filter(subsetUSAnoPEInoH, 
                            (C_usa1_NV > C_usa1_avg/weak_depth & !(C_usa2_NV > C_usa2_avg/weak_depth | C_usa3_NV > C_usa3_avg/weak_depth | C_usa4_NV > C_usa4_avg/weak_depth | C_usa5_NV > C_usa5_avg/weak_depth)))
        uniqU2noPEInoH <- filter(subsetUSAnoPEInoH, 
                            (C_usa2_NV > C_usa2_avg/weak_depth & !(C_usa1_NV > C_usa1_avg/weak_depth | C_usa3_NV > C_usa3_avg/weak_depth | C_usa4_NV > C_usa4_avg/weak_depth | C_usa5_NV > C_usa5_avg/weak_depth)))
        uniqU3noPEInoH <- filter(subsetUSAnoPEInoH, 
                            (C_usa3_NV > C_usa3_avg/weak_depth & !(C_usa2_NV > C_usa2_avg/weak_depth | C_usa1_NV > C_usa1_avg/weak_depth | C_usa4_NV > C_usa4_avg/weak_depth | C_usa5_NV > C_usa5_avg/weak_depth)))
        uniqU4noPEInoH <- filter(subsetUSAnoPEInoH, 
                            (C_usa4_NV > C_usa4_avg/weak_depth & !(C_usa2_NV > C_usa2_avg/weak_depth | C_usa3_NV > C_usa3_avg/weak_depth | C_usa1_NV > C_usa1_avg/weak_depth | C_usa5_NV > C_usa5_avg/weak_depth)))
        uniqU5noPEInoH <- filter(subsetUSAnoPEInoH, 
                            (C_usa5_NV > C_usa5_avg/weak_depth & !(C_usa2_NV > C_usa2_avg/weak_depth | C_usa3_NV > C_usa3_avg/weak_depth | C_usa4_NV > C_usa4_avg/weak_depth | C_usa1_NV > C_usa1_avg/weak_depth)))
    # look at all maine, etc
        multipleU1234noPEInoH <- filter(multipleUSAnoPEInoH, 
                                (C_usa1_NV > C_usa1_avg/stringent_depth) & 
                                (C_usa2_NV > C_usa2_avg/stringent_depth) & 
                                (C_usa3_NV > C_usa3_avg/stringent_depth) & 
                                (C_usa4_NV > C_usa4_avg/stringent_depth) & 
                                !(C_usa5_NV > C_usa5_avg/stringent_depth))
        multipleU1235noPEInoH <- filter(multipleUSAnoPEInoH, 
                                (C_usa1_NV > C_usa1_avg/stringent_depth) & 
                                (C_usa2_NV > C_usa2_avg/stringent_depth) & 
                                (C_usa3_NV > C_usa3_avg/stringent_depth) & 
                                !(C_usa4_NV > C_usa4_avg/stringent_depth) & 
                                (C_usa5_NV > C_usa5_avg/stringent_depth))
        multipleU1245noPEInoH <- filter(multipleUSAnoPEInoH, 
                                (C_usa1_NV > C_usa1_avg/stringent_depth) & 
                                (C_usa2_NV > C_usa2_avg/stringent_depth) & 
                                !(C_usa3_NV > C_usa3_avg/stringent_depth) & 
                                (C_usa4_NV > C_usa4_avg/stringent_depth) & 
                                (C_usa5_NV > C_usa5_avg/stringent_depth))
        multipleU1345noPEInoH <- filter(multipleUSAnoPEInoH, 
                                (C_usa1_NV > C_usa1_avg/stringent_depth) & 
                                !(C_usa2_NV > C_usa2_avg/stringent_depth) & 
                                (C_usa3_NV > C_usa3_avg/stringent_depth) & 
                                (C_usa4_NV > C_usa4_avg/stringent_depth) & 
                                (C_usa5_NV > C_usa5_avg/stringent_depth))
        multipleU2345noPEInoH <- filter(multipleUSAnoPEInoH, 
                                !(C_usa1_NV > C_usa1_avg/stringent_depth) & 
                                (C_usa2_NV > C_usa2_avg/stringent_depth) & 
                                (C_usa3_NV > C_usa3_avg/stringent_depth) & 
                                (C_usa4_NV > C_usa4_avg/stringent_depth) & 
                                (C_usa5_NV > C_usa5_avg/stringent_depth))          

        multipleU123noPEInoH <- filter(multipleUSAnoPEInoH, 
                                (C_usa1_NV > C_usa1_avg/stringent_depth) & 
                                (C_usa2_NV > C_usa2_avg/stringent_depth) & 
                                (C_usa3_NV > C_usa3_avg/stringent_depth) & 
                                !(C_usa4_NV > C_usa4_avg/stringent_depth) & 
                                !(C_usa5_NV > C_usa5_avg/stringent_depth))
        multipleU124noPEInoH <- filter(multipleUSAnoPEInoH, 
                                (C_usa1_NV > C_usa1_avg/stringent_depth) & 
                                (C_usa2_NV > C_usa2_avg/stringent_depth) & 
                                !(C_usa3_NV > C_usa3_avg/stringent_depth) & 
                                (C_usa4_NV > C_usa4_avg/stringent_depth) & 
                                !(C_usa5_NV > C_usa5_avg/stringent_depth))
        multipleU134noPEInoH <- filter(multipleUSAnoPEInoH, 
                                (C_usa1_NV > C_usa1_avg/stringent_depth) & 
                                !(C_usa2_NV > C_usa2_avg/stringent_depth) & 
                                (C_usa3_NV > C_usa3_avg/stringent_depth) & 
                                (C_usa4_NV > C_usa4_avg/stringent_depth) & 
                                !(C_usa5_NV > C_usa5_avg/stringent_depth))
        multipleU234noPEInoH <- filter(multipleUSAnoPEInoH, 
                                !(C_usa1_NV > C_usa1_avg/stringent_depth) & 
                                (C_usa2_NV > C_usa2_avg/stringent_depth) & 
                                (C_usa3_NV > C_usa3_avg/stringent_depth) & 
                                (C_usa4_NV > C_usa4_avg/stringent_depth) & 
                                !(C_usa5_NV > C_usa5_avg/stringent_depth))

        multipleU12noPEInoH <- filter(multipleUSAnoPEInoH, 
                                (C_usa1_NV > C_usa1_avg/stringent_depth) & 
                                (C_usa2_NV > C_usa2_avg/stringent_depth) & 
                                !(C_usa3_NV > C_usa3_avg/stringent_depth) & 
                                !(C_usa4_NV > C_usa4_avg/stringent_depth) & 
                                !(C_usa5_NV > C_usa5_avg/stringent_depth))
        multipleU13noPEInoH <- filter(multipleUSAnoPEInoH, 
                                (C_usa1_NV > C_usa1_avg/stringent_depth) & 
                                !(C_usa2_NV > C_usa2_avg/stringent_depth) & 
                                (C_usa3_NV > C_usa3_avg/stringent_depth) & 
                                !(C_usa4_NV > C_usa4_avg/stringent_depth) & 
                                !(C_usa5_NV > C_usa5_avg/stringent_depth))
        multipleU14noPEInoH <- filter(multipleUSAnoPEInoH, 
                                (C_usa1_NV > C_usa1_avg/stringent_depth) & 
                                !(C_usa2_NV > C_usa2_avg/stringent_depth) & 
                                !(C_usa3_NV > C_usa3_avg/stringent_depth) & 
                                (C_usa4_NV > C_usa4_avg/stringent_depth) & 
                                !(C_usa5_NV > C_usa5_avg/stringent_depth))
        multipleU23noPEInoH <- filter(multipleUSAnoPEInoH, 
                                !(C_usa1_NV > C_usa1_avg/stringent_depth) & 
                                (C_usa2_NV > C_usa2_avg/stringent_depth) & 
                                (C_usa3_NV > C_usa3_avg/stringent_depth) & 
                                !(C_usa4_NV > C_usa4_avg/stringent_depth) & 
                                !(C_usa5_NV > C_usa5_avg/stringent_depth))
        multipleU24noPEInoH <- filter(multipleUSAnoPEInoH, 
                                !(C_usa1_NV > C_usa1_avg/stringent_depth) & 
                                (C_usa2_NV > C_usa2_avg/stringent_depth) & 
                                !(C_usa3_NV > C_usa3_avg/stringent_depth) & 
                                (C_usa4_NV > C_usa4_avg/stringent_depth) & 
                                !(C_usa5_NV > C_usa5_avg/stringent_depth))
        multipleU34noPEInoH <- filter(multipleUSAnoPEInoH, 
                                !(C_usa1_NV > C_usa1_avg/stringent_depth) & 
                                !(C_usa2_NV > C_usa2_avg/stringent_depth) & 
                                (C_usa3_NV > C_usa3_avg/stringent_depth) & 
                                (C_usa4_NV > C_usa4_avg/stringent_depth) & 
                                !(C_usa5_NV > C_usa5_avg/stringent_depth))

        multipleU1noPEInoH <- filter(multipleUSAnoPEInoH, C_usa1_NV > C_usa1_avg/stringent_depth)
        multipleU2noPEInoH <- filter(multipleUSAnoPEInoH, C_usa2_NV > C_usa2_avg/stringent_depth)
        multipleU3noPEInoH <- filter(multipleUSAnoPEInoH, C_usa3_NV > C_usa3_avg/stringent_depth)
        multipleU4noPEInoH <- filter(multipleUSAnoPEInoH, C_usa4_NV > C_usa4_avg/stringent_depth)
        multipleU5noPEInoH <- filter(multipleUSAnoPEInoH, C_usa5_NV > C_usa5_avg/stringent_depth)
    # PEI sublineage
        allPEInoUSAnoH <- filter(anyPEInoUSAnoH, 
                            (C_pei1_NV > C_pei1_avg/weak_depth &
                            C_pei2_NV > C_pei2_avg/weak_depth &
                            C_pei3_NV > C_pei3_avg/weak_depth))
        subsetPEInoUSAnoH <- filter(anyPEInoUSAnoH, 
                            !(C_pei1_NV > C_pei1_avg/weak_depth &
                            C_pei2_NV > C_pei2_avg/weak_depth &
                            C_pei3_NV > C_pei3_avg/weak_depth))                        
        multiplePEInoUSAnoH <- filter(subsetPEInoUSAnoH, 
                            !((C_pei1_NV > C_pei1_avg/weak_depth & !(C_pei2_NV > C_pei2_avg/weak_depth | C_pei3_NV > C_pei3_avg/weak_depth)) |
                            (C_pei2_NV > C_pei2_avg/weak_depth & !(C_pei1_NV > C_pei1_avg/weak_depth | C_pei3_NV > C_pei3_avg/weak_depth)) |
                            (C_pei3_NV > C_pei3_avg/weak_depth & !(C_pei2_NV > C_pei2_avg/weak_depth | C_pei1_NV > C_pei1_avg/weak_depth))))  
        uniqPEInoUSAnoH <- filter(subsetPEInoUSAnoH, 
                            (C_pei1_NV > C_pei1_avg/weak_depth & !(C_pei2_NV > C_pei2_avg/weak_depth | C_pei3_NV > C_pei3_avg/weak_depth)) |
                            (C_pei2_NV > C_pei2_avg/weak_depth & !(C_pei1_NV > C_pei1_avg/weak_depth | C_pei3_NV > C_pei3_avg/weak_depth)) |
                            (C_pei3_NV > C_pei3_avg/weak_depth & !(C_pei2_NV > C_pei2_avg/weak_depth | C_pei1_NV > C_pei1_avg/weak_depth)))  
        uniqP1noUSAnoH <- filter(subsetPEInoUSAnoH, 
                            (C_pei1_NV > C_pei1_avg/weak_depth & !(C_pei2_NV > C_pei2_avg/weak_depth | C_pei3_NV > C_pei3_avg/weak_depth))) 
        uniqP2noUSAnoH <- filter(subsetPEInoUSAnoH, 
                            (C_pei2_NV > C_pei2_avg/weak_depth & !(C_pei1_NV > C_pei1_avg/weak_depth | C_pei3_NV > C_pei3_avg/weak_depth))) 
        uniqP3noUSAnoH <- filter(subsetPEInoUSAnoH, 
                            (C_pei3_NV > C_pei3_avg/weak_depth & !(C_pei2_NV > C_pei2_avg/weak_depth | C_pei1_NV > C_pei1_avg/weak_depth))) 
        multipleP12noUSAnoH <- filter(multiplePEInoUSAnoH, 
                            (C_pei1_NV > C_pei1_avg/stringent_depth) & (C_pei2_NV > C_pei2_avg/stringent_depth) & !(C_pei3_NV > C_pei3_avg/stringent_depth))
        multipleP13noUSAnoH <- filter(multiplePEInoUSAnoH, 
                            (C_pei1_NV > C_pei1_avg/stringent_depth) & !(C_pei2_NV > C_pei2_avg/stringent_depth) & (C_pei3_NV > C_pei3_avg/stringent_depth))
        multipleP23noUSAnoH <- filter(multiplePEInoUSAnoH, 
                            !(C_pei1_NV > C_pei1_avg/stringent_depth) & (C_pei2_NV > C_pei2_avg/stringent_depth) & (C_pei3_NV > C_pei3_avg/stringent_depth))                        
        multipleP1noUSAnoH <- filter(multiplePEInoUSAnoH, C_pei1_NV > C_pei1_avg/stringent_depth)
        multipleP2noUSAnoH <- filter(multiplePEInoUSAnoH, C_pei2_NV > C_pei2_avg/stringent_depth)
        multipleP3noUSAnoH <- filter(multiplePEInoUSAnoH, C_pei3_NV > C_pei3_avg/stringent_depth)

    #bins for helmsman run
    keeps <- c("allCanyH", "allCnoH", "allPEInoUSAnoH", "allUSAnoPEInoH",
    "H_P", "H_R", "H_RP", "H_RU", "H_RUP", "H_U", "H_UP", 
    "multiplePEInoUSAnoH", "multipleUSAnoPEInoH",
    "multipleP12noUSAnoH", "multipleP13noUSAnoH", "multipleP23noUSAnoH",
    "multipleU1234noPEInoH", "multipleU1235noPEInoH", "multipleU123noPEInoH", "multipleU1245noPEInoH", "multipleU124noPEInoH", "multipleU12noPEInoH", "multipleU1345noPEInoH", "multipleU134noPEInoH", "multipleU13noPEInoH", "multipleU14noPEInoH", "multipleU2345noPEInoH", "multipleU234noPEInoH", "multipleU23noPEInoH", "multipleU24noPEInoH", "multipleU34noPEInoH",
    "subsetCanyH", "subsetCPandUnoH",
    "uniqP1noUSAnoH", "uniqP2noUSAnoH", "uniqP3noUSAnoH", "uniqU1noPEInoH", "uniqU2noPEInoH", "uniqU3noPEInoH", "uniqU4noPEInoH", "uniqU5noPEInoH")

    helmsman_build <- function(input, name){
        #name <- deparse(substitute(input))
        helmsman <- input %>% select(CHROM,POS,REF,ALT) %>%
            mutate(id = paste(name,region, sep="."))
        #dndscv <- helmsman %>% select(id, CHROM,POS,REF,ALT)
        write.table(helmsman, file = paste0("helmsman.txt"), append = TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep = '\t')
        #write.table(dndscv,file = paste0("dndscv.txt"), append = TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep = '\t')
    }
    for(i in keeps){
        helmsman_build(get(i), i)
    }
}
