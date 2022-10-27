# Original file here: C:\Users\shart\Metzger Lab Dropbox\Sam_data\Mya_genome\revision\purity_nuclear.R
library(tidyverse)
library(gridExtra)

# Load full genome data and calculate thresholds
    setwd("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run")

    # Load data to calculate thresholds
        snvs.nr <- read.table("Somatypus_SNVs_final_NR.txt", header=T, check.names=F)
        snvs.nv <- read.table("Somatypus_SNVs_final_NV.txt", header=T, check.names=F)
        snvs.metadata <- read.table("Somatypus_SNVs_final_Metadata.txt", header=T)

    # Define sample names and combine into one table
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

    #Set thresholds
        min_cutoff=3
        stringent_depth=6 # was 8
        weak_depth=16

# Load LOH regions data
    setwd("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/validation")
    for(region in c("bothCN2")){ # , "bothCN4"
            print(paste("Starting:", region, "snvs"))
            # Load dataSNVs_bothCN2
                snvs.nr <- read.table(paste("SNVs", region, "NR.txt", sep = "_"), header=T, check.names=F)
                snvs.nv <- read.table(paste("SNVs", region, "NV.txt", sep = "_"), header=T, check.names=F)
                snvs.metadata <- read.table(paste("SNVs", region, "Metadata.txt", sep = "_"), header=T)
            # Define sample names and combine into one table
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
            # Save table for each loop
                assign(paste(region, "snvs", sep = "_"), cbind(snvs.metadata2,snvs.nv,snvs.nr))
    }

# remove extra cols, calculate VAFs
for(region in c("bothCN2")){ # , "bothCN4"
    snvs2 <- get(paste(region, "snvs", sep = "_")) %>%
        # select(-c(C_pei4_NV, C_usa6_NV, C_usa7_NV, C_pei4T_NV, C_usa6T_NV, C_usa7T_NV,
        #           C_pei4_NR, C_usa6_NR, C_usa7_NR, C_pei4T_NR, C_usa6T_NR, C_usa7T_NR)) %>%
        mutate(pei_VAF = (C_pei1_NV + C_pei2_NV + C_pei3_NV) / (C_pei1_NR + C_pei2_NR + C_pei3_NR),
               usa_VAF = (C_usa1_NV + C_usa2_NV + C_usa3_NV + C_usa4_NV + C_usa5_NV) / (C_usa1_NR + C_usa2_NR + C_usa3_NR + C_usa4_NR + C_usa5_NR),
               Href_VAF = H_ref_NV / H_ref_NR,
               Husa_VAF = H_usa_NV / H_usa_NR,
               Hpei_VAF = H_pei_NV / H_pei_NR,
               pei1_VAF = C_pei1_NV / C_pei1_NR,
               pei2_VAF = C_pei2_NV / C_pei2_NR,
               pei3_VAF = C_pei3_NV / C_pei3_NR,
               pei4_VAF = C_pei4_NV / C_pei4_NR,
               usa1_VAF = C_usa1_NV / C_usa1_NR,
               usa2_VAF = C_usa2_NV / C_usa2_NR,
               usa3_VAF = C_usa3_NV / C_usa3_NR,
               usa4_VAF = C_usa4_NV / C_usa4_NR,
               usa5_VAF = C_usa5_NV / C_usa5_NR,
               usa6_VAF = C_usa6_NV / C_usa6_NR,
               usa7_VAF = C_usa7_NV / C_usa7_NR,
               pei2T_VAF = C_pei2T_NV / C_pei2T_NR,
               pei3T_VAF = C_pei3T_NV / C_pei3T_NR,
               pei4T_VAF = C_pei4T_NV / C_pei4T_NR,
               usa3T_VAF = C_usa3T_NV / C_usa3T_NR,
               usa4T_VAF = C_usa4T_NV / C_usa4T_NR,
               usa5T_VAF = C_usa5T_NV / C_usa5T_NR,
               usa6T_VAF = C_usa6T_NV / C_usa6T_NR,
               usa7T_VAF = C_usa7T_NV / C_usa7T_NR 
               )
        assign(paste(region, "snvs2", sep = "_"), snvs2)
}
head(bothCN2_snvs2)

#BINNING - same as for signature extraction
    # Is SNV found in any cancer
        anyC <- filter(bothCN2_snvs2, C_pei1_NV > C_pei1_avg/stringent_depth |
                            C_pei2_NV > C_pei2_avg/stringent_depth |
                            C_pei3_NV > C_pei3_avg/stringent_depth |
                            C_usa1_NV > C_usa1_avg/stringent_depth |
                            C_usa2_NV > C_usa2_avg/stringent_depth |
                            C_usa3_NV > C_usa3_avg/stringent_depth |
                            C_usa4_NV > C_usa4_avg/stringent_depth |
                            C_usa5_NV > C_usa5_avg/stringent_depth)
    # Is SNV found in any healthy clam    
        anyHanyC <- filter(anyC, H_ref_NV > min_cutoff |
                            H_usa_NV > min_cutoff |
                            H_pei_NV > min_cutoff)
        noHanyC <- filter(anyC, !(H_ref_NV > min_cutoff |
                            H_usa_NV > min_cutoff |
                            H_pei_NV > min_cutoff))
        anyH <- filter(bothCN2_snvs2, H_ref_NV > min_cutoff |
                            H_usa_NV > min_cutoff |
                            H_pei_NV > min_cutoff)                                       
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
                            C_usa5_NV > C_usa5_avg/weak_depth)) %>% mutate(bin = "allUSAnoPEInoH")
        allPEInoUSAnoH <- filter(subsetCnoH, (C_pei1_NV > C_pei1_avg/weak_depth &
                            C_pei2_NV > C_pei2_avg/weak_depth &
                            C_pei3_NV > C_pei3_avg/weak_depth) &
                            !(C_usa1_NV > C_usa1_avg/weak_depth |
                            C_usa2_NV > C_usa2_avg/weak_depth |
                            C_usa3_NV > C_usa3_avg/weak_depth |
                            C_usa4_NV > C_usa4_avg/weak_depth |
                            C_usa5_NV > C_usa5_avg/weak_depth)) %>% mutate(bin = "allPEInoUSAnoH")

        allUSAnoPEIanyH <- filter(subsetCanyH, !(C_pei1_NV > C_pei1_avg/weak_depth |
                            C_pei2_NV > C_pei2_avg/weak_depth |
                            C_pei3_NV > C_pei3_avg/weak_depth) &
                            (C_usa1_NV > C_usa1_avg/weak_depth &
                            C_usa2_NV > C_usa2_avg/weak_depth &
                            C_usa3_NV > C_usa3_avg/weak_depth &
                            C_usa4_NV > C_usa4_avg/weak_depth &
                            C_usa5_NV > C_usa5_avg/weak_depth)) %>% mutate(bin = "allUSAnoPEIanyH")
        allPEInoUSAanyH <- filter(subsetCanyH, (C_pei1_NV > C_pei1_avg/weak_depth &
                            C_pei2_NV > C_pei2_avg/weak_depth &
                            C_pei3_NV > C_pei3_avg/weak_depth) &
                            !(C_usa1_NV > C_usa1_avg/weak_depth |
                            C_usa2_NV > C_usa2_avg/weak_depth |
                            C_usa3_NV > C_usa3_avg/weak_depth |
                            C_usa4_NV > C_usa4_avg/weak_depth |
                            C_usa5_NV > C_usa5_avg/weak_depth)) %>% mutate(bin = "allPEInoUSAanyH")

homoz_cutoff <- 0.8
homozT_cutoff <- 0.9
allCnoH_homoz <- filter(allCnoH,
                            pei1_VAF > homoz_cutoff,
                            pei2_VAF > homoz_cutoff,
                            pei3_VAF > homoz_cutoff,
                            usa1_VAF > homoz_cutoff,
                            usa2_VAF > homoz_cutoff,
                            usa3_VAF > homoz_cutoff,
                            usa4_VAF > homoz_cutoff,
                            usa5_VAF > homoz_cutoff)
allCnoH_homoz2 <- filter(allCnoH,
                            pei1_VAF > homoz_cutoff,
                            pei2_VAF > homoz_cutoff,
                            pei3_VAF > homoz_cutoff,
                            usa1_VAF > homoz_cutoff,
                            usa2_VAF > homoz_cutoff,
                            usa3_VAF > homoz_cutoff,
                            usa4_VAF > homoz_cutoff,
                            usa5_VAF > homoz_cutoff,
                            pei2T_VAF < homozT_cutoff,
                            pei3T_VAF < homozT_cutoff,
                            usa3T_VAF < homozT_cutoff,
                            usa4T_VAF < homozT_cutoff,
                            usa5T_VAF < homozT_cutoff)

samplesALL <- c("Href", " Husa", "Hpei", "usa1","usa2","usa3","usa4","usa5","usa6","usa7","pei1","pei2","pei3","pei4")
samplesC <- c("usa1","usa2","usa3","usa4","usa5","pei1","pei2","pei3")
samplesU <- c("usa1","usa2","usa3","usa4","usa5")
samplesP <- c("pei1","pei2","pei3")
samplesH <- c("Husa", "Hpei")
samplesT <- c("usa3","usa4","usa5","pei2","pei3")
samplesTU <- c("usa3","usa4","usa5")

for(sample in samplesC){
    print(sample)
    median(allCnoH_homoz[,(paste0(sample,"_VAF"))]) %>%  print()
}
for(sample in samplesC){
    print(sample)
    mean(allCnoH_homoz[,(paste0(sample,"_VAF"))]) %>%  print()
    mean(allCnoH_homoz2[,(paste0(sample,"_VAF"))]) %>%  print()
}

for(sample in samplesC){
    print(sample)
    median(allCnoH_homoz[,(paste0(sample,"_VAF"))]) %>%  print()
    median(allCnoH_homoz2[,(paste0(sample,"_VAF"))]) %>%  print()
}

mean_corr <- (mean(filter(anyH, Husa_VAF > 0.8)[,"Husa_VAF"]) + mean(filter(anyH, Hpei_VAF > 0.8)[,"Hpei_VAF"])) / 2

for(sample in samples){
    print(sample)
    median(allCnoH_homoz[,paste0(sample, "_VAF")]) %>% print()
    mean(allCnoH_homoz[,paste0(sample, "_VAF")]) %>% print()
    sd(allCnoH_homoz[,paste0(sample, "_VAF")]) %>% print()
    (mean(allCnoH_homoz[,paste0(sample, "_VAF")]) + mean_corr) %>% print()
}

for(sample in samplesH){
    plot1 <- filter(anyH, get(paste0(sample, "_VAF")) > 0.8)
    print(sample)
    mean(plot1[,paste0(sample, "_VAF")]) %>% print()
}

    vafnames <- NULL
    vafs <- NULL
    for(sample in samplesU){
        vafs1 <- allCnoH_homoz[,paste0(sample, "_VAF")]
        vafnames1 <- rep(sample, nrow(allCnoH_homoz))
        vafs <- c(vafs, vafs1)
        vafnames <- c(vafnames, vafnames1)
    }
    for(sample in samplesP){
        vafs1 <- allCnoH_homoz[,paste0(sample, "_VAF")]
        vafnames1 <- rep(sample, nrow(allCnoH_homoz))
        vafs <- c(vafs, vafs1)
        vafnames <- c(vafnames, vafnames1)
    }
    for(sample in samplesH){
        anyHsubset <- filter(anyH, get(paste0(sample, "_VAF")) > 0.8)
        print(nrow(anyHsubset))
        vafs1 <- anyHsubset[,paste0(sample, "_VAF")]
        vafnames1 <- rep(sample, nrow(anyHsubset))
        print(nrow(anyHsubset))
        vafs <- c(vafs, vafs1)
        vafnames <- c(vafnames, vafnames1)
    }
    purity <- data.frame(vafs,vafnames) 
        head(purity, n=10)
    purity_labs  <-group_by(purity, vafnames) %>% summarize(median = median(vafs))

    vafnames <- NULL
    vafs <- NULL
    for(sample in samplesTU){
        allCnoH_homoz3 <- allCnoH_homoz2 %>%
            filter(get(paste0(sample, "T_VAF"))>0)
        vafs1 <- allCnoH_homoz3[,paste0(sample, "T_VAF")]
        vafnames1 <- rep(sample, nrow(allCnoH_homoz3))
        vafs <- c(vafs, vafs1)
        vafnames <- c(vafnames, vafnames1)
    }
    purityT <- data.frame(vafs,vafnames) 
    purityT_labs  <-group_by(purityT, vafnames) %>% summarize(median = median(vafs))

# pdf("purity_plot_tests.pdf")
# ggplot(purity, aes(x = vafnames, y = vafs))+
#         #geom_violin(fill = c(rep("black",2),rep("red",3),rep("blue",5)))+
#         geom_violin()+
#         xlab("samples")+
#         ylab("VAF for homozygous cancer-specific SNVs (CN2 only)")+
#         theme_classic()
# ggplot(purityT, aes(x = vafnames, y = vafs))+
#         #geom_violin(fill = c(rep("black",2),rep("red",3),rep("blue",5)))+
#         geom_violin()+
#         xlab("samples")+
#         ylab("VAF for homozygous cancer-specific SNVs in tissue (CN2 only)")+
#         theme_classic()
# ggplot(purity, aes(x = vafnames, y = vafs))+
#         geom_boxplot(outlier.shape = NA, fill = c(rep("black",2),rep("red",3),rep("blue",5)))+
#         geom_text(data = purity_labs, aes(x=vafnames, y=median, label=round(median, digits=2)), vjust = -0.75)+
#         xlab("samples")+
#         ylab("VAF for homozygous cancer-specific SNVs (CN2 only)")+
#         theme_classic()
# ggplot(purityT, aes(x = vafnames, y = vafs))+
#         geom_boxplot(outlier.shape = NA, fill = c(rep("red",2),rep("blue",3)))+
#         geom_text(data = purityT_labs, aes(x=vafnames, y=median, label=round(median, digits=2)), vjust = -0.75)+
#         xlab("samples")+
#         ylab("VAF for homozygous cancer-specific SNVs in tissue (CN2 only)")+
#         theme_classic()
# dev.off()

# pdf("CN2_VAF_hist_withH.pdf")
# for(sample in samples){
#     plot1 <- ggplot(allCnoH_homoz, aes(get(paste0(sample, "_VAF")))) +
#         geom_histogram(binwidth = 0.01)+
#         xlim(0.8,1.02)+
#         xlab(paste(sample," VAF (allCnoH"))+
#         theme_classic()
#     print(plot1)
# }
# for(sample in samplesH){
#     plot1 <- filter(anyH, get(paste0(sample, "_VAF")) > 0.8) %>%
#     ggplot(aes(get(paste0(sample, "_VAF")))) +
#         geom_histogram(binwidth = 0.01)+
#         xlim(0.8,1.02)+
#         xlab(paste(sample," VAF (anyH)"))+
#         theme_classic()
#     print(plot1)
# }
# dev.off()



######### compare with mito
# From: C:\Users\shart\Metzger Lab Dropbox\Sam_data\Mya_genome\mitochondria\mito_somatypus_analysis.r
    snvsmt <- read.delim("/ssd3/Mar_genome_analysis/bwa_mapping/mito/all_samples/somatypus/Somatypus_SNVs_final.counts", header = TRUE)

#not in dloop region, correct for name differences
    snvsmt <- filter(snvsmt, pos < 12060 | pos > 12971) %>%
        select(Href_f,Husa_f,Hpei_f,Cpei0_f,Cpei1_f,Cpei3_f,Cusa0a_f,Cusa0b_f,Cusa2_f,Cusa4_f,Cusa5_f,Tpei1_f,Tpei3_f,Tusa2_f,Tusa4_f,Tusa5_f)
    colnames(snvsmt) <- c("Href_f","Husa_f","Hpei_f","pei1_f","pei2_f","pei3_f","usa1_f","usa2_f","usa3_f","usa4_f","usa5_f","pei2T_f","pei3T_f","usa3T_f","usa4T_f","usa5T_f")
# Filter for snvs not found in healthies
    snvs_noH <- filter(snvsmt, Href_f < 0.5, Husa_f < 0.5, Hpei_f < 0.5)

    #######
    vafs <- NULL
    vafnames <- NULL
    for(sample in samplesC){
        subset <- filter(snvs_noH, get(paste0(sample,"_f")) > 0.5)
        vafs1 <- subset[,paste0(sample,"_f")]
        vafnames1 <- rep(sample, nrow(subset))
        vafs <- c(vafs, vafs1)
        vafnames <- c(vafnames, vafnames1)
    }
    vafs1 <- NULL
    vafnames1 <- NULL
    for(sample in c("Href", "Husa", "Hpei")){
        subset <- filter(snvsmt, get(paste0(sample,"_f")) > 0.5)
        vafs1 <- subset[,paste0(sample,"_f")]
        vafnames1 <- rep(sample, nrow(subset))
        vafs <- c(vafs, vafs1)
        vafnames <- c(vafnames, vafnames1)
    }
    purityM <- data.frame(vafs,vafnames) 
        head(purityM, n=10)
   purityM_labs  <-group_by(purityM, vafnames) %>% summarize(median = median(vafs))

    vafnames <- NULL
    vafs <- NULL
    for(sample in samplesTU){
        subset <- filter(snvs_noH, get(paste0(sample,"_f")) > 0.5)
        vafs1 <- subset[,paste0(sample,"T_f")]
        vafnames1 <- rep(sample, nrow(subset))
        vafs <- c(vafs, vafs1)
        vafnames <- c(vafnames, vafnames1)
    }
    purityMT <- data.frame(vafs,vafnames) 
    purityMT_labs <- group_by(purityT, vafnames) %>% summarize(median = median(vafs))
purityMT_labs 


p1 <- ggplot(purity, aes(x = vafnames, y = vafs))+
        geom_boxplot(outlier.shape = NA, fill = c(rep("black",2),rep("red",3),rep("blue",5)))+
        geom_text(data = purity_labs, aes(x=vafnames, y=1, label=round(median, digits=3)), vjust = -0.75, size=2)+
        ylim(0.79,1) +
        xlab(NA)+
        ylab("SAMPLE PURITY\n(VAF for homozygous CN2\ncancer-specific SNVs)")+
        theme_classic()+
        theme(axis.title.x=element_blank(),
            axis.text=element_text(size=6)
        )
p2 <-ggplot(purityT, aes(x = vafnames, y = vafs))+
        #geom_boxplot(outlier.shape = NA, fill = c(rep("red",2),rep("blue",3)))+
        geom_boxplot(outlier.shape = NA, fill = rep("blue",3))+
        geom_text(data = purityT_labs, aes(x=vafnames, y=median, label=round(median, digits=3)), vjust = -1.5, size=2)+
        ylim(0,1) +
        xlab(NA)+
        ylab("CANCER DISSEMINATION\n(VAF for homozygous CN2\ncancer-specific SNVs in paired tissue)")+
        theme_classic()+
        theme(axis.title.x=element_blank(),
            axis.text=element_text(size=6)
        )
p3 <-ggplot(purityM, aes(x = vafnames, y = vafs))+
        geom_boxplot(outlier.shape = NA, fill = c(rep("black",3),rep("red",3),rep("blue",5)))+
        geom_text(data = purityM_labs, aes(x=vafnames, y=1, label=round(median, digits=3)), vjust = -0.75, size=2)+
        ylim(0.79,1) +
        xlab(NA)+
        ylab("SAMPLE PURITY\n(VAF for mitochondrial\ncancer-specific SNVs)")+
        theme_classic()+
        theme(axis.title.x=element_blank(),
            axis.text=element_text(size=6)
        )
p4 <-ggplot(purityMT, aes(x = vafnames, y = vafs))+
        #geom_boxplot(outlier.shape = NA, fill = c(rep("red",2),rep("blue",3)))+
        geom_boxplot(outlier.shape = NA, fill = rep("blue",3))+
        geom_text(data = purityMT_labs, aes(x=vafnames, y=median, label=round(median, digits=3)), vjust = -1.5, size=2)+
        ylim(0,1) +
        xlab(NA)+
        ylab("CANCER DISSEMINATION\n(VAF for mitochondrial\ncancer-specificSNVs in paired tissue)")+
        theme_classic()+
        theme(axis.title.x=element_blank(),
            axis.text=element_text(size=6)
        )
pdf("purity_plots.pdf")
grid.arrange(p3,p4,p1,p2, nrow = 2, widths=c(2,1))
dev.off()

# pdf("purity_tissue_hists.pdf")
# for(sample in samplesT){
#     allCnoH_homoz3 <- allCnoH_homoz2 %>%
#         filter(get(paste0(sample, "T_VAF"))>0)
#     plot1 <- ggplot(allCnoH_homoz3, aes(get(paste0(sample, "T_VAF")))) +
#         geom_histogram(binwidth = 0.01)+
#         xlab(paste(sample," VAF (allCnoH"))+
#         theme_classic()
#     print(plot1)
# }
# for(sample in samplesT){
#     subset <- filter(snvs_noH, get(paste0(sample,"_f")) > 0.5)
#     plot1 <- ggplot(subset, aes(get(paste0(sample, "T_f")))) +
#         geom_histogram(binwidth = 0.01)+
#         xlab(paste(sample," VAF (mito"))+
#         theme_classic()
#     print(plot1)
# }
# dev.off()

# Estimating mtDNA copies per cell
    snvsmt2 <- read.delim("/ssd3/Mar_genome_analysis/bwa_mapping/mito/all_samples/somatypus/Somatypus_SNVs_final.counts", header = TRUE) %>%
        filter(pos < 12060 | pos > 12971) %>%
            select(Href_t,Husa_t,Hpei_t,Cpei0_t,Cpei1_t,Cpei3_t,Cusa0a_t,Cusa0b_t,Cusa2_t,Cusa4_t,Cusa5_t,Tpei1_t,Tpei3_t,Tusa2_t,Tusa4_t,Tusa5_t)
        colnames(snvsmt2) <- c("Href_t","Husa_t","Hpei_t","pei1_t","pei2_t","pei3_t","usa1_t","usa2_t","usa3_t","usa4_t","usa5_t","pei2T_t","pei3T_t","usa3T_t","usa4T_t","usa5T_t")

    colMeans(snvsmt2)
    #colMeans(snvsmt2cancer[,30:40])

    healthy_cov <- mean(colMeans(snvsmt2[,1:3]))
    colMeans(snvsmt2)/healthy_cov


            samples_cov <- c("H_ref","H_usa","H_pei",
            "C_pei1","C_pei2","C_pei3", # C_pei4 has high host contamination
            "C_usa1","C_usa2","C_usa3","C_usa4","C_usa5", # C_usa6 & 7 have high host contamination
            "C_pei2T","C_pei3T","C_usa3T","C_usa4T","C_usa5T")
            genome_cov <- c()
            for(sample in samples_cov){
                print(sample)
                print(get(paste0(sample,"_avg")))
                genome_cov <- append(genome_cov,get(paste0(sample,"_avg")))
            }

    genome_vs_mt <- data.frame(samples = colnames(snvsmt2), mt_cov = colMeans(snvsmt2), g_cov = genome_cov) %>%
        mutate(mt_vs_g = mt_cov / g_cov) %>%
        print()



