# Original file here: C:\Users\shart\Metzger Lab Dropbox\Sam_data\Mya_genome\LOH\LOH_revision_validation.R

library(tidyverse)
library(gridExtra)

# Load full genome data and calculate thresholds
    setwd("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run")

    # Load data if running for all SNVs
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
    for(region in c("USA", "PEI")){
        for(loh in c("LOH")){
            print(paste("Starting:", region, loh, "snvs", sep = "_"))
            # Load data
                snvs.nr <- read.table(paste("SNVs", region, loh, "NR.txt", sep = "_"), header=T, check.names=F)
                snvs.nv <- read.table(paste("SNVs", region, loh, "NV.txt", sep = "_"), header=T, check.names=F)
                snvs.metadata <- read.table(paste("SNVs", region, loh, "Metadata.txt", sep = "_"), header=T)
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
                assign(paste(region, loh, "snvs", sep = "_"), cbind(snvs.metadata2,snvs.nv,snvs.nr))
        }
    }

# Check number SNV positions in each set
    for(file in c("USA_LOH_snvs", "PEI_LOH_snvs", "snvs")){
         nrow(get(file))
    }

# merge into one file (adding column with the LOH and region status)
    USA_LOH_snvs <- USA_LOH_snvs %>%
        mutate(USA_LOH = "LOH")
    PEI_LOH_snvs <- PEI_LOH_snvs %>%
        mutate(PEI_LOH = "LOH")
        
    snvs_LOH <- left_join(snvs, USA_LOH_snvs) %>%
        left_join(PEI_LOH_snvs) %>%
        mutate(PEI_LOH = case_when(is.na(PEI_LOH) ~ "nonLOH", !is.na(PEI_LOH) ~ "LOH"),
               USA_LOH = case_when(is.na(USA_LOH) ~ "nonLOH", !is.na(USA_LOH) ~ "LOH"))
    
# remove extra cols, calculate VAF for cancer samples   
    snvs_LOH2 <- snvs_LOH %>%
        select(-c(C_pei4_NV, C_usa6_NV, C_usa7_NV, C_pei2T_NV, C_pei4T_NV, C_pei3T_NV, C_usa6T_NV, C_usa3T_NV, C_usa7T_NV, C_usa4T_NV,
                  C_pei4_NR, C_usa6_NR, C_usa7_NR, C_pei2T_NR, C_pei4T_NR, C_pei3T_NR, C_usa6T_NR, C_usa3T_NR, C_usa7T_NR, C_usa4T_NR)) %>%
        mutate(pei_VAF = (C_pei1_NV + C_pei2_NV + C_pei3_NV) / (C_pei1_NR+ C_pei2_NR + C_pei3_NR),
               usa_VAF = (C_usa1_NV + C_usa2_NV + C_usa3_NV + C_usa4_NV + C_usa5_NV) / (C_usa1_NR + C_usa2_NR + C_usa3_NR + C_usa4_NR + C_usa5_NR))

# save RDS file to save time loading later
    saveRDS(snvs_LOH2, file = "snvs_LOH.RDS")
    # readRDS(file = "snvs_LOH.RDS")

#BINNING - same as for signature extraction
    # Is SNV found in any cancer
        anyC <- filter(snvs_LOH2, C_pei1_NV > C_pei1_avg/stringent_depth |
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

    sublin_specific <- rbind(allUSAnoPEInoH, allPEInoUSAnoH, allUSAnoPEIanyH, allPEInoUSAanyH)
    usa_specific <- rbind(allUSAnoPEInoH, allUSAnoPEIanyH)
    pei_specific <- rbind(allPEInoUSAnoH, allPEInoUSAanyH)

# how much of genome are we excluding with LOH calls
    USA_LOH_genome_size <- read.table("/ssd3/Mar_genome_analysis/LOH/july_2021/output/new/USA_hetero_LOH_10_counts.bed") %>% mutate(V4 = V3-V2) %>% .$V4 %>%as.numeric() %>% sum()
    PEI_LOH_genome_size <- read.table("/ssd3/Mar_genome_analysis/LOH/july_2021/output/new/PEI_hetero_LOH_10_counts.bed") %>% mutate(V4 = V3-V2) %>% .$V4 %>%as.numeric() %>% sum()
    BOTH_LOH_genome_size <- read.table("/ssd3/Mar_genome_analysis/LOH/july_2021/output/new/BOTH_SUBLINEAGES_LOH_10_counts.bed") %>% mutate(V4 = V3-V2) %>% .$V4 %>%as.numeric() %>% sum()
    full_genome_size <- read.table("/ssd3/Mar_genome_analysis/genomes/Mar.3.4.6.p1_Q30Q30A.fasta.fai") %>% .$V2 %>%as.numeric() %>% sum()
    USAnonLOH_genome_size = full_genome_size - USA_LOH_genome_size 
    PEInonLOH_genome_size = full_genome_size - PEI_LOH_genome_size
    BOTHnonLOH_genome_size = full_genome_size - BOTH_LOH_genome_size
    USA_LOH_genome_size # 98,151,754
    PEI_LOH_genome_size # 155,373,786 
    BOTH_LOH_genome_size # 237,564,314
    USAnonLOH_genome_size # 1,117,916,537
    PEInonLOH_genome_size # 1,060,694,505
    BOTHnonLOH_genome_size # 978,503,977



genome_LOH_plot <- data.frame(regions = c("Full genome", "Not LOH in either", "LOH in either", "LOH in PEI", "LOH in USA"),bp = c(full_genome_size, BOTHnonLOH_genome_size, BOTH_LOH_genome_size, PEI_LOH_genome_size, USA_LOH_genome_size))
genome_LOH_plot <- data.frame(regions = c("Full genome", "LOH in PEI", "LOH in USA"),bp = c(full_genome_size, PEI_LOH_genome_size, USA_LOH_genome_size))

anyH_usaLOH <- rbind(allUSAnoPEIanyH,allPEInoUSAanyH) %>%
    count(bin, USA_LOH) %>% 
    mutate(mu_per_mb = case_when(USA_LOH == "LOH" ~ (n / USA_LOH_genome_size * 1000000),
                                 USA_LOH == "nonLOH"  ~ (n / USAnonLOH_genome_size * 1000000)),
           LOH = case_when(USA_LOH == "LOH" ~ "usaLOH",
                               USA_LOH == "nonLOH"  ~ "usa_nonLOH")) %>%
    select(-USA_LOH)
                                        
anyH_usaLOH
anyH_peiLOH <- rbind(allUSAnoPEIanyH,allPEInoUSAanyH) %>% 
    count(bin, PEI_LOH) %>%
    mutate(mu_per_mb = case_when(PEI_LOH == "LOH" ~ (n / PEI_LOH_genome_size * 1000000),
                                 PEI_LOH == "nonLOH"  ~ (n / PEInonLOH_genome_size * 1000000)),
           LOH = case_when(PEI_LOH == "LOH" ~ "peiLOH",
                               PEI_LOH == "nonLOH"  ~ "pei_nonLOH")) %>%                                 
    select(-PEI_LOH)           
plot_output <- rbind(anyH_peiLOH, anyH_usaLOH)

plot_output[9,] <- list("allPEInoUSAanyH", (plot_output[1,"n"]+plot_output[2,"n"]), (plot_output[1,2]+plot_output[2,2])/full_genome_size*1000000, "full_genome")
plot_output[10,] <- list("allUSAnoPEIanyH", (plot_output[3,2]+plot_output[4,2]), ((plot_output[3,2]+plot_output[4,2])/full_genome_size*1000000), "full_genome")



    p1 <- ggplot(genome_LOH_plot, aes(x = regions, y=bp/1000000)) + 
        geom_bar(stat = "identity",fill = c("black", "red","blue"))+ # c("grey","grey","grey", "red","blue")
        geom_text(aes(label=round(bp/1000000)), vjust=-0.25, size =3)+
        ylab("Megabases")+
        theme_classic()   +
        theme(axis.title.x=element_blank(),
              axis.text=element_text(size=6) )
    p2 <- ggplot(plot_output[c(3:6,9,10),], aes(x = bin, y=mu_per_mb, fill = LOH)) + 
        geom_bar(stat = "identity", position = "dodge", width = 0.5)+
        geom_text(aes(label=round(mu_per_mb)),position=position_dodge(width=0.5), vjust=-0.25, size =3)+
        ylab("SNVs per Mb")+
        theme_classic()+
        theme(axis.title.x=element_blank(),
              axis.text=element_text(size=6))

pdf("LOH_supp_plot.pdf", width = 7, height = 3.5)
grid.arrange(p1,p2, nrow = 1, widths=c(1,2))
dev.off()



# Merge VAFs and call homozygous
    allCanyH_P <- allCanyH %>% 
        mutate(VAF = pei_VAF, sub = "PEI") %>%
        select(USA_LOH,PEI_LOH, VAF, sub)
    allCanyH_U <- allCanyH %>% 
        mutate(VAF = usa_VAF, sub = "USA") %>%
        select(USA_LOH,PEI_LOH, VAF, sub)
    allCanyH_UP <- rbind(allCanyH_P,allCanyH_U) %>%
        mutate(homozygous = case_when(VAF > 0.8 ~ TRUE,
                                    VAF <= 0.8 ~ FALSE))
    count(allCanyH_UP, sub, PEI_LOH,homozygous)
    count(allCanyH_UP, sub, USA_LOH,homozygous)
# >     count(allCanyH_UP, sub, PEI_LOH,homozygous)
#   sub PEI_LOH homozygous       n
# 1 PEI     LOH      FALSE  200838
# 2 PEI     LOH       TRUE  513851
# 3 PEI  nonLOH      FALSE 4224011
# 4 PEI  nonLOH       TRUE 3129253
# 5 USA     LOH      FALSE  361938
# 6 USA     LOH       TRUE  352751
# 7 USA  nonLOH      FALSE 4128385
# 8 USA  nonLOH       TRUE 3224879
# >     count(allCanyH_UP, sub, USA_LOH,homozygous)
#   sub USA_LOH homozygous       n
# 1 PEI     LOH      FALSE  264133
# 2 PEI     LOH       TRUE  265582
# 3 PEI  nonLOH      FALSE 4160716
# 4 PEI  nonLOH       TRUE 3377522
# 5 USA     LOH      FALSE  154142
# 6 USA     LOH       TRUE  375573
# 7 USA  nonLOH      FALSE 4336181
# 8 USA  nonLOH       TRUE 3202057

PEI_LOH_HOMOZ <- 513851/(513851+200838) # 0.7189855
USA_LOH_HOMOZ <- 375573/(375573+154142) # 0.7090096
PEI_nonLOH_HOMOZ <- 3129253/(3129253+4224011) # 0.4255597
USA_nonLOH_HOMOZ <- 3202057/(3202057+4336181) # 0.4247753

pdf("LOH_allCanyH_overlay.pdf")
    ggplot(plot_output, aes(x = bin, y=mu_per_mb, fill = LOH)) + 
        geom_bar(stat = "identity", position = "dodge", width = 0.5)+
        geom_text(aes(label=round(mu_per_mb)),position=position_dodge(width=0.5), vjust=-0.25)+
        theme_classic()
        plot1 <- ggplot(allCanyH_UP)+
            geom_freqpoly(aes(x=VAF, color = sub), binwidth = 0.01)+
            geom_hline(yintercept=0)+
            scale_color_manual(values=c("red", "blue"))+
            xlim(0, 1.01)+
            ggtitle("allCanyH for PEI LOH") +
            theme_classic() +
            facet_grid(rows = vars(PEI_LOH), scales = "free")
        plot2 <- ggplot(allCanyH_UP)+
            geom_freqpoly(aes(x=VAF, color = sub), binwidth = 0.01)+
            geom_hline(yintercept=0)+
            scale_color_manual(values=c("red", "blue"))+
            xlim(0, 1.01)+
            ggtitle("allCanyH for USA LOH") +
            theme_classic() +
            facet_grid(rows = vars(USA_LOH), scales = "free")
    plot1
    plot2
dev.off()
