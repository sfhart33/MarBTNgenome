# Original file here: C:\Users\shart\Metzger Lab Dropbox\Sam_data\Mya_genome\LOH\somatypus_output_LOH_upstream.r

library(tidyverse)
library(zoo)


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

setwd("/ssd3/Mar_genome_analysis/LOH/july_2021")

# Get high confidence germline alleles: found in a healthy clam and in all cancer samples (with high stringency value)
    #Set thresholds
    min_cutoff=3
    stringent_depth=6 # was 8
    weak_depth=16
    allCanyH_stringent <- filter(snvs, 
                        (H_ref_NV > H_pei_avg/stringent_depth |
                        H_usa_NV > H_pei_avg/stringent_depth |
                        H_pei_NV > H_pei_avg/stringent_depth) &
                        C_pei1_NV > C_pei1_avg/stringent_depth &
                        C_pei2_NV > C_pei2_avg/stringent_depth &
                        C_pei3_NV > C_pei3_avg/stringent_depth &
                        C_usa1_NV > C_usa1_avg/stringent_depth &
                        C_usa2_NV > C_usa2_avg/stringent_depth &
                        C_usa3_NV > C_usa3_avg/stringent_depth &
                        C_usa4_NV > C_usa4_avg/stringent_depth &
                        C_usa5_NV > C_usa5_avg/stringent_depth)
    nrow(allCanyH_stringent)

#plot to determine threshold at which to consider an allele homozygous
    allCanyH_freq <- mutate(allCanyH_stringent, pei1 = C_pei1_NV/C_pei1_NR,
                            pei2 = C_pei2_NV/C_pei2_NR,
                            pei3 = C_pei3_NV/C_pei3_NR,
                            usa1 = C_usa1_NV/C_usa1_NR,
                            usa2 = C_usa2_NV/C_usa2_NR,
                            usa3 = C_usa3_NV/C_usa3_NR,
                            usa4 = C_usa4_NV/C_usa4_NR,
                            usa5 = C_usa5_NV/C_usa5_NR ) %>%
                    mutate(usaAVG = (usa1+usa2+usa3+usa4+usa5)/5,
                            peiAVG = (pei1+pei2+pei3)/3)
    pdf("frequencies_stringent.pdf")
    for(i in c("pei1","pei2","pei3","usa1","usa2","usa3","usa4","usa5","usaAVG","peiAVG")){
        plot1 <- ggplot(allCanyH_freq, aes(get(i)))+
            geom_histogram(binwidth=0.01)+
            xlab("Allele freq")+
            ylab("count")+
            xlim(0,1) +
            theme_classic() +
                theme(axis.text=element_text(size=12,face="bold"),
                axis.title=element_text(size=16,face="bold"),
                text=element_text(size=18,face="bold")) +
            ggtitle(i)
        print(plot1)
    }
    dev.off()
	
# Transform table to call homozygous/heterozygous sites and potential LOH sites when discordant between USA and PEI
    threshold=0.8
    allCanyH_LOH <- mutate(allCanyH_stringent, 
                            pei1 = C_pei1_NV/C_pei1_NR,
                            pei2 = C_pei2_NV/C_pei2_NR,
                            pei3 = C_pei3_NV/C_pei3_NR,
                            usa1 = C_usa1_NV/C_usa1_NR,
                            usa2 = C_usa2_NV/C_usa2_NR,
                            usa3 = C_usa3_NV/C_usa3_NR,
                            usa4 = C_usa4_NV/C_usa4_NR,
                            usa5 = C_usa5_NV/C_usa5_NR ) %>%
                    mutate(P_homoz = ifelse(pei1>threshold & pei2>threshold & pei3>threshold,
                                        1,0),
                            P_hetero = ifelse(pei1<threshold & pei2<threshold & pei3<threshold,
                                        1,0),
                            P_LOH = ifelse(pei1>threshold & pei2>threshold & pei3>threshold & 
                                        usa1<threshold & usa3<threshold & usa5<threshold, #& usa2<threshold & usa4<threshold ONLY USE THREE USAs to not bias the calls
                                        1,0),
                            U_homoz = ifelse(usa1>threshold & usa3>threshold & usa5>threshold, #& usa2>threshold  & usa4>threshold
                                        1,0),
                            U_hetero = ifelse(usa1<threshold & usa3<threshold & usa5<threshold, #& usa2<threshold  & usa4<threshold
                                        1,0),
                            U_LOH = ifelse(pei1<threshold & pei2<threshold & pei3<threshold & 
                                        usa1>threshold & usa3>threshold & usa5>threshold, #& usa2>threshold  & usa4>threshold
                                        1,0),
                                ) %>%
                    select(CHROM,POS,P_homoz,U_homoz,P_hetero,U_hetero,P_LOH, U_LOH)
    head(allCanyH_LOH)				

# Filter for just heterozygous snps in each sublineage
    USA_hetero <- filter(allCanyH_LOH, U_hetero == 1) %>%
        select(CHROM,POS,P_LOH)
    PEI_hetero <- filter(allCanyH_LOH, P_hetero == 1) %>%
        select(CHROM,POS,U_LOH)
    head(USA_hetero)

# Calculate windows of 100 SNVs, counting occurences
    bins=100
    bins2=50
    for(i in 0:16){
        print(paste0("Mar.3.4.6.p1_scaffold", i), quote = FALSE)
        allC_freqX <- filter(allCanyH_LOH, CHROM == paste0("Mar.3.4.6.p1_scaffold",i))
        USA_heteroX <- filter(USA_hetero, CHROM == paste0("Mar.3.4.6.p1_scaffold",i))
        PEI_heteroX <- filter(PEI_hetero, CHROM == paste0("Mar.3.4.6.p1_scaffold",i))
        allC_freqX <-  allC_freqX %>%
            mutate(p_homoz_count = rollapply(P_homoz[1:(nrow(allC_freqX)+bins-1)],bins,sum,align="left"),
                u_homoz_count = rollapply(U_homoz[1:(nrow(allC_freqX)+bins-1)],bins,sum,align="left"),
                p_hetero_count = rollapply(P_hetero[1:(nrow(allC_freqX)+bins-1)],bins,sum,align="left"),
                u_hetero_count = rollapply(U_hetero[1:(nrow(allC_freqX)+bins-1)],bins,sum,align="left"),
                p_loh_count = rollapply(P_LOH[1:(nrow(allC_freqX)+bins-1)],bins,sum,align="left"),
                u_loh_count = rollapply(U_LOH[1:(nrow(allC_freqX)+bins-1)],bins,sum,align="left"),
                start = POS-1,
                end = lead(POS[1:(nrow(allC_freqX))], (bins-1))+1) 
        USA_heteroX <-  USA_heteroX %>%
            mutate(loh_count = rollapply(P_LOH[1:(nrow(USA_heteroX)+bins2-1)],bins2,sum,align="left"),
                    start = POS-1,
                    end = lead(POS[1:(nrow(USA_heteroX))], (bins2-1))+1) 
        PEI_heteroX <-  PEI_heteroX %>%
            mutate(loh_count = rollapply(U_LOH[1:(nrow(PEI_heteroX)+bins2-1)],bins2,sum,align="left"),
                    start = POS-1,
                    end = lead(POS[1:(nrow(PEI_heteroX))], (bins2-1))+1) 
            if(i == 0){
            allCanyH_LOH2 <- allC_freqX
            USA_hetero2 <- USA_heteroX
            PEI_hetero2 <- PEI_heteroX
        }
        if(i > 0){
            allCanyH_LOH2 <- rbind(allCanyH_LOH2, allC_freqX)
            USA_hetero2 <- rbind(USA_hetero2, USA_heteroX)
            PEI_hetero2 <- rbind(PEI_hetero2, PEI_heteroX)
        }
    }


# plot to check distribution is even between usa and pei (why I only use 3 usa samples)
    pdf("homozygousity_PEIvsUSA.pdf")
        ggplot(allCanyH_LOH2, aes(p_homoz_count - u_homoz_count))+
            geom_histogram(binwidth=1)+
            xlab("p_homoz_count-u_homoz_count")+
            ylab("count")+
            xlim(-100,100) +
            theme_classic() +
                theme(axis.text=element_text(size=12,face="bold"),
                axis.title=element_text(size=16,face="bold"),
                text=element_text(size=18,face="bold")) +
            ggtitle("p_homoz_count-u_homoz_count")
    dev.off()


# plot to determine best cutoff thresholds
pdf("homoz_hetero_loh_counts_stringent.pdf")
    for(i in c("p_homoz_count","u_homoz_count","p_hetero_count","u_hetero_count","p_loh_count","u_loh_count")){
        plot1 <- ggplot(allCanyH_LOH2, aes(get(i)))+
            geom_histogram(binwidth=1)+
            xlab(i)+
            ylab("count")+
            xlim(1,100) +
            theme_classic() +
                theme(axis.text=element_text(size=12,face="bold"),
                axis.title=element_text(size=16,face="bold"),
                text=element_text(size=18,face="bold")) +
            ggtitle(i)
        print(plot1)
    }
    ggplot(allCanyH_LOH2, aes(u_loh_count/p_hetero_count))+
        geom_histogram(binwidth=0.01)+
        xlab("fraction of PEI hetero SNVs that are homoz in USA")+
        ylab("count")+
        xlim(-0.01,1.01) +
        theme_classic() +
            theme(axis.text=element_text(size=12,face="bold"),
            axis.title=element_text(size=16,face="bold"),
            text=element_text(size=18,face="bold")) +
        ggtitle("fraction of PEI hetero SNVs that are homoz in USA")
    ggplot(allCanyH_LOH2, aes(p_loh_count/u_hetero_count))+
        geom_histogram(binwidth=0.01)+
        xlab("fraction of USA hetero SNVs that are homoz in PEI")+
        ylab("count")+
        xlim(-0.01,1.01) +
        theme_classic() +
            theme(axis.text=element_text(size=12,face="bold"),
            axis.title=element_text(size=16,face="bold"),
            text=element_text(size=18,face="bold")) +
        ggtitle("fraction of USA hetero SNVs that are homoz in PEI")
    ggplot(allCanyH_LOH2, aes(u_loh_count/p_hetero_count))+
        geom_histogram(binwidth=0.01)+
        xlab("fraction of PEI hetero SNVs that are homoz in USA")+
        ylab("count")+
        xlim(0.01,1.01) +
        theme_classic() +
            theme(axis.text=element_text(size=12,face="bold"),
            axis.title=element_text(size=16,face="bold"),
            text=element_text(size=18,face="bold")) +
        ggtitle("fraction of PEI hetero SNVs that are homoz in USA")
    ggplot(allCanyH_LOH2, aes(p_loh_count/u_hetero_count))+
        geom_histogram(binwidth=0.01)+
        xlab("fraction of USA hetero SNVs that are homoz in PEI")+
        ylab("count")+
        xlim(0.01,1.01) +
        theme_classic() +
            theme(axis.text=element_text(size=12,face="bold"),
            axis.title=element_text(size=16,face="bold"),
            text=element_text(size=18,face="bold")) +
        ggtitle("fraction of USA hetero SNVs that are homoz in PEI")
    dev.off()

# Plot heterozygous windows for # homozygous in the other
    pdf("hetero_to_homoz_fractions.pdf")
    ggplot(USA_hetero2, aes(loh_count/bins2))+
        geom_histogram(binwidth=1/bins2)+
        xlab("% PEI homoz per 50 USA hetero")+
        ylab("count")+
        # xlim(0.01,1.01) +
        theme_classic() +
            theme(axis.text=element_text(size=12,face="bold"),
            axis.title=element_text(size=16,face="bold"),
            text=element_text(size=18,face="bold")) +
        ggtitle("% PEI homoz per 50 USA hetero")
    ggplot(PEI_hetero2, aes(loh_count/bins2))+
        geom_histogram(binwidth=1/bins)+
        xlab("% USA homoz per 50 PEI hetero")+
        ylab("count")+
        # xlim(0.01,1.01) +
        theme_classic() +
            theme(axis.text=element_text(size=12,face="bold"),
            axis.title=element_text(size=16,face="bold"),
            text=element_text(size=18,face="bold")) +
        ggtitle("% USA homoz per 50 PEI hetero")
    dev.off()

# Preious testing of various calling methods
    # ratio_threshold=0
    # for(ratio_threshold in (1:19)/20){
    #     P_LOH <- filter(allCanyH_LOH2,
    #                     p_loh_count/u_hetero_count > ratio_threshold)
    #     print(nrow(P_LOH)/nrow(allCanyH_LOH2))
    # }
    # for(ratio_threshold in (1:19)/20){
    #     U_LOH <- filter(allCanyH_LOH2,
    #                     u_loh_count/p_hetero_count > ratio_threshold)
    #     print(nrow(U_LOH)/nrow(allCanyH_LOH2))
    # }
    # homoz_threshold=0.7 # 0.7
    # loh_threshold=5 # 5
    # for(loh_threshold in 0:30){
    #     P_LOH <- filter(allCanyH_LOH2,
    #                     #p_homoz_count > bins*homoz_threshold &
    #                     p_homoz_count > u_homoz_count &
    #                     p_loh_count > loh_threshold &
    #                     p_loh_count > u_loh_count
    #                     ) %>%
    #         select(CHROM,start,end,p_homoz_count,u_homoz_count,p_loh_count,u_loh_count)
    #     print(nrow(P_LOH)/nrow(allCanyH_LOH2))
    # }
    # for(loh_threshold in 0:30){ # (25:45)/50
    #     U_LOH <- filter(allCanyH_LOH2,
    #                     #u_homoz_count > bins*homoz_threshold &
    #                     u_homoz_count > p_homoz_count &
    #                     u_loh_count > loh_threshold &
    #                     u_loh_count > p_loh_count
    #                     ) %>%
    #         select(CHROM,start,end,p_homoz_count,u_homoz_count,p_loh_count,u_loh_count)
    #     print(nrow(U_LOH)/nrow(allCanyH_LOH2))
    # }

# Test using hetero calls
    #threshold=0
    for(threshold in c(1:5,(1:5)*10)){
        P_LOH <- filter(USA_hetero2,
                        loh_count > threshold)
        print(nrow(P_LOH)/nrow(USA_hetero2))
    }
    for(threshold in c(1:5,(1:5)*10)){
        U_LOH <- filter(PEI_hetero2,
                        loh_count > threshold)
        print(nrow(U_LOH)/nrow(PEI_hetero2))
    }

# Output bed file of all likely LOH regions
    LOH_call_to_bed <- function(loh_threshold){
        allCanyH_LOH2 %>% 
            filter(u_loh_count > loh_threshold-1) %>%
            select(CHROM,start,end) %>%
            write.table(file = paste0("USA_LOH_",loh_threshold,"_counts.bed"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep = '\t')
        allCanyH_LOH2 %>% 
            filter(p_loh_count > loh_threshold-1) %>%
            select(CHROM,start,end) %>%
            write.table(file = paste0("PEI_LOH_",loh_threshold,"_counts.bed"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep = '\t')	
    }
    LOH_call_to_bed2 <- function(loh_threshold){
        USA_hetero2 %>% 
            filter(loh_count > loh_threshold-1) %>%
            select(CHROM,start,end) %>%
            write.table(file = paste0("PEI_LOH_",loh_threshold,"_hetero_counts.bed"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep = '\t')	
        PEI_hetero2 %>% 
            filter(loh_count > loh_threshold-1) %>%
            select(CHROM,start,end) %>%
            write.table(file = paste0("USA_LOH_",loh_threshold,"_hetero_counts.bed"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep = '\t')	
    }

    # for(i in 0:50){
    #     LOH_call_to_bed(i)
    # }

    for(i in 1:51){ # c(1:5,(1:5)*10)
            LOH_call_to_bed2(i)
    }

