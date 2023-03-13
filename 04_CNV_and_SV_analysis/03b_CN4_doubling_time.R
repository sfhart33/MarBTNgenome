# File notes found here: C:\Users\shart\Metzger Lab Dropbox\Sam_data\Mya_genome\revision\doubling_timing.r    
    setwd("/ssd3/Mar_genome_analysis/revision/WGD")
    library(tidyverse)



# allCanyH germline SNVs generated here: /ssd3/Mar_genome_analysis/CNV_calling/FREQ/README
    PEI_GERM <- read.delim("/ssd3/Mar_genome_analysis/CNV_calling/FREQ/OUTPUT/PEI_allC_anyH.txt")
    USA_GERM <- read.delim("/ssd3/Mar_genome_analysis/CNV_calling/FREQ/OUTPUT/USA_allC_anyH.txt")
    GERM <- inner_join(PEI_GERM, USA_GERM, by = c("chr", "pos"))

# Mode function
    getmode <- function(v) {
        uniqv <- unique(v)
        uniqv[which.max(tabulate(match(v, uniqv)))]
    }
    ratiocalc <- function(x) {
        low <- length(x[x < 0.375])
        high <- length(x[x > 0.375 & x < 0.625])
        total <- length(x)
        high/total
    }
# PEI and USA somatic sub-lineage specific SNVs generated in CNV_calling_pipeline

    PEI_medians <- PEI_SNVs %>%
        mutate(segment = paste0(chr,"_",round(as.numeric(start+1), -5))) %>%
        filter(CN != "8plus", pei_f < 0.9) %>%
        group_by(segment) %>%
        summarize(PEI_median = median(pei_f, na.rm = TRUE),
                  PEI_mean = mean(pei_f, na.rm = TRUE),
                  PEI_mode = 0.25*getmode(round(pei_f*4,1)), # rounded to nearest 0.025
                  PEI_ratio = ratiocalc(pei_f),
                  PEI_count = n())
    USA_medians <- USA_SNVs %>%
        mutate(segment = paste0(chr,"_",round(as.numeric(start+1), -5))) %>%
        filter(CN != "8plus", usa_f < 0.9) %>%
        group_by(segment) %>%
        summarize(USA_median = median(usa_f, na.rm = TRUE),
                  USA_mean = mean(usa_f, na.rm = TRUE),
                  USA_mode = 0.25*getmode(round(usa_f*4,1)), # rounded to nearest 0.025
                  USA_ratio = ratiocalc(usa_f),
                  USA_count = n())
    GERM_medians <- GERM %>%
        mutate(segment = paste0(chr,"_",round(as.numeric(pos), -5))) %>%
        filter(USA_CN != "8plus", PEI_CN != "8plus", PEI_freq < 0.9, USA_freq < 0.9) %>%
        group_by(segment) %>%
        summarize(PEI_allC_median = median(PEI_freq, na.rm = TRUE),
                  USA_allC_median = median(USA_freq, na.rm = TRUE),
                  PEI_allC_mean = mean(PEI_freq, na.rm = TRUE),
                  PEI_allC_mode = 0.25*getmode(round(PEI_freq*4,1)), # rounded to nearest 0.025
                  USA_allC_mean = mean(USA_freq, na.rm = TRUE),
                  USA_allC_mode = 0.25*getmode(round(USA_freq*4,1)), # rounded to nearest 0.025
                  allC_count = n(),
                  USA_CN = median(as.numeric(USA_CN), na.rm = TRUE),
                  PEI_CN = median(as.numeric(PEI_CN), na.rm = TRUE))
    head(GERM_medians)
    filter(GERM_medians, USA_CN == 4, PEI_CN == 4) %>% nrow()

    SOMATIC <- full_join(PEI_medians, USA_medians, by = "segment")

    ALL <- full_join(GERM_medians, SOMATIC)
    ALL_CN4 <- filter(ALL, USA_CN == 4, PEI_CN == 4) 
    as.data.frame(ALL_CN4) %>% head()
# Shared CN4 regions
    PEI_SNVs_CN4 <- PEI_SNVs %>%
        mutate(segment = paste0(chr,"_",round(as.numeric(start+1), -5))) %>%
        filter(segment %in% ALL_CN4$segment, CN == 4)
    USA_SNVs_CN4 <- USA_SNVs %>%
        mutate(segment = paste0(chr,"_",round(as.numeric(start+1), -5))) %>%
        filter(segment %in% ALL_CN4$segment, CN == 4)
    GERM_SNVs_CN4 <- GERM %>%
        mutate(segment = paste0(chr,"_",round(as.numeric(pos), -5))) %>%
        filter(segment %in% ALL_CN4$segment, PEI_CN == 4, USA_CN == 4, PEI_freq < 0.9, USA_freq < 0.9)

    
    pdf("Shared CN4 regions.pdf")
        # ggplot(ALL_CN4) +
        #     geom_freqpoly(aes(USA_allC_median), color="blue", binwidth=0.025) +
        #     geom_freqpoly(aes(PEI_allC_median), color="red", binwidth=0.025) +
        #     geom_freqpoly(aes(USA_median), color="blue", binwidth=0.025, linetype = "dotted") +
        #     geom_freqpoly(aes(PEI_median), color="red", binwidth=0.025, linetype = "dotted") +
        #     theme_classic() +
        #     theme(axis.text=element_text(size=12,face="bold"),
        #             axis.title=element_text(size=16,face="bold")) +
        #     scale_x_continuous(name="Median allele freq (100kB windows)", breaks=c(0,0.25,0.5,0.75), limits = c(0,0.9))
        # ggplot(ALL_CN4) +
        #     geom_freqpoly(aes(USA_allC_mean), color="blue", binwidth=0.025) +
        #     geom_freqpoly(aes(PEI_allC_mean), color="red", binwidth=0.025) +
        #     geom_freqpoly(aes(USA_mean), color="blue", binwidth=0.025, linetype = "dotted") +
        #     geom_freqpoly(aes(PEI_mean), color="red", binwidth=0.025, linetype = "dotted") +
        #     theme_classic() +
        #     theme(axis.text=element_text(size=12,face="bold"),
        #             axis.title=element_text(size=16,face="bold")) +
        #     scale_x_continuous(name="Mean allele freq (100kB windows)", breaks=c(0,0.25,0.5,0.75), limits = c(0,0.9))
        # ggplot(ALL_CN4) +
        #     geom_freqpoly(aes(USA_allC_mode), color="blue", binwidth=0.025) +
        #     geom_freqpoly(aes(PEI_allC_mode), color="red", binwidth=0.025) +
        #     geom_freqpoly(aes(USA_mode), color="blue", binwidth=0.025, linetype = "dotted") +
        #     geom_freqpoly(aes(PEI_mode), color="red", binwidth=0.025, linetype = "dotted") +
        #     theme_classic() +
        #     theme(axis.text=element_text(size=12,face="bold"),
        #             axis.title=element_text(size=16,face="bold")) +
        #     scale_x_continuous(name="Mode allele freq (100kB windows)", breaks=c(0,0.25,0.5,0.75), limits = c(0,0.9))
        ggplot() +
            geom_freqpoly(data = GERM_SNVs_CN4, aes(USA_freq, y = stat(density)), color="blue", binwidth=0.01) +
            geom_freqpoly(data = GERM_SNVs_CN4, aes(PEI_freq, y = stat(density)), color="red", binwidth=0.01) +
            geom_freqpoly(data = USA_SNVs_CN4, aes(usa_f, y = stat(density)), color="blue", binwidth=0.01, linetype = "dotted") +
            geom_freqpoly(data = PEI_SNVs_CN4, aes(pei_f, y = stat(density)), color="red", binwidth=0.01, linetype = "dotted") +
            theme_classic() +
            theme(axis.text=element_text(size=12,face="bold"),
                    axis.title=element_text(size=16,face="bold")) +
            scale_x_continuous(name="Variant allele frequency", breaks=c(0,0.25,0.5,0.75), limits = c(0,0.89))      
        ggplot() +
            geom_freqpoly(data = GERM_SNVs_CN4, aes(USA_freq), color="blue", binwidth=0.01) +
            geom_freqpoly(data = GERM_SNVs_CN4, aes(PEI_freq), color="red", binwidth=0.01) +
            theme_classic() +
            theme(axis.text=element_text(size=12,face="bold"),
                    axis.title=element_text(size=16,face="bold")) +
            scale_x_continuous(name="Germline variant allele frequency", breaks=c(0,0.25,0.5,0.75), limits = c(0,0.89))      
        ggplot() +
            geom_freqpoly(data = USA_SNVs_CN4, aes(usa_f), color="blue", binwidth=0.01) +
            geom_freqpoly(data = PEI_SNVs_CN4, aes(pei_f), color="red", binwidth=0.01) +
            theme_classic() +
            theme(axis.text=element_text(size=12,face="bold"),
                    axis.title=element_text(size=16,face="bold")) +
            scale_x_continuous(name="Somatic mutation allele frequency", breaks=c(0,0.25,0.5,0.75), limits = c(0,0.89))      
        ggplot(ALL_CN4) +
            geom_freqpoly(aes(USA_ratio), color="blue", binwidth=0.05) +
            geom_freqpoly(aes(PEI_ratio), color="red", binwidth=0.05) +
            theme_classic() +
            theme(axis.text=element_text(size=12,face="bold"),
                    axis.title=element_text(size=16,face="bold")) +
            scale_x_continuous(name="Fraction pre-duplication somatic mutations (2/4 frequency)", breaks=c(0,0.2,0.4,0.6,0.8,1), limits = c(0,1))     
        ggplot(GERM_SNVs_CN4[sample(nrow(GERM_SNVs_CN4), 100000), ], aes(USA_freq, PEI_freq)) +
            #geom_density_2d_filled(alpha = 0.5) +
            geom_point(size=0.1) +
            theme_classic() +
            xlab("Variant allele frequency in USA")+
            ylab("Variant allele frequency in PEI")+
            theme(axis.text=element_text(size=12,face="bold"),
                    axis.title=element_text(size=16,face="bold")) +
            scale_x_continuous(breaks=c(0,0.25,0.5,0.75), limits = c(0,0.9)) +
            scale_y_continuous(breaks=c(0,0.25,0.5,0.75), limits = c(0,0.9))
    dev.off()

# calculate timing estimate
usa_low <- filter(USA_SNVs_CN4, usa_f < 0.375) %>% nrow()
usa_high <- filter(USA_SNVs_CN4, usa_f > 0.375, usa_f < 0.625) %>% nrow()
usa_total <- nrow(USA_SNVs_CN4)
pei_low <- filter(PEI_SNVs_CN4, pei_f <  0.375) %>% nrow()
pei_high <- filter(PEI_SNVs_CN4, pei_f >  0.375, pei_f < 0.625) %>% nrow()
pei_total <- nrow(PEI_SNVs_CN4)
usa_low/usa_total # 0.83
pei_low/pei_total # 0.71
usa_high/usa_total # 0.14
pei_high/pei_total # 0.26
