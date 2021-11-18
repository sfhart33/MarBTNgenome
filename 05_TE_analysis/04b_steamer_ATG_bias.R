#Original file here: C:\Users\shart\Metzger Lab Dropbox\Sam_data\Mya_genome\steamer\atg_bias.r

library(tidyverse)
library(ggseqlogo)
setwd("/ssd3/Mar_genome_analysis/steamer/final_pipeline/logos")

data(ggseqlogo_sample)

# GC content of the genome WARNING - OLD VERSION OF GENOME - NEW VERSION SHOULD BE SIMILAR BUT NEED TO CHECK
    GC = 0.176471649081384
    AT = 0.323528350918616

# Read info from .sh script and restructure
    downstream_all <- read.table("steamer_bias_down.txt",sep="\t", header=T, check.names=F) %>%
        mutate(freq = count/total)
    upstream_all <- read.table("steamer_bias_up.txt",sep="\t", header=T, check.names=F) %>%
        mutate(freq = count/total)
    genome_all <- read.table("steamer_bias_genome.txt",sep="\t", header=T, check.names=F) %>%
        mutate(freq = count/total)
    head(genome_all)
    downstream_all$exp_freq <- rep(c(AT,AT,GC,GC),nrow(downstream_all)/4)
    upstream_all$exp_freq <- rep(c(AT,AT,GC,GC),nrow(upstream_all)/4)
    genome_all$exp_freq <- rep(c(AT,AT,GC,GC),nrow(genome_all)/4)
    downstream_all <- mutate(downstream_all, dif = freq-exp_freq)
    upstream_all <- mutate(upstream_all, dif = freq-exp_freq)
    genome_all <- mutate(genome_all, dif = freq-exp_freq)

# Plotting tests
    pdf("logo_plot_tests.pdf")
    for(updown in c("up","down")){
        for(samples in c("MELC-2E11", "MELC-A9", "PEI-DF490", "FFM-19G1", "FFM-20B2", "FFM-22F10", "MELC-A11_S1", "NYTC-C9_S2", "DF-488", "DN-HL03", "PEI-DN08_S3")){
            ffm_test <- filter(get(paste0(updown,"stream_all")), sample == samples)
            nrow(ffm_test) %>% print()
            test_data <- data.frame(A=filter(ffm_test, nuc == "A")$dif,
                                        T=filter(ffm_test, nuc == "T")$dif,
                                        G=filter(ffm_test, nuc == "G")$dif,
                                        C=filter(ffm_test, nuc == "C")$dif) %>%
                                        t()
            #colnames(test_data) <- 1:20
            
            plot1 <- ggseqlogo(test_data, method = 'custom' ) +
                        ggtitle(paste0(samples,": ",updown,"stream")) +
                        xlab("Position relative to Steamer insertion") +
                        ylab("Probability vs expected")
            print(plot1)
        }
    }
    dev.off()
    pdf("logo_plot_genome.pdf")
    for(samples in c("all_sites", "anyCnoH", "allCnoH_allPEI_allUSA", "allCnoH", "allUSA", "allPEI")){
        print(samples)
        sample_plot <- filter(genome_all, sample == samples)
        #nrow(sample_plot) %>% print()
        test_data <- data.frame(A=filter(sample_plot, nuc == "A")$dif,
                                    T=filter(sample_plot, nuc == "T")$dif,
                                    G=filter(sample_plot, nuc == "G")$dif,
                                    C=filter(sample_plot, nuc == "C")$dif) %>%
                                    t()
        plot1 <- ggseqlogo(test_data, method = 'custom' ) +
                    ggtitle(paste0(samples,": retrieved from genome")) +
                    xlab("Position (Steamer 5bp dup: 16-20)") +
                    ylab("Probability vs expected")
        print(plot1)
    }
    dev.off()


# Final plot for paper
    steamer_to_plot <- filter(genome_all, sample == "anyCnoH", bp > 9, bp < 26) %>%
        mutate(bp2 = bp-10)
    pdf("logo_plot_final.pdf", width=4, height=2)
        sample_plot <- steamer_to_plot
        test_data <- data.frame(A=filter(sample_plot, nuc == "A")$dif,
                                    T=filter(sample_plot, nuc == "T")$dif,
                                    G=filter(sample_plot, nuc == "G")$dif,
                                    C=filter(sample_plot, nuc == "C")$dif) %>%
                                    t()
        plot1 <- ggseqlogo(test_data, method = 'custom' ) +
                    ggtitle("anyCnoH") +
                    xlab("Position (Steamer 5bp dup: 7-11)") +
                    ylab("Probability vs expected") +
                    ylim(-0.4,0.4)+
                    theme_classic()
        print(plot1)
    dev.off()


