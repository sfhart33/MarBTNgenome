# Original file here: C:\Users\shart\Metzger Lab Dropbox\Sam_data\Mya_genome\CNV_final\CNV_calling_pipeline.r

# Load packages
    library(cn.mops)
    library(tidyverse)
    #library(vcfR)
    library(mixtools)
    #library(karyoploteR)
    #library(fastseg)
    library(zoo)
    library(bedr)
    library(gridExtra)

# Use cn.mops to count reads
    setwd("/ssd3/Mar_genome_analysis/bwa_mapping/Mar.3.4.6.p1/all_samples")
    BAM_FILES <- list.files(pattern=".bam$")
    CASES <- (BAM_FILES[c(1,2,3,4,5,7,8,9,11,13,14)]) # exclude host contamination samples and tissue samples
    CASES
    # CASES_USA <- (BAM_FILES[c(3,4,6,8)])
    # CASES_PEI <- (BAM_FILES[c(1,2,10)])
    # CONTROL <- (BAM_FILES[c(7,9)])
    # REFERENCE <- (BAM_FILES[5])
    # ##########################################
    bam_CASES <- getReadCountsFromBAM(CASES,
                                        WL = 1000,
                                        parallel = 70)
    # bam_CASES_100 <- getReadCountsFromBAM(CASES,
    #                                     WL = 100,
    #                                     parallel = 70)
    setwd("/ssd3/Mar_genome_analysis/CNV_calling/FINAL")
    saveRDS(bam_CASES, "ReadCountsFromBAM_1kB.rds") # to save
    #saveRDS(bam_CASES_100, "ReadCountsFromBAM_100bp.rds") # to save
    bam_CASES <- readRDS("ReadCountsFromBAM_1kB.rds") # to load
    bam_counts <- as.data.frame(bam_CASES)
    colnames(bam_counts)
    colnames(bam_counts) <- c("chr", "start","end","width","strand","REF","MELCA9","DF490","DF488","DNHL03","DN08","FFM19G1","FFM20B2","FFM22F10","MELCA11","NYTCC9")
    nrow(bam_counts) # 1216076
    head(bam_counts)

# Reference displays two peaks
    # pdf("reference_reads_mapped.pdf")
    #     ggplot(bam_counts, aes(REF)) +
    #         geom_histogram(binwidth = 1) +
    #         xlim(10,1000)+
    #         theme(axis.text=element_text(size=12,face="bold"),
    #             axis.title=element_text(size=16,face="bold"),
    #             text=element_text(size=18,face="bold"),
    #             legend.title = element_blank())
    # dev.off()

# Pre-processing data: normalize by reference, coverage to get logR
    # SPLIT LABELS FROM DATA
        bam_labels <- bam_counts[,1:5]
        bam_counts <- bam_counts[,6:16]

    # ELIMINATE REGIONS WITH LOW COVERAGE IN REFERNCE ANIMAL
        ref_threshold <- mean(bam_counts$REF)/4
        filter(bam_counts,REF > ref_threshold) %>% nrow()/nrow(bam_counts) # 0.9255301
        bam_counts[!(bam_counts$REF > ref_threshold),] <- c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
        bam_counts[bam_counts<1]=1
        head(bam_counts)
        
    # PLOT RAW COUNTS
        pdf("nonreference_reads_mapped.pdf")
        for(i in c("REF","MELCA9","DF490","DF488","DNHL03","DN08","FFM19G1","FFM20B2","FFM22F10","MELCA11","NYTCC9")){
            print(i)
            plot1 <- ggplot(bam_counts, aes(get(i))) +
                geom_histogram(binwidth = 1) +
                xlim(10,1000)+
                xlab(i)+
                theme_classic() +
                theme(axis.text=element_text(size=12,face="bold"),
                    axis.title=element_text(size=16,face="bold"),
                    text=element_text(size=18,face="bold"),
                    legend.title = element_blank())
            print(plot1)
        }
            # ggplot(bam_counts, aes(FFM20B2)) +
            # geom_histogram(binwidth = 1) +
            # xlim(10,1000)
            # ggplot(bam_counts, aes(MELCA11)) +
            # geom_histogram(binwidth = 1) +
            # xlim(10,1000)
            # ggplot(bam_counts, aes(NYTCC9)) +
            # geom_histogram(binwidth = 1) +
            # xlim(10,1000)
            # ggplot(bam_counts, aes(DF488)) +
            # geom_histogram(binwidth = 1) +
            # xlim(10,1000)
            # ggplot(bam_counts, aes(DNHL03)) +
            # geom_histogram(binwidth = 1) +
            # xlim(10,1000)
            # ggplot(bam_counts, aes(DN08)) +
            # geom_histogram(binwidth = 1) +
            # xlim(10,1000)
            # ggplot(bam_counts, aes(REF)) +
            # geom_histogram(binwidth = 1) +
            # xlim(10,1000)
            # ggplot(bam_counts, aes(MELCA9)) +
            # geom_histogram(binwidth = 1) +
            # xlim(10,1000)
            # ggplot(bam_counts, aes(DF490)) +
            # geom_histogram(binwidth = 1) +
            # xlim(10,1000)
        dev.off()

    # NORMALIZE TO REF COVERAGE - NEW
        bam_counts[,2:11] <- bam_counts[,2:11]/bam_counts[,1]
    # NORMALIZE TO AVERAGE COVERAGE AND COVERT TO LOG2
        coverage_all <- as.data.frame((t(t(as.matrix(bam_counts))/colMeans(bam_counts, na.rm=TRUE)))) %>%
            mutate(usaAVG = (FFM19G1+FFM20B2+FFM22F10+MELCA11+NYTCC9)/5,
                   peiAVG = (DF488+DNHL03+DN08)/3)
        coverage_all_logR <- log2(coverage_all)
        head(coverage_all_logR)

    # PLOT LOG2 HISTOGRAM OR FREQPOLY
        pdf("logR_reads_mapped.pdf")
            binsize=0.1
            ggplot(coverage_all_logR) +
                geom_freqpoly(aes(DF488), binwidth = binsize, color = "grey") +
                geom_freqpoly(aes(DNHL03), binwidth = binsize, color = "grey") +
                geom_freqpoly(aes(DN08), binwidth = binsize, color = "grey") +
                geom_freqpoly(aes(peiAVG), binwidth = binsize, color = "black") +
                xlab("LogR")+
                xlim(-3,3)+
                ggtitle("PEI sublineage, logR")+
                theme_classic() +
                theme(axis.text=element_text(size=12,face="bold"),
                    axis.title=element_text(size=16,face="bold"),
                    text=element_text(size=18,face="bold"),
                    legend.title = element_blank())
            ggplot(coverage_all_logR) +
                geom_freqpoly(aes(FFM19G1), binwidth = binsize, color = "grey") +
                geom_freqpoly(aes(FFM20B2), binwidth = binsize, color = "grey") +
                geom_freqpoly(aes(FFM22F10), binwidth =binsize, color = "grey") +
                geom_freqpoly(aes(MELCA11), binwidth =binsize, color = "grey") +
                geom_freqpoly(aes(NYTCC9), binwidth = binsize, color = "grey") +
                geom_freqpoly(aes(usaAVG), binwidth = binsize, color = "black") +
                xlab("LogR")+
                xlim(-3,3)+
                ggtitle("USA sublineage, logR")+
                theme_classic() +
                theme(axis.text=element_text(size=12,face="bold"),
                    axis.title=element_text(size=16,face="bold"),
                    text=element_text(size=18,face="bold"),
                    legend.title = element_blank())
            ggplot(coverage_all_logR) +
                #geom_freqpoly(aes(REF), binwidth = binsize, color = "grey") +
                geom_freqpoly(aes(MELCA9), binwidth = binsize) +
                geom_freqpoly(aes(DF490), binwidth = binsize) +
                xlab("LogR")+
                xlim(-3,3)+
                ggtitle("Healthy clams logR")+
                theme_classic() +
                theme(axis.text=element_text(size=12,face="bold"),
                    axis.title=element_text(size=16,face="bold"),
                    text=element_text(size=18,face="bold"),
                    legend.title = element_blank())
        dev.off()

    # ADD LABELS AND CHANCE SCAF->CHR
        coverage_logR <- cbind(bam_labels,coverage_all_logR) %>%
            separate(col = chr,
                    sep = "scaffold",
                    into = c(NA,"chr"),
                    remove = TRUE,
                    convert = TRUE) %>%
            mutate(chr = chr+1) %>%
            dplyr::select(-width, -strand)
        head(coverage_logR )

    # ROUND TO NEAREST logR OF 0.01 FOR MODE CALLING
        #coverage_logR[, 4:15] <- round(coverage_logR[,4:15], digits = 2)

    #Use read counts to call a peak - mode or median
        # getmode <- function(v) {
        #     uniqv <- unique(v[!is.na(v)])
        #     uniqv[which.max(tabulate(match(v, uniqv)))]
        # }
        # # Call if homozygous deletion, if not call mode
        # cutoff=-2
        # mode_or_deletion <- function(v) {
        #     v<-v[!is.na(v)]
        #     deletion <- v<(cutoff)
        #     if(length(deletion[deletion== FALSE])>length(deletion[deletion== TRUE])){
        #         v <- v[!is.na(v)] # exclude NAs
        #         v <- v[v>(cutoff)] # exclude low values
        #         uniqv <- unique(v) # determine unique values
        #         return(uniqv[which.max(tabulate(match(v, uniqv)))]) # call mode and return
        #         }
        #     if(length(deletion[deletion== TRUE])>=length(deletion[deletion== FALSE])){
        #         return(-10)
        #     }
        # }
        SH_median <- function(input){
            output <- median(input, na.rm=TRUE)
            return(output)
        }
        # TESTING
        # v=c(0, 0.1, 0.1,0.3, 0.2, -3, -3, -3, NA, NA, NA, NA, NA)
        # mode_or_deletion(v)
        # v=c(0, 0.1, -5,-3, -2, -3, NA, NA, NA, NA,] NA)
        # mode_or_deletion(v)
        # 

    # LOOP THROUGH CHROMOSOMES
    pdf("sample_by_chr_coverage.pdf")
    for(i in 1:17){
        coverage_logR_chrX <- filter(coverage_logR, chr == i)
        # Exclude missing regions - call seperately if using mode call
            # coverage_logR_chrX_values <-coverage_logR_chrX[,4:15]
            # coverage_logR_chrX_values[coverage_logR_chrX_values<(-2)]=NA
            # coverage_logR_chrX <- cbind(coverage_logR_chrX[,1:3], coverage_logR_chrX_values)
            # coverage_logR_chrX
            # min(coverage_logR_chrX$MELCA11)
        bins = 100
        # COLLAPSE DOWN TO NEW 100kB REGIONS
            end_pos <- coverage_logR_chrX[nrow(coverage_logR_chrX),3]
            coverage_logR_chrX <- mutate(coverage_logR_chrX, end = lead(end, (bins-1)))
            coverage_logR_chrX_out = coverage_logR_chrX[seq(1, nrow(coverage_logR_chrX), bins), ]
            coverage_logR_chrX_out[nrow(coverage_logR_chrX_out),3] <- end_pos 
        for(column  in 4:16){
            coverage_logR_chrX_out[,column] <- rollapply(coverage_logR_chrX[1:(nrow(coverage_logR_chrX)), column],
                                                    bins,SH_median, by = bins, align = "left", partial = TRUE)
        }
        # Old methods
                                                        #bins,mode_or_deletion, by = bins, align = "left", partial = TRUE)
            # CALCULATE MODES FOR #bins*windowsize bp WINDOWS
            # for(column  in 4:15){
            #   coverage_logR_chrX[,column] <- rollapply(coverage_logR_chrX[1:(nrow(coverage_logR_chrX)), column],
            #                                            bins,getmode, align = "left", partial = TRUE)
            # }
            # for(column  in 4:15){
            #   coverage_logR_chrX_out[,column] <- rollapply(coverage_logR_chrX[seq(1, nrow(coverage_logR_chrX), bins), column],
            #                                            bins,getmode, align = "left", partial = TRUE)
            # }
            
        # print(coverage_logR_chrX)
        # print(coverage_logR_chrX_out)
        # PLOT SUMMARY SAMPLES
            test_plot <- ggplot(coverage_logR_chrX_out)+
                geom_segment(aes(x=start,xend=end,y=peiAVG, yend=peiAVG), size = 2, color = "black")
            print(test_plot)
            test_plot <- ggplot(coverage_logR_chrX_out)+
                geom_segment(aes(x=start,xend=end,y=usaAVG, yend=usaAVG), size = 2, color = "grey")
            print(test_plot)
            test_plot <- ggplot(coverage_logR_chrX_out)+
                geom_segment(aes(x=start,xend=end,y=MELCA9, yend=MELCA9), size = 2, color = "gold")
            print(test_plot)
        # COMBINE CHROMOSOMES INTO ONE OUTPUT
            if(i == 1){
                coverage_logR_1Mb_genomewide <- coverage_logR_chrX_out
            }
            if(i > 1){
                coverage_logR_1Mb_genomewide <- rbind(coverage_logR_1Mb_genomewide, coverage_logR_chrX_out)
            }
            print(paste("Done with:", i))                
    }

    # genome-wide plots
        binsize = 0.05
        ggplot(coverage_logR_1Mb_genomewide) +
            geom_freqpoly(aes(DF488), binwidth = binsize, color = "grey") +
            geom_freqpoly(aes(DNHL03), binwidth = binsize, color = "grey") +
            geom_freqpoly(aes(DN08), binwidth = binsize, color = "grey") +
            geom_freqpoly(aes(peiAVG), binwidth = binsize, color = "black") +
            xlab("LogR")+
            xlim(-3,3)+
            ggtitle("PEI sublineage, 1Mb bins, logR")
        ggplot(coverage_logR_1Mb_genomewide) +
            geom_freqpoly(aes(FFM19G1), binwidth = binsize, color = "grey") +
            geom_freqpoly(aes(FFM20B2), binwidth = binsize, color = "grey") +
            geom_freqpoly(aes(FFM22F10), binwidth = binsize, color = "grey") +
            geom_freqpoly(aes(MELCA11), binwidth = binsize, color = "grey") +
            geom_freqpoly(aes(NYTCC9), binwidth = binsize, color = "grey") +
            geom_freqpoly(aes(usaAVG), binwidth = binsize, color = "black") +
            xlab("LogR")+
            xlim(-3,3)+
            ggtitle("USA sublineage, 1Mb bins, logR")
        ggplot(coverage_logR_1Mb_genomewide) +
            geom_freqpoly(aes(REF), binwidth = binsize, color = "grey") +
            geom_freqpoly(aes(MELCA9), binwidth = binsize) +
            geom_freqpoly(aes(DF490), binwidth = binsize) +
            xlab("LogR")+
            xlim(-3,3)+
            ggtitle("Healthy clams, 1Mb bins logR")
        #coverage_logR_1Mb_genomewide                   

    dev.off()

    # calculate correction for logR -> ploidy
        pdf("logR_to_ploidy_check.pdf")
        binsize = 0.1
        for(peiP in (60:80)/20){
        test_plot <- ggplot(coverage_logR_1Mb_genomewide) +
            geom_freqpoly(aes(peiP*2^DF488), binwidth = binsize, color = "grey") +
            geom_freqpoly(aes(peiP*2^DNHL03), binwidth = binsize, color = "grey") +
            geom_freqpoly(aes(peiP*2^DN08), binwidth = binsize, color = "grey") +
            geom_freqpoly(aes(peiP*2^peiAVG), binwidth = binsize, color = "black") +
            xlim(0,8)+
            ggtitle(paste("PEI sublineage, ploidy:",peiP))
        print(test_plot)
        }
        #peiP = 3.6 # Best fit
        for(usaP in (50:70)/20){
        test_plot <- ggplot(coverage_logR_1Mb_genomewide) +
            geom_freqpoly(aes(usaP*2^FFM19G1), binwidth = binsize, color = "grey") +
            geom_freqpoly(aes(usaP*2^FFM20B2), binwidth = binsize, color = "grey") +
            geom_freqpoly(aes(usaP*2^FFM22F10), binwidth = binsize, color = "grey") +
            geom_freqpoly(aes(usaP*2^MELCA11), binwidth = binsize, color = "grey") +
            geom_freqpoly(aes(usaP*2^NYTCC9), binwidth = binsize, color = "grey") +
            geom_freqpoly(aes(usaP*2^usaAVG), binwidth = binsize, color = "black") +
            xlim(0,8)+
            ggtitle(paste("USA sublineage, ploidy:",usaP))
        print(test_plot)
        }
        #usaP = 3.3 # Best fit
        for(helP in (30:40)/20){
        test_plot <- ggplot(coverage_logR_1Mb_genomewide) +
            geom_freqpoly(aes(helP*2^REF), binwidth = 0.05, color = "grey") +
            xlim(0,4)+
            ggtitle(paste("Healthy  REF clam, ploidy:",helP))
        print(test_plot)
        }
        #refP = 1.6 # Best fit
        for(helP in (30:40)/20){
        test_plot <- ggplot(coverage_logR_1Mb_genomewide) +
            geom_freqpoly(aes(helP*2^MELCA9), binwidth = 0.05) +
            xlim(0,4)+
            ggtitle(paste("Healthy USA clam, ploidy:",helP))
        print(test_plot)
        }
        #helP = 1.9 # Best fit
        for(helP in (30:40)/20){
        test_plot <- ggplot(coverage_logR_1Mb_genomewide) +
            geom_freqpoly(aes(helP*2^DF490), binwidth = 0.05) +
            xlim(0,4)+
            ggtitle(paste("Healthy PEI clam, ploidy:",helP))
        print(test_plot)
        }
        #helP = 1.9 # Best fit
        dev.off()
    # Correction for: MAPQ0, 1000bp bins, refnorm
        peiP=3.6
        usaP=3.3
        HrefP=1.6
        HusaP=1.9
        HpeiP=1.9

# average at end for comparison (makes little difference)
    ploidy <- coverage_logR_1Mb_genomewide %>%
    mutate(usaAVGend = (FFM19G1+FFM20B2+FFM22F10+MELCA11+NYTCC9)/5,
            peiAVGend = (DF488+DNHL03+DN08)/3)
    head(ploidy)
    ploidy[,c(10:15,17)] <- usaP*2^(ploidy[,c(10:15,17)])
    ploidy[,c(7:9,16,18)] <- peiP*2^(ploidy[,c(7:9,16,18)])
    ploidy[,4] <- HrefP*2^(ploidy[,4])
    ploidy[,5] <- HusaP*2^(ploidy[,5])
    ploidy[,6] <- HpeiP*2^(ploidy[,6])
    head(ploidy)
    pdf("avg_start_vs_end.pdf")
        ggplot(ploidy) +
            geom_density_2d_filled(aes(usaAVG,usaAVGend), alpha = 0.5) +
            geom_point(aes(jitter(usaAVG, factor = 1),jitter(usaAVGend, factor = 1)), size = 0.1) +
            theme_classic() +
            theme(axis.text=element_text(size=12,face="bold"),
                    axis.title=element_text(size=16,face="bold"),
                    text=element_text(size=18,face="bold")) +
            ggtitle("USAcancer averaged at end or start comparison")+
            xlim(0,6)+
            ylim(0,6)
        ggplot(ploidy) +
            geom_density_2d_filled(aes(peiAVG,peiAVGend), alpha = 0.5) +
            geom_point(aes(jitter(peiAVG, factor = 1),jitter(peiAVGend, factor = 1)), size = 0.1) +
            theme_classic() +
            theme(axis.text=element_text(size=12,face="bold"),
                    axis.title=element_text(size=16,face="bold"),
                    text=element_text(size=18,face="bold")) +
            ggtitle("PEIcancer averaged at end or start comparison")+
            xlim(0,6)+
            ylim(0,6)
        ggplot(ploidy) +
            geom_density_2d_filled(aes(peiAVG,usaAVG), alpha = 0.5) +
            geom_point(aes(jitter(peiAVG, factor = 1),jitter(usaAVG, factor = 1)), size = 0.1) +
            theme_classic() +
            theme(axis.text=element_text(size=12,face="bold"),
                    axis.title=element_text(size=16,face="bold"),
                    text=element_text(size=18,face="bold")) +
            ggtitle("PEI vs USA averaged start")+
            xlim(0,6)+
            ylim(0,6)
        ggplot(ploidy) +
            geom_density_2d_filled(aes(peiAVGend,usaAVGend), alpha = 0.5) +
            geom_point(aes(jitter(peiAVG, factor = 1),jitter(usaAVG, factor = 1)), size = 0.1) +
            theme_classic() +
            theme(axis.text=element_text(size=12,face="bold"),
                    axis.title=element_text(size=16,face="bold"),
                    text=element_text(size=18,face="bold")) +
            ggtitle("PEI vs USA averaged end")+
            xlim(0,6)+
            ylim(0,6)
    dev.off()

    # smooth based off surrounding calls for ploidy call        
        # SMOOTHING PARAMETERS
            bin=21
            middle=(bin+1)/2
            rightmid=middle-(middle-1)/2
            leftmid=middle+(middle-1)/2
            sd_max=1

        #function to pick whether to smooth right, left or centered
        # ploidy_remove_na <- function(input){
        #   if(is.na(input[2])){
        #     if(is.na(input[1])){
        #       output = 0
        #     } else{
        #       output = input[1]
        #     }
        #   } else{
        #     output = input[2]
        #   }
        #   return(output)
        # }

        ploidy_smoothing <- function(input){
        if(sd(input[rightmid:leftmid]) < sd_max &&
                sd(input[rightmid:leftmid]) < sd(input[1:middle]) &&
                sd(input[rightmid:leftmid]) < sd(input[middle:bin])){
                output = median(input[rightmid:leftmid])
                } else{
                    if(sd(input[1:middle]) < sd_max &&
                    sd(input[1:middle]) < sd(input[middle:bin])){
                    output = median(input[1:middle])
                    } else{
                    if(sd(input[middle:bin]) < sd_max){
                        output = median(input[middle:bin])
                    } else{
                        output = input[middle]
                    }
                    }
                }
        return(output)
        }
        ploidy_smoothing_right <- function(input){
                if(sd(input[1:middle]) < sd_max){
                output = median(input[1:middle])
                } else{
                output = input[middle]
                }
        return(output)
        }
        ploidy_smoothing_left <- function(input){
                if(sd(input[1:middle]) < sd_max){
                output = median(input[1:middle])
                } else{
                output = input[1]
                }
        return(output)
        }
        head(ploidy)
        ploidy2 <- ploidy
        #ploidy2[, 4:18] <- round(ploidy2[, 4:18]) # round before smoothing?
        ploidy2[is.na(ploidy2)] <- 0
        head(ploidy2)
        filter(ploidy2, is.na(REF)) %>% nrow()
        # LOOP THRU CHR
        for(i in 1:17){
            ploidy_chrX <- filter(ploidy2, chr == i)
            ploidy_chrX2 <- ploidy_chrX
            for(column  in 4:18){
                start=(bin+1)/2
                end=(nrow(ploidy_chrX2)-((bin-1)/2))
                # ploidy_chrX2[, column] <- rollapply(ploidy_chrX[, column],
                #                                        1,ploidy_remove_na, align = "center", partial = TRUE)

                ploidy_chrX2[start:end, column] <- rollapply(ploidy_chrX[, column],
                                                        bin,ploidy_smoothing, align = "center")
                ploidy_chrX2[1:(start-1), column] <- rollapply(ploidy_chrX[1:(bin-1), column],
                                                        start,ploidy_smoothing_left, align = "left")
                ploidy_chrX2[(end+1):nrow(ploidy_chrX2), column] <- rollapply(ploidy_chrX[(end-(bin-1)/2+1):nrow(ploidy_chrX2), column],
                                                        start,ploidy_smoothing_right, align = "right")
            }
        # COMBINE CHROMOSOMES INTO ONE OUTPUT
            if(i == 1){
                ploidy_genome <- ploidy_chrX2
            }
            if(i > 1){
                ploidy_genome <- rbind(ploidy_genome, ploidy_chrX2)
            }
        }
        head(ploidy_genome, n=100)
        ploidy_genome_unrounded <- ploidy_genome
        ploidy_genome[, 4:18] <- round(ploidy_genome[, 4:18])

# FINAL GENOME-WIDE PLOT
        genome_full <- read.delim("/ssd3/Mar_genome_analysis/genomes/Mar.3.4.6.p1_Q30Q30A.fasta.fai",header=FALSE, col.names=c("scaf","end","A","B","C"))
        genome_full$chr <- 1:17
    genome_plot <- function(chr_col,full_ploidy,ploidy_call,maintitle,y_limit){
        colourTotal="gray87"
        colourRound="black"
        chr.segs=genome_full$end
        chr.names=genome_full$chr
        par(mar = c(0.5,5,5,0.5), cex = 0.4, cex.main=3, cex.axis = 2.5)
        ticks=seq(0, y_limit, 1)
        plot(c(1,length(full_ploidy)), c(0,y_limit+1), type = "n", xaxt = "n", yaxt="n", main = maintitle, xlab = "", ylab = "", xaxs="i")
        axis(side = 2, at = ticks)
        abline(h=ticks, col="lightgrey", lty=1)
        # ploidy plot
        # plot sample calls
        A_rle<-rle(full_ploidy)
        start=0
        for(i in 1:length(A_rle$values)){
            val<-A_rle$values[i]
            size<-A_rle$lengths[i]
            rect(start, (val-0.05), (start+size-1), (val+0.05), col=colourTotal, border=colourTotal)
            start=start+size
        }
        A_rle<-rle(ploidy_call)
        start=0
        # plot total copy number rounded
        for(i in 1:length(A_rle$values)){
            val<-A_rle$values[i]
            size<-A_rle$lengths[i]
            rect(start, (val-0.05), (start+size-1), (val+0.05), col=colourRound, border=NA) #colourRound
            start=start+size
        }
        # chr names
        chrk_tot_len = 0
        abline(v=0,lty=1,col="lightgrey")
        for (i in chr.names){
            #chrk = chr.segs[i];
            chrk_tot_len_prev = chrk_tot_len
            chrk_tot_len = chrk_tot_len + sum(chr_col==i)
            vpos = chrk_tot_len;
            tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
            text(tpos,y_limit+1,i, pos = 1, cex = 4)
            abline(v=vpos,lty=1,col="lightgrey")
        }
    }

    pdf("CN_genomewide.pdf", width=24, height=6)
    #genome_plot(ploidy$chr,ploidy$peiAVGend,ploidy_genome$peiAVGend,"PEI (end avg)",7)
    genome_plot(ploidy$chr,ploidy$peiAVG,ploidy_genome$peiAVG,"PEI (start avg)",7)
    #genome_plot(ploidy$chr,ploidy$usaAVGend,ploidy_genome$usaAVGend,"USA (end avg)",7)
    genome_plot(ploidy$chr,ploidy$usaAVG,ploidy_genome$usaAVG,"USA (start avg)",7)
    genome_plot(ploidy$chr,ploidy$MELCA9,ploidy_genome$MELCA9,"MELCA9",7)
    genome_plot(ploidy$chr,ploidy$DF490,ploidy_genome$DF490,"DF490",7)
    #genome_plot(ploidy$chr,ploidy$REF,ploidy_genome$REF,"REF",3)
    dev.off()

    pdf("CN_freqpolys.pdf", width=6, height=6)
        binsize = 0.1
        ggplot(ploidy) +
            geom_freqpoly(aes(DF488), binwidth = binsize, color = "grey") +
            geom_freqpoly(aes(DNHL03), binwidth = binsize, color = "grey") +
            geom_freqpoly(aes(DN08), binwidth = binsize, color = "grey") +
            geom_freqpoly(aes(peiAVG), binwidth = binsize, color = "black") +
            xlab("LogR")+
            #xlim(0,8)+
            scale_x_continuous(expand = c(0, 0), limits =c(0,8)) +
            scale_y_continuous(expand = c(0, 0)) +
            theme_classic() +
            theme(axis.text=element_text(size=12,face="bold"),
                axis.title=element_text(size=16,face="bold"),
                text=element_text(size=14,face="bold")) +
            coord_flip() +
            ggtitle("PEI sublineage")
        ggplot(ploidy) +
            geom_freqpoly(aes(FFM19G1), binwidth = binsize, color = "grey") +
            geom_freqpoly(aes(FFM20B2), binwidth = binsize, color = "grey") +
            geom_freqpoly(aes(FFM22F10), binwidth = binsize, color = "grey") +
            geom_freqpoly(aes(MELCA11), binwidth = binsize, color = "grey") +
            geom_freqpoly(aes(NYTCC9), binwidth = binsize, color = "grey") +
            geom_freqpoly(aes(usaAVG), binwidth = binsize, color = "black") +
            xlab("LogR")+
            #xlim(0,8)+
            scale_x_continuous(expand = c(0, 0), limits =c(0,8)) +
            scale_y_continuous(expand = c(0, 0)) +
            theme_classic() +
            theme(axis.text=element_text(size=12,face="bold"),
                axis.title=element_text(size=16,face="bold"),
                text=element_text(size=14,face="bold")) +
            coord_flip() +
            ggtitle("USA sublineage")
        ggplot(ploidy) +
            #geom_freqpoly(aes(REF), binwidth = binsize, color = "grey") +
            geom_freqpoly(aes(MELCA9), binwidth = binsize) +
            geom_freqpoly(aes(DF490), binwidth = binsize) +
            xlab("LogR")+
            #xlim(0,3)+
            scale_x_continuous(expand = c(0, 0), limits =c(0,8)) +
            scale_y_continuous(expand = c(0, 0)) +
            theme_classic() +
            theme(axis.text=element_text(size=12,face="bold"),
                axis.title=element_text(size=16,face="bold"),
                text=element_text(size=14,face="bold")) +
            coord_flip() +
            ggtitle("Healthy clams")   
    dev.off()
    # ploidy_noH0 <- filter(ploidy, MELCA9>1.5 & DF490>1.5)
    # pdf("CN_freqpolys_noH0-1.pdf", width=6, height=6)
    #     binsize = 0.1
    #     ggplot(ploidy_noH0) +
    #         geom_freqpoly(aes(DF488), binwidth = binsize, color = "grey") +
    #         geom_freqpoly(aes(DNHL03), binwidth = binsize, color = "grey") +
    #         geom_freqpoly(aes(DN08), binwidth = binsize, color = "grey") +
    #         geom_freqpoly(aes(peiAVG), binwidth = binsize, color = "black") +
    #         xlab("LogR")+
    #         scale_x_continuous(expand = c(0, 0), limits =c(0,8)) +
    #         scale_y_continuous(expand = c(0, 0)) +
    #         theme_classic() +
    #         theme(axis.text=element_text(size=12,face="bold"),
    #             axis.title=element_text(size=16,face="bold"),
    #             text=element_text(size=14,face="bold")) +
    #         coord_flip() +
    #         ggtitle(paste("PEI sublineage, 1Mb mode bins, set logR(0)=:",peiP))
    #     ggplot(ploidy_noH0) +
    #         geom_freqpoly(aes(FFM19G1), binwidth = binsize, color = "grey") +
    #         geom_freqpoly(aes(FFM20B2), binwidth = binsize, color = "grey") +
    #         geom_freqpoly(aes(MELCA11), binwidth = binsize, color = "grey") +
    #         geom_freqpoly(aes(NYTCC9), binwidth = binsize, color = "grey") +
    #         geom_freqpoly(aes(usaAVG), binwidth = binsize, color = "black") +
    #         xlab("LogR")+
    #         scale_x_continuous(expand = c(0, 0), limits =c(0,8)) +
    #         scale_y_continuous(expand = c(0, 0)) +
    #         theme_classic() +
    #         theme(axis.text=element_text(size=12,face="bold"),
    #             axis.title=element_text(size=16,face="bold"),
    #             text=element_text(size=14,face="bold")) +
    #         coord_flip() +
    #         ggtitle(paste("USA sublineage, 1Mb mode bins, set logR(0)=:",usaP))
    #     ggplot(ploidy_noH0) +
    #         #geom_freqpoly(aes(REF), binwidth = binsize, color = "grey") +
    #         geom_freqpoly(aes(MELCA9), binwidth = binsize) +
    #         geom_freqpoly(aes(DF490), binwidth = binsize) +
    #         xlab("LogR")+
    #         scale_x_continuous(expand = c(0, 0), limits =c(0,4)) +
    #         scale_y_continuous(expand = c(0, 0)) +
    #         theme_classic() +
    #         theme(axis.text=element_text(size=12,face="bold"),
    #             axis.title=element_text(size=16,face="bold"),
    #             text=element_text(size=14,face="bold")) +
    #         coord_flip() +
    #         ggtitle(paste("Healthy clams, 1Mb mode bins, set logR(0)=:","X"))   
    # dev.off()

# comparison plots
    for(A in c("MELCA9","DF490","DF488","DNHL03","DN08","FFM19G1","FFM20B2","FFM22F10","MELCA11","NYTCC9")){
        print(A)
        for(B in c("MELCA9","DF490","DF488","DNHL03","DN08","FFM19G1","FFM20B2","FFM22F10","MELCA11","NYTCC9")){
            if(A != B){
            print(paste0("   vs   ", B))
            ploidy_regression <- lm(get(A) ~ get(B), data=ploidy)
            print(paste0("   ", summary(ploidy_regression)$r.squared)) 
            }

        }  
    }
    #PEI vs PEI all ~0.98
    #PEI vs USA all ~0.52-0.56
    #PEI vs USA avg 0.563
    #USA vs USA all ~0.94-0.98

    # (Too big a file)
    # pdf("densityplot_comparisons_fullgenome.pdf") #, width=4, height=4)
    # for(A in c("MELCA9","DF490","DF488","DNHL03","DN08","FFM19G1","FFM20B2","FFM22F10","MELCA11","NYTCC9")){
    #     for(B in c("MELCA9","DF490","DF488","DNHL03","DN08","FFM19G1","FFM20B2","FFM22F10","MELCA11","NYTCC9")){
    #         if(A != B){
    #             ploidy_regression <- lm(get(A) ~ get(B), data=ploidy)
    #             plot1 <- ggplot(ploidy, aes(get(A),get(B))) +
    #                 #geom_density_2d_filled(aes(,), alpha = 0.5) +
    #                 geom_point(, size = 0.1, alpha = 0.1) +
    #                 theme_classic() +
    #                 theme(axis.text=element_text(size=12,face="bold"),
    #                         axis.title=element_text(size=16,face="bold"),
    #                         text=element_text(size=18,face="bold")) +
    #                 ggtitle(paste(A,"vs",B,"r.squared:",round(summary(ploidy_regression)$r.squared,3)))+
    #                 xlim(0,6)+
    #                 ylim(0,6)+
    #                 xlab(A)+
    #                 ylab(B) 
    #             print(plot1)
    #         }
    #     }  
    # }
    # dev.off()

    plot_CN_comparison <- function(A,B){
        ploidy_regression <- lm(get(A) ~ get(B), data=ploidy)
        plot1 <- ggplot(ploidy, aes(get(A),get(B))) +
            #geom_density_2d_filled(aes(,), alpha = 0.5) +
            geom_point(, size = 0.5, alpha = 0.5) +
            theme_classic() +
            theme(axis.text=element_text(size=12,face="bold"),
                    axis.title=element_text(size=16,face="bold"),
                    text=element_text(size=18,face="bold")) +
            ggtitle(paste(A,"vs",B,"r.squared:",round(summary(ploidy_regression)$r.squared,3)))+
            xlim(0,6)+
            ylim(0,6)+
            xlab(A)+
            ylab(B) 
        print(plot1)
    }
    pdf("densityplot_comparisons_fullgenome.pdf") #, width=4, height=4)
    plot_CN_comparison("usaAVG", "peiAVG")
    plot_CN_comparison("DF488", "DN08")
    plot_CN_comparison("FFM19G1", "NYTCC9")
    dev.off()

# Pie chart
    # piechart_ploidy <- data.frame(ploidy_genome$peiAVG,
    #                             ploidy_genome$usaAVG,
    #                             ploidy_genome$MELCA9,
    #                             ploidy_genome$DF490,
    #                             ploidy_genome$REF
    #                             ) %>%
    #                             filter(H_USA == 0 & H_PEI == 0)
    # colnames(piechart_ploidy) <- c("PEI", "USA", "H_USA", "H_PEI", "H_REF")
    # piechart_ploidy[piechart_ploidy >= 7]  = "7+"

    # pdf("CNV_piecharts.pdf", width=4, height=4)
    #     ggplot(piechart_ploidy, aes(x="", y="", fill=PEI)) +
    #         geom_bar(stat="identity", width=1) +
    #         coord_polar("y", start=0) +
    #         theme_void()+
    #         scale_fill_brewer(palette="Set1")
    #     ggplot(piechart_ploidy, aes(x="", y="", fill=USA)) +
    #         geom_bar(stat="identity", width=1) +
    #         coord_polar("y", start=0)+
    #         theme_void()+
    #         scale_fill_brewer(palette="Set1")
    #     ggplot(piechart_ploidy, aes(x="", y="", fill=H_PEI)) +
    #         geom_bar(stat="identity", width=1) +
    #         coord_polar("y", start=0)+
    #         theme_void()+
    #         scale_fill_brewer(palette="Set1")
    #     ggplot(piechart_ploidy, aes(x="", y="", fill=H_USA)) +
    #         geom_bar(stat="identity", width=1) +
    #         coord_polar("y", start=0)+
    #         theme_void()+
    #         scale_fill_brewer(palette="Set1")
    #     ggplot(piechart_ploidy, aes(x="", y="", fill=H_REF)) +
    #         geom_bar(stat="identity", width=1) +
    #         coord_polar("y", start=0)+
    #         theme_void()+
    #         scale_fill_brewer(palette="Set1")
    # dev.off()
            
# bed output for healthy freq check
    head(ploidy_genome)
    ploidy_round <- ploidy_genome
    ploidy_round[,1] <- paste("Mar.3.4.6.p1_scaffold", (ploidy_round[,1]-1), sep = "")
    # ploidy_round_noH0 <- filter(ploidy_round, MELCA9>1.5 | DF490>1.5 )
    head(ploidy_round)
    CN_to_bed <- function(input, column, max, name){
        for(i in 0:max){
            CN_subset <- filter(input, get(column) == i) %>%
                mutate(start = start-1, end = end-1) %>%
                dplyr::select(chr,start,end)
            write.table(CN_subset, file = paste(name,"_CN",i,".bed",sep=""),
                        sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
        }
        CN_subset <- filter(input, get(column) > max) %>%
                mutate(start = start-1, end = end-1) %>%
                dplyr::select(chr,start,end)
            write.table(CN_subset, file = paste(name,"_CN",(max+1),"plus.bed",sep=""),
                        sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    }
    setwd("/ssd3/Mar_genome_analysis/CNV_calling/FINAL/output_bed")
    CN_to_bed(ploidy_round,"usaAVGend",7,"USA")
    #CN_to_bed(ploidy_round_noH0,"usaAVGend",6,"USA_noH0-1")
    CN_to_bed(ploidy_round,"peiAVGend",7,"PEI")
    #CN_to_bed(ploidy_round_noH0,"peiAVGend",6,"PEI_noH0-1")
