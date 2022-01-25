# Original files here: C:\Users\shart\Metzger Lab Dropbox\Sam_data\Mya_genome\LOH\somatypus_output_LOH_testing.r

library(sigfit)
library(tidyverse)
library(lsa)
data(cosmic_signatures_v2)
data(cosmic_signatures_v3)

# Load helmsman data and transform to sigfit input
    helmsman_output_file <- "/ssd3/Mar_genome_analysis/LOH/july_2021/output/new/helmsman/subtype_count_matrix.txt"
    helmsman <- read.table(helmsman_output_file,
                                    sep="\t",
                    skip = 1)
    rownames(helmsman) <- helmsman$V1
    helmsman <- select(helmsman, -V1)
    colnames(helmsman) <- colnames(cosmic_signatures_v2)
    helmsman <- helmsman[order(row.names(helmsman)), ]
    head(helmsman)

# Load mutational oppertunities and extracted signatures
    mutational_oppertunities_many <- readRDS("/ssd3/Mar_genome_analysis/mut_sig/trinuc_freq/mutational_oppertunities_all.rds")
    mutopps <- mutational_oppertunities_many[rep(5, nrow(helmsman)),]
    # Latest full extraction finalized
        sig_extraction <- readRDS(file = "/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/sigfit/signatures_4_extracted_full.rds")
    # Prelim extraction
        # sig_extraction <- readRDS(file = "/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/sigfit/signatures_4_prelim_extracted.rds")
    # Previous runs done with this extraction:
        # sig_extraction <- readRDS(file = "/ssd3/Mar_genome_analysis/mut_sig/sigfit/final_bins/corrected_mutopps/4sigs_final_extraction_genome-gene-cds_mutoppcorr.rds")
        # signatures4 <- retrieve_pars(sig_extraction, par = "signatures")
        signatures4 <- sig_extraction

# Fit to signatures
    setwd("/ssd3/Mar_genome_analysis/LOH/july_2021/output/new")
    helmsman_fit <- fit_signatures(counts = helmsman,  
                                    signatures = signatures4,
                                    opportunities =  mutopps,
                                    iter = 2000, 
                                    warmup = 1000, 
                                    chains = 1, 
                                    seed = 1756)
    saveRDS(helmsman_fit, "LOH_cutoffs_sigfit.rds")
    #helmsman_fit <- readRDS("LOH_cutoffs_sigfit.rds")
    #plot_exposures(mcmc_samples = helmsman_fit , pdf_path = "helmsman_fit_exposures.pdf")
    #plot_reconstruction(mcmc_samples = helmsman_fit , pdf_path = "helmsman_fit_reconstruction.pdf")

# View SigS exposures
    output_data <- retrieve_pars(helmsman_fit, par = "exposures")
    output_data <- output_data$mean 
    #output_data <- select(output_data, "sigS")
    #output_data <- select(output_data, "Signature A")
    output_data
    output_data$name = rownames(output_data)
    rownames(output_data) <- NULL
    output_data

#Transform to plot and determine ideal # of SNVs to cutoff at
    output_data <- output_data %>%
        separate(name, c("LOH",NA,"count",NA,"rest"), sep = "_") %>%
        #separate(name, c("LOH",NA,"count","rest"), sep = "_") %>%
        separate(rest, c(NA,"SNVs"), sep = ".merge.")
    output_data$count <- as.integer(output_data$count)
    output_data <- arrange(output_data, count)
    output_data

# Correction
    null_data <- filter(output_data, count==51)
    null_pei_sigS <- null_data[1,4]
    null_usa_sigS <- null_data[2,4]
    null_pei_sig40 <- null_data[1,3]
    null_usa_sig40 <- null_data[1,3]
    output_data <- filter(output_data, count<51)


    USAloh <- filter(output_data, SNVs == "USAloh") %>%
        select(count, sig40, sigS) %>%
        rename(USAloh_sig40 = sig40, USAloh_sigS = sigS)
    USAsom <- filter(output_data, SNVs == "USAsom") %>%
        select(sig40, sigS) %>%
        rename(USAsom_sig40 = sig40, USAsom_sigS = sigS)
    #USAsom[1,] <- c(NA,NA) #c(0,0)
    PEIloh <- filter(output_data, SNVs == "PEIloh") %>%
        select(sig40, sigS) %>%
        rename(PEIloh_sig40 = sig40, PEIloh_sigS = sigS)
    PEIsom <- filter(output_data, SNVs == "PEIsom") %>%
        select(sig40, sigS) %>%
        rename(PEIsom_sig40 = sig40, PEIsom_sigS = sigS)
    #PEIsom[1, ] <- c(NA,NA) #c(0,0)
    full_data <- cbind(USAloh,USAsom,PEIloh,PEIsom)
    full_data <- mutate(full_data,
                        PEIdif_sigS = PEIsom_sigS - PEIloh_sigS,
                        USAdif_sigS = USAsom_sigS - USAloh_sigS, 
                        PEIdif_sig40 =  PEIloh_sig40 - PEIsom_sig40,
                        USAdif_sig40 = USAloh_sig40 - USAsom_sig40)
    # null_data <- filter(full_data, count==0)
    # null_pei_sigS <- null_data[1,7]
    # null_usa_sigS <- null_data[1,3]
    # null_pei_sig40 <- null_data[1,6]
    # null_usa_sig40 <- null_data[1,2]
    full_data <- filter(full_data, count>0)


# Load data about size of LOH regions
    LOH_size <- read.table("LOH.count.test.txt", header=T,fill=TRUE) %>%
        separate(Sample, c("region",NA,"count",NA), sep = "_") %>%
        arrange(as.integer(count))
    LOH_size_USA <- filter(LOH_size, region=="USA")
    LOH_size_PEI <- filter(LOH_size, region=="PEI")

    output_data2 <- filter(output_data, count >0)

# pdf("optimal_cutoff_plotb.pdf")
# ggplot(output_data2)+
#     geom_line(aes(x=count,y=sigS,color=factor(SNVs)), size=1, linetype="solid")+ #, color="red"
#     scale_color_manual(values=c("red","red2","blue","blue2"))+
#     geom_abline(aes(intercept=null_pei_sigS, slope=0), color="red", linetype="dotted")+
#     geom_abline(aes(intercept=null_usa_sigS, slope=0), color="blue", linetype="dotted")+
#     xlab("# discordant homozygous alleles in other sublineage")+
#     ylab("Fraction of SigS exposure")+
#     theme_classic() +
#         theme(axis.text=element_text(size=12,face="bold"),
#         axis.title=element_text(size=16,face="bold"),
#         text=element_text(size=16,face="bold")) +
#     ggtitle("PEI_SNVs in USA_LOH (red)
# USA_SNVs in PEI_LOH (blue)")
# dev.off()

pdf("optimal_cutoff_plot_hetero.pdf")
ggplot(full_data)+
    # geom_line(aes(count,PEIsom_sigS), size=1, color="red", linetype="solid")+
    # geom_line(aes(count,USAsom_sigS), size=1, color="blue", linetype="solid")+
    # geom_line(aes(count,PEIloh_sigS), size=1, color="red", linetype="dashed")+
    # geom_line(aes(count,USAloh_sigS), size=1, color="blue", linetype="dashed")+
    geom_point(aes(count,PEIsom_sigS), size=2, fill="red", shape = 21)+
    geom_point(aes(count,USAsom_sigS), size=2, fill="blue", shape = 21)+
    geom_point(aes(count,PEIloh_sigS), size=2, fill="red", shape = 22)+
    geom_point(aes(count,USAloh_sigS), size=2, fill="blue", shape = 22)+
    geom_abline(aes(intercept=null_pei_sigS, slope=0), color="red", linetype="dashed")+
    geom_abline(aes(intercept=null_usa_sigS, slope=0), color="blue", linetype="dashed")+
    geom_vline(xintercept=10, color="black", linetype="dashed")+
    xlab("# discordant homozygous alleles to call LOH")+
    ylab("Fraction of SNVs attributed to SigS")+
    theme_classic() +
        theme(axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=16,face="bold"),
        text=element_text(size=16,face="bold")) #+
    #ggtitle("PEI_SNVs in USA_LOH (red)\nUSA_SNVs in PEI_LOH (blue)")
ggplot(full_data)+
    # geom_line(aes(count,PEIdif_sigS), size=1, color="red", linetype="solid")+
    # geom_line(aes(count,USAdif_sigS), size=1, color="blue", linetype="solid")+
    geom_point(aes(count,PEIdif_sigS), size=2, fill="red", shape = 21)+
    geom_point(aes(count,USAdif_sigS), size=2, fill="blue", shape = 21)+
    geom_line(aes(count,(USAdif_sigS+PEIdif_sigS)/2), size=1, color="black", linetype="solid")+
    geom_vline(xintercept=10, color="black", linetype="dashed")+
    xlab("# discordant homozygous alleles to call LOH")+
    ylab("Difference in SigS fraction: LOH vs nonLOH")+
    theme_classic() +
        theme(axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=16,face="bold"),
        text=element_text(size=16,face="bold")) #+
    #ggtitle("PEI_SNVs in USA_LOH (red)\nUSA_SNVs in PEI_LOH (blue)\nAverage(black)")
# ggplot(full_data)+
#     geom_point(aes(PEIsom_sigS,PEIsom_sig40), size=2, color="red")+
#     geom_point(aes(USAsom_sigS,USAsom_sig40), size=2, color="blue")+
#     xlab("sigS (each point is one cutoff value)")+
#     ylab("sig40")+
#     theme_classic() +
#         theme(axis.text=element_text(size=12,face="bold"),
#         axis.title=element_text(size=16,face="bold"),
#         text=element_text(size=16,face="bold")) +
#     ggtitle("SigS negatively coorelates with sig40")
# ggplot(full_data)+
#     geom_line(aes(count,PEIsom_sig40), size=1, color="red", linetype="solid")+
#     geom_line(aes(count,USAsom_sig40), size=1, color="blue", linetype="solid")+
#     geom_line(aes(count,PEIloh_sig40), size=1, color="red", linetype="dashed")+
#     geom_line(aes(count,USAloh_sig40), size=1, color="blue", linetype="dashed")+
#     geom_abline(aes(intercept=null_pei_sig40, slope=0), color="red", linetype="dotted")+
#     geom_abline(aes(intercept=null_usa_sig40, slope=0), color="blue", linetype="dotted")+
#     xlab("# discordant homozygous alleles in other sublineage")+
#     ylab("Fraction of sig40 exposure")+
#     theme_classic() +
#         theme(axis.text=element_text(size=12,face="bold"),
#         axis.title=element_text(size=16,face="bold"),
#         text=element_text(size=16,face="bold")) +
#     ggtitle("PEI_SNVs in USA_LOH (red)
# USA_SNVs in PEI_LOH (blue)")
# ggplot(full_data)+
#     geom_line(aes(count,PEIdif_sig40), size=1, color="red", linetype="solid")+
#     geom_line(aes(count,USAdif_sig40), size=1, color="blue", linetype="solid")+
#     geom_line(aes(count,(USAdif_sig40+PEIdif_sig40)/2), size=1, color="black", linetype="dotted")+
#     xlab("# discordant homozygous alleles in other sublineage")+
#     ylab("Difference in sig40 exposure: LOH vs nonLOH")+
#     theme_classic() +
#         theme(axis.text=element_text(size=12,face="bold"),
#         axis.title=element_text(size=16,face="bold"),
#         text=element_text(size=16,face="bold")) +
#     ggtitle("PEI_SNVs in USA_LOH (red)
# USA_SNVs in PEI_LOH (blue)
# Average(black)")
ggplot()+
    geom_point(data=LOH_size_PEI, aes(as.integer(count),bp/1215282629), size=2, fill="red", shape=21)+
    geom_point(data=LOH_size_USA, aes(as.integer(count),bp/1215282629), size=2, fill="blue", shape=21)+
    geom_vline(xintercept=10, color="black", linetype="dashed")+
    geom_hline(yintercept=0, color="black", linetype="solid")+
    geom_hline(yintercept=LOH_size_USA[10,4]/1215282629, color="blue", linetype="dashed")+
    geom_hline(yintercept=LOH_size_PEI[10,4]/1215282629, color="red", linetype="dashed")+
    xlab("# discordant homozygous alleles to call LOH")+
    ylab("Fraction of genome in LOH regions")+
    theme_classic() +
        theme(axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=16,face="bold"),
        text=element_text(size=16,face="bold")) #+
    #ggtitle("PEI_LOH (red), USA_LOH (blue)")
# ggplot()+
#     geom_line(data=LOH_size_PEI, aes(as.integer(count),bp/1215282629), size=1, color="red", linetype="solid")+
#     geom_line(data=LOH_size_USA, aes(as.integer(count),bp/1215282629), size=1, color="blue", linetype="solid")+
#     xlab("# discordant homozygous alleles in other sublineage")+
#     ylab("Fraction of genome in LOH regions")+
#     xlim(0,20)+
#     ylim(0,0.2)+
#     theme_classic() +
#         theme(axis.text=element_text(size=12,face="bold"),
#         axis.title=element_text(size=16,face="bold"),
#         text=element_text(size=16,face="bold")) +
#     ggtitle("PEI_LOH (red), USA_LOH (blue)")
dev.off()




######################################## 8/12/21

setwd("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/bins")

# Load abberent counts data (allPEInoUSAanyH & allUSAnoPEIanyH)
    abberent_counts <- read.table("LOH_abberent_SNVs_counts.txt", header=T)
    genome <- read.table("/ssd3/Mar_genome_analysis/genomes/Mar.3.4.6.p1_Q30Q30A.fasta.fai", header=F)
    genomesize <- sum(genome$V2)
    abberent_counts <- abberent_counts %>% 
    mutate(LOH_FRAC = LOH_SIZE/genomesize,
           nonLOH_SIZE = genomesize - LOH_SIZE
           ) %>%
    mutate(LOH_SNVs_freq =  LOH_SNVs/LOH_SIZE,
           nonLOH_SNVs_freq = nonLOH_SNVs/nonLOH_SIZE
           )

# count baseline values   
    allPEInoUSAanyH_regular_count <- read.table("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/bins/allPEInoUSAanyH_regular.bed", header=F) %>% nrow()/genomesize
    allPEInoUSAanyH_stringent_count <- read.table("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/bins/allPEInoUSAanyH_stringent.bed", header=F) %>% nrow()/genomesize
    allUSAnoPEIanyH_regular_count <- read.table("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/bins/allUSAnoPEIanyH_regular.bed", header=F) %>% nrow()/genomesize
    allUSAnoPEIanyH_stringent_count <- read.table("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/bins/allUSAnoPEIanyH_stringent.bed", header=F) %>% nrow()/genomesize

pdf("LOH_abberent_SNVs.pdf")
ggplot(abberent_counts)+
    geom_line(aes(x=count,y=LOH_SNVs_freq,color=factor(paste(snv_region,stringency))), size=1,linetype="dotted")+ #, color="red" linetype="solid", linetype=factor(stringency)
    geom_line(aes(x=count,y=nonLOH_SNVs_freq,color=factor(paste(snv_region,stringency))), size=1,linetype="solid")+
    scale_color_manual("SNVs, Stringeny",values=c("red","red2","blue","blue2"))+
    geom_abline(aes(intercept=allPEInoUSAanyH_regular_count, slope=0), color="red", linetype="dashed")+
    geom_abline(aes(intercept=allPEInoUSAanyH_stringent_count, slope=0), color="red2", linetype="dashed")+
    geom_abline(aes(intercept=allUSAnoPEIanyH_regular_count, slope=0), color="blue", linetype="dashed")+
    geom_abline(aes(intercept=allUSAnoPEIanyH_stringent_count, slope=0), color="blue2", linetype="dashed")+
    xlab("# discordant homozygous alleles in other sublineage/50")+
    ylab("Frequency of abberent SNVs")+
    ylim(0,NA) +
    theme_classic() +
        theme(axis.text=element_text(size=10,face="bold"),
        axis.title=element_text(size=10,face="bold"),
        text=element_text(size=10,face="bold")) +
    ggtitle("Frequency of abberent SNVs inLOH (dotted) and nonLOH (solid)")
ggplot(abberent_counts)+
    geom_line(aes(x=count,y=LOH_SNVs_freq,color=factor(paste(snv_region,stringency))), size=1,linetype="solid")+
    scale_color_manual("SNVs, Stringeny",values=c("red","red2","blue","blue2"))+
    geom_abline(aes(intercept=allPEInoUSAanyH_regular_count, slope=0), color="red", linetype="dashed")+
    geom_abline(aes(intercept=allPEInoUSAanyH_stringent_count, slope=0), color="red2", linetype="dashed")+
    geom_abline(aes(intercept=allUSAnoPEIanyH_regular_count, slope=0), color="blue", linetype="dashed")+
    geom_abline(aes(intercept=allUSAnoPEIanyH_stringent_count, slope=0), color="blue2", linetype="dashed")+
    xlab("# discordant homozygous alleles in other sublineage/50")+
    ylab("Frequency of abberent SNVs")+
    ylim(0,NA) +
    theme_classic() +
        theme(axis.text=element_text(size=10,face="bold"),
        axis.title=element_text(size=10,face="bold"),
        text=element_text(size=10,face="bold")) +
    ggtitle("Frequency of abberent SNVs in LOH")
ggplot(abberent_counts)+
    geom_line(aes(x=count,y=nonLOH_SNVs_freq,color=factor(paste(snv_region,stringency))), size=1,linetype="solid")+
    scale_color_manual("SNVs, Stringeny",values=c("red","red2","blue","blue2"))+
    geom_abline(aes(intercept=allPEInoUSAanyH_regular_count, slope=0), color="red", linetype="dashed")+
    geom_abline(aes(intercept=allPEInoUSAanyH_stringent_count, slope=0), color="red2", linetype="dashed")+
    geom_abline(aes(intercept=allUSAnoPEIanyH_regular_count, slope=0), color="blue", linetype="dashed")+
    geom_abline(aes(intercept=allUSAnoPEIanyH_stringent_count, slope=0), color="blue2", linetype="dashed")+
    xlab("# discordant homozygous alleles in other sublineage/50")+
    ylab("Frequency of abberent SNVs")+
    ylim(0,NA) +
    theme_classic() +
        theme(axis.text=element_text(size=10,face="bold"),
        axis.title=element_text(size=10,face="bold"),
        text=element_text(size=10,face="bold")) +
    ggtitle("Frequency of abberent SNVs in nonLOH")
dev.off()


######################################################################################################
# after creating dnds files for all and running dndscv
# /ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/bins/dnds/README


library(tidyverse)
setwd("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/bins/dnds/LOH")
dnds_files <- list.files(pattern = ".dnds.rds")

testing <- readRDS("allPEInoUSAnoH_regular_37_LOH.dnds.rds")
testing$globaldnds

count=0
for(dnds_file in dnds_files){
    dnds_name <- str_split(dnds_file,".dnds")[[1]][1]
    print(paste("Loading",dnds_name), quote = FALSE)
    dnds_output <- readRDS(dnds_file)
    print(dnds_output$globaldnds[1,])
    if(count == 0){
        dnds_summary <- data.frame(dnds_output$globaldnds[1,], stringsAsFactors = FALSE) %>%
		mutate(name = dnds_name)
        count=count+1
    }
    if(count > 0){
        dnds_value <- data.frame(dnds_output$globaldnds[1,], stringsAsFactors = FALSE) %>%
            mutate(name = dnds_name)
        dnds_summary <- rbind(dnds_summary,dnds_value)
    }
    print("     done", quote = FALSE)	
}

print(dnds_summary)
saveRDS(dnds_summary, "dnds_global_summary.rds")

names <- c("name", "mle", "cilow", "cihigh")
testing <- data.frame(name="allPEInoUSAnoH_regular_1_LOH", "mle"=0.8313108, "cilow"=0.7766297, "cihigh"=0.8898419)

dnds_summary_new <- dnds_summary %>%
        separate(name, c("SNVs","stringency","count","LOH"), sep = "_") %>%
        separate(SNVs, c("SNVs",NA,NA), sep = "no") %>%
        separate(SNVs, c(NA,"SNVs"), sep = "ll")

# load threshold values
    normal_dnds <- readRDS("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/bins/dnds/outputs/dnds_global_summary.rds")
    pei_starting_point <- normal_dnds[3,2]
    usa_starting_point <- normal_dnds[4,2]

# make table to calculate difference in dnds
    dnds_LOH <- filter(dnds_summary_new, LOH == "LOH")
    dnds_nonLOH <- filter(dnds_summary_new, LOH == "nonLOH")
    dnds_BOTH <- full_join(dnds_LOH, dnds_nonLOH, by = c("SNVs", "stringeny","count")


pdf("dnds_by_LOH_thresholding.pdf", width=10, height=10)
# filter(dnds_summary_new, stringency == "stringent") %>%
#     ggplot(aes(x=as.integer(count),y=mle,ymin=cilow,ymax=cihigh,color=as.factor(SNVs),shape=as.factor(LOH)))+
#         geom_pointrange()+  
#         geom_line()+
#         scale_color_manual("SNVs from (excluding other LOH)",values=c("red","blue"))+
#         geom_abline(aes(intercept=1, slope=0), color="black", linetype="dashed")+
#         geom_vline(aes(xintercept = 10), color="black", linetype="dashed")+
#         geom_vline(aes(xintercept = 25), color="black", linetype="dashed")+
#         geom_abline(aes(intercept=pei_starting_point, slope=0), color="red", linetype="dashed")+
#         geom_abline(aes(intercept=usa_starting_point, slope=0), color="blue", linetype="dashed")+
#         xlab("Threshold for calling LOH out of 50")+
#         ylab("dNdS (corrected: 1 = neutral)")+
#         theme_classic() +
#             theme(axis.text=element_text(size=12,face="bold"),
#             axis.title=element_text(size=16,face="bold"),
#             text=element_text(size=18,face="bold")) +
#         ggtitle("dNdS - determining optimal LOH cutoff (stringent SNVs)")
filter(dnds_summary_new, stringency == "regular") %>%
    ggplot(aes(x=as.integer(count),y=mle,ymin=cilow,ymax=cihigh,color=as.factor(SNVs),shape=as.factor(LOH)))+
        geom_pointrange()+  
        geom_line()+
        scale_color_manual("SNVs from (excluding other LOH)",values=c("red","blue"))+
        geom_abline(aes(intercept=1, slope=0), color="black", linetype="dashed")+
        geom_vline(aes(xintercept = 10), color="black", linetype="dashed")+
        geom_vline(aes(xintercept = 25), color="black", linetype="dashed")+
        geom_abline(aes(intercept=pei_starting_point, slope=0), color="red", linetype="dashed")+
        geom_abline(aes(intercept=usa_starting_point, slope=0), color="blue", linetype="dashed")+
        xlab("Threshold for calling LOH out of 50")+
        ylab("dNdS (corrected: 1 = neutral)")+
        theme_classic() +
            theme(axis.text=element_text(size=12,face="bold"),
            axis.title=element_text(size=16,face="bold"),
            text=element_text(size=18,face="bold")) +
        ggtitle("dNdS - determining optimal LOH cutoff")
dev.off()

