# Original file here: C:\Users\shart\Metzger Lab Dropbox\Sam_data\Mya_genome\SNPs\somatypus_output_sigS_downstream.r

library(sigfit)
library(tidyverse)
library(lsa)
data(cosmic_signatures_v2)
data(cosmic_signatures_v3)

setwd("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/sigfit")

# Load helmsman output
    helmsman_output_file <- "/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/helmsman/subtype_count_matrix.txt"
    helmsman <- read.table(helmsman_output_file,
                                    sep="\t",
                    skip = 1)
    rownames(helmsman) <- helmsman$V1
    helmsman <- select(helmsman, -V1)
    colnames(helmsman) <- colnames(cosmic_signatures_v2)
    helmsman <- helmsman[order(row.names(helmsman)), ]
    head(helmsman)


#bins from helmsman run
    keeps <- c("allCanyH", "allCnoH", "allPEInoUSAnoH", "allUSAnoPEInoH",
    "H_P", "H_R", "H_RP", "H_RU", "H_RUP", "H_U", "H_UP", 
    "multiplePEInoUSAnoH", "multipleUSAnoPEInoH",
    "multipleP12noUSAnoH", "multipleP13noUSAnoH", "multipleP23noUSAnoH",
    "multipleU1234noPEInoH", "multipleU1235noPEInoH", "multipleU123noPEInoH", "multipleU1245noPEInoH", "multipleU124noPEInoH", "multipleU12noPEInoH", "multipleU1345noPEInoH", "multipleU134noPEInoH", "multipleU13noPEInoH", "multipleU14noPEInoH", "multipleU2345noPEInoH", "multipleU234noPEInoH", "multipleU23noPEInoH", "multipleU24noPEInoH", "multipleU34noPEInoH",
    "subsetCanyH", "subsetCPandUnoH",
    "uniqP1noUSAnoH", "uniqP2noUSAnoH", "uniqP3noUSAnoH", "uniqU1noPEInoH", "uniqU2noPEInoH", "uniqU3noPEInoH", "uniqU4noPEInoH", "uniqU5noPEInoH")


# correct for missing rows due to no mutations
    for(i in keeps){
        for(j in c("CDS","exon","five_prime_UTR","gene","genome","three_prime_UTR")){
            name<-paste(i,j,sep=".")
            if(i == "allCanyH" & j=="CDS"){
                helmsman2 <- subset(helmsman, rownames(helmsman) == name)
            }
            else{
                if(nrow(subset(helmsman, rownames(helmsman) == name))>0){
                    helmsman2 <- rbind(helmsman2,subset(helmsman, rownames(helmsman) == name))
                } 
                else{
                    print(paste(i,j,"didn't have any mutations"))
                    helmsman2[nrow(helmsman2) + 1,] <- rep(0,96)
                    }
            }
        }
    }

# Load mutational oppertunity data
    mutational_oppertunities_many <- readRDS("/ssd3/Mar_genome_analysis/mut_sig/trinuc_freq/mutational_oppertunities_all.rds")
    mutopps <- mutational_oppertunities_many[rep(1:6, nrow(helmsman2)/6),]

# # EXTRACT SIGNATURES - TEST BEST ONE
#     helmsman_extr <- extract_signatures(counts = helmsman2,
#                                         nsignatures = 2:7,
#                                         opportunities =  mutopps,
#                                         iter = 1000, 
#                                         seed = 1756)
# # Estimated best number of signatures: 3
# # used 4 signatures since this allows us to discern a signature1-like signature with distinct NCpG>NTpG

# # Plot outputs
#     setwd("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/sigfit")
#     pdf("goodness_of_fit_generegions.pdf")
#     plot_gof(helmsman_extr)
#     dev.off()
#     signatures_3_prelim <- retrieve_pars(helmsman_extr[[3]], par = "signatures")
#     signatures_4_prelim <- retrieve_pars(helmsman_extr[[4]], par = "signatures")
#     plot_spectrum(signatures_3_prelim, pdf_path = "extract_3sigs_generegions_prelim.pdf")
#     plot_spectrum(signatures_4_prelim, pdf_path = "extract_4sigs_generegions_prelim.pdf")

# # Compare to known signatures, rename and save RDS file
#     match_signatures(signatures_4_prelim, cosmic_signatures_v2)
#     # 1 => 5, 2 => 1, 3 => 30, 4 => 9
#     match_signatures(signatures_4_prelim, cosmic_signatures_v3)
#     # 1 => 5, 2 => 1, 3 => 45, 4 => 12 
#     # note these are offset due to multiple sub-signatures (eg. 7a,7b,7c), really coorespond to 5,1,40,9
#     signames <- c("sig5", "sig1","sig40", "sigS")
#     rownames(signatures_4_prelim$mean) <- signames
#     rownames(signatures_4_prelim$upper_95) <- signames
#     rownames(signatures_4_prelim$lower_95) <- signames 
#     saveRDS(signatures_4_prelim, file = "signatures_4_prelim_extracted.rds")


# Run full extraction
    helmsman_extr_4full <- extract_signatures(counts = helmsman2,
                                            nsignatures = 4,
                                            opportunities =  mutopps,
                                            iter = 10000, 
                                            seed = 1756)
    signatures_4 <- retrieve_pars(helmsman_extr_4full, par = "signatures")
    # saveRDS(helmsman_extr_4full, file = "signatures_4_extracted_full.rds")
    # signatures_4 <- readRDS(file = "/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/sigfit/signatures_4_extracted_full.rds")
    # plot_spectrum(signatures_4, pdf_path = "extract_4sigs_generegions.pdf")

############################################################################################################################################## 8/10/21 left off
# Match known signatures
    match_signatures(signatures_4, cosmic_signatures_v2)
    # 1 => 5, 2 => 1, 3 => 30, 4 => 9
    match_signatures(signatures_4, cosmic_signatures_v3)
    # 1 => 5, 2 => 1, 3 => 45, 4 => 12
    # note these are offset due to multiple sub-signatures (eg. 7a,7b,7c), really coorespond to 5,1,40,9
    cosine(as.numeric(signatures_4$mean[1,]), as.numeric(cosmic_signatures_v2[5,])) 
    # 0.8746322
    cosine(as.numeric(signatures_4$mean[2,]), as.numeric(cosmic_signatures_v2[1,])) 
    # 0.9119926
    cosine(as.numeric(signatures_4$mean[3,]), as.numeric(cosmic_signatures_v2[8,])) 
    # 0.7234261
    cosine(as.numeric(signatures_4$mean[4,]), as.numeric(cosmic_signatures_v2[9,])) 
    # 0.7108884
    cosine(as.numeric(signatures_4$mean[1,]), as.numeric(cosmic_signatures_v3[5,])) 
    # 0.8378334
    cosine(as.numeric(signatures_4$mean[2,]), as.numeric(cosmic_signatures_v3[1,])) 
    # 0.8115095
    cosine(as.numeric(signatures_4$mean[3,]), as.numeric(cosmic_signatures_v3[45,])) 
    # 0.7519232
    cosine(as.numeric(signatures_4$mean[4,]), as.numeric(cosmic_signatures_v3[12,])) 
    # 0.6984928

# Rename and save RDS file
    signames <- c("sig5", "sig1","sig40", "sigS")
    rownames(signatures_4$mean) <- signames
    rownames(signatures_4$upper_95) <- signames
    rownames(signatures_4$lower_95) <- signames 
    plot_spectrum(signatures_4, pdf_path = "extract_4sigs_generegions.pdf")
    saveRDS(signatures_4, file = "signatures_4_extracted_full.rds")
    #signatures_4 <- readRDS(file = "/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/sigfit/signatures_4_extracted_full.rds")

# fit data for the generegions
    helmsman_4sigs_fit <- fit_signatures(counts = helmsman2,  
                                    signatures = signatures_4,
                                    opportunities =  mutopps,
                                    iter = 2000, 
                                    warmup = 1000, 
                                    chains = 1, 
                                    seed = 1756)
    plot_exposures(mcmc_samples = helmsman_4sigs_fit , pdf_path = "exposures_4sigs_generegions.pdf")
    plot_reconstruction(mcmc_samples = helmsman_4sigs_fit, pdf_path = "reconstruction_4sigs_generegions.pdf")
    saveRDS(helmsman_4sigs_fit, file = "fit_4sigs_generegions.rds")
    # helmsman_4sigs_fit <- readRDS(file = "/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/sigfit/fit_4sigs_generegions.rds")

# export fitting outputs to play with in excel
    helmsman_exposures <- retrieve_pars(helmsman_4sigs_fit, par = "exposures")
    sig_exp <- helmsman_exposures$mean
    sig_expU <- helmsman_exposures$upper_95
    sig_expL <- helmsman_exposures$lower_95
    signames <- c("sig5", "sig1","sig40", "sigS")
    colnames(sig_exp) <- signames
    colnames(sig_expU) <- c("sig5u", "sig1u","sig40u", "sigSu")
    colnames(sig_expL) <- c("sig5l", "sig1l","sig40l", "sigSl")
    counts <- rowSums(helmsman2)
    output_4sigs <- cbind(sig_exp, sig_expU, sig_expL, counts)
    write.table(output_4sigs, "fit_4sigs_generegions.tsv", sep="\t")

# Plot raw spectra

    genome_only <- helmsman2[(1:(nrow(helmsman2)/6)*6-1),]
    genome_only[nrow(genome_only) + 1,] <- helmsman2[c("H_P.genome","H_U.genome","H_R.genome","H_RP.genome","H_RU.genome","H_RUP.genome","H_UP.genome"),] %>% colSums()
    rownames(genome_only)[rownames(genome_only) == "42"] <- "anyH.genome"

    plot_spectrum(genome_only, pdf_path = "genome_raw_spectra.pdf")
    helmsman_freq <- genome_only/mutational_oppertunities_many[rep(5, nrow(genome_only)),]
    helmsman_freq <- helmsman_freq/rowSums(helmsman_freq)
    helmsman_freq <- filter(helmsman_freq, !is.na(get("TTC>TGC")))
    plot_spectrum(helmsman_freq, pdf_path = "genome_raw_spectra_mutoppscorrected.pdf")

    anyH_raw <- helmsman2[c("H_P.genome","H_U.genome","H_R.genome","H_RP.genome","H_RU.genome","H_RUP.genome","H_UP.genome"),] %>% colSums()
    H_U_raw <- helmsman2[c("H_U.genome","H_RU.genome","H_RUP.genome","H_UP.genome"),] %>% colSums()
    H_P_raw <- helmsman2[c("H_P.genome","H_RP.genome","H_RUP.genome","H_UP.genome"),] %>% colSums()
    rownames(anyH_raw) <- "anyH"
    rbind(anyH_raw,H_U_raw,H_P_raw)


    exp_names <- tibble::rownames_to_column(helmsman_exposures$mean, var = "rowname") %>%
        separate(rowname, c("sample","region"), sep = "\\.") %>%
        #filter(B == "genome") %>%
        select(sample,region)
    exp_mean <- helmsman_exposures$mean
    exp_upper <- helmsman_exposures$upper_95%>%
        rename(sigS_U=sigS,sig1_U=sig1,sig5_U=sig5,sig40_U=sig40)
    exp_lower <- helmsman_exposures$lower_95 %>%
        rename(sigS_L=sigS,sig1_L=sig1,sig5_L=sig5,sig40_L=sig40)
    mutation_counts <- rowSums(helmsman2) %>% as.data.frame() 
    colnames(mutation_counts) <- "count"
    exposures_total_B <- cbind(exp_names,mutation_counts,exp_mean,exp_upper,exp_lower) 
    exposures_genome <- filter(exposures_total_B,region=="genome")


    # pdf("sig1_test.pdf",width=30, height=10)
    # for(i in c("sigS","sig1")){
    # plot1 <- ggplot(exposures_genome) + 
    #     geom_col(aes(x=sample,y=get(i)))+
    #     geom_linerange(aes(x=sample,y=get(i),ymin=get(paste0(i,"_L")),ymax=get(paste0(i,"_U"))))+
    #     ylab(paste(i, "exposure"))+
    #     theme_classic() +
    #         theme(axis.text=element_text(size=12,face="bold"),
    #         axis.title=element_text(size=16,face="bold"),
    #         text=element_text(size=18,face="bold")) +
    #         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    #     ggtitle(paste(i, "exposure"))
    #     print(plot1)
    # }
    # dev.off()


# Load LOH data, select only desired rows and fit to signatures
    helmsman_output_file_LOH <- "/ssd3/Mar_genome_analysis/LOH/july_2021/output/new/helmsman/subtype_count_matrix.txt"
    helmsman_LOH <- read.table(helmsman_output_file_LOH,
                                    sep="\t",
                    skip = 1)
    rownames(helmsman_LOH) <- helmsman_LOH$V1
    helmsman_LOH <- select(helmsman_LOH, -V1)
    colnames(helmsman_LOH) <- colnames(cosmic_signatures_v2)
    helmsman_LOH <- helmsman_LOH[order(row.names(helmsman_LOH)), ]
    helmsman_LOH <- helmsman_LOH[c("PEI_LOH_10_hetero_counts.bed.merge.USAloh", "PEI_LOH_10_hetero_counts.bed.merge.USAsom","USA_LOH_10_hetero_counts.bed.merge.PEIloh", "USA_LOH_10_hetero_counts.bed.merge.PEIsom"),]
    head(helmsman_LOH)
    mutopps_LOH <- mutational_oppertunities_many[rep(5, nrow(helmsman_LOH)),]
    helmsman_fit_LOH <- fit_signatures(counts = helmsman_LOH,  
                                    signatures = signatures_4,
                                    opportunities =  mutopps_LOH,
                                    iter = 2000, 
                                    warmup = 1000, 
                                    chains = 1, 
                                    seed = 1756)
    helmsman_exposures_LOH <- retrieve_pars(helmsman_fit_LOH, par = "exposures")
    exp_names <- tibble::rownames_to_column(helmsman_exposures_LOH$mean, var = "rowname") %>%
        separate(rowname, c("region",NA,NA,"sample"), sep = "\\.") %>%
        select(sample,region)
    exp_mean <- helmsman_exposures_LOH$mean
    exp_upper <- helmsman_exposures_LOH$upper_95%>%
        rename(sigS_U=sigS,sig1_U=sig1,sig5_U=sig5,sig40_U=sig40)
    exp_lower <- helmsman_exposures_LOH$lower_95 %>%
        rename(sigS_L=sigS,sig1_L=sig1,sig5_L=sig5,sig40_L=sig40)
    mutation_counts <- rowSums(helmsman_LOH) %>% as.data.frame() 
    colnames(mutation_counts) <- "count"
    exposures_total_LOH <- cbind(exp_names,mutation_counts,exp_mean,exp_upper,exp_lower) 

    exposures_genome_and_LOH <- rbind(exposures_genome,exposures_total_LOH)
    saveRDS(exposures_genome_and_LOH, "/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/sigfit/exposures_genome_and_LOH.rds")
    #exposures_genome_and_LOH <- readRDS("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/sigfit/exposures_genome_and_LOH.rds")

# Plot noLOH spectra
    plot_spectrum(helmsman_LOH, pdf_path = "noLOH_raw_spectra.pdf")
    LOH_freq <- helmsman_LOH/mutational_oppertunities_many[rep(5, nrow(helmsman_LOH)),]
    LOH_freq <- LOH_freq/rowSums(LOH_freq)
    plot_spectrum(LOH_freq, pdf_path = "noLOH_raw_spectra_mutoppscorrected.pdf")

exposures_genome_plot <- exposures_genome_and_LOH[c("H_P.genome","H_U.genome","H_RUP.genome",
                                            "allCanyH.genome","allCnoH.genome","allPEInoUSAnoH.genome","allUSAnoPEInoH.genome","USA_LOH_10_hetero_counts.bed.merge.PEIsom","PEI_LOH_10_hetero_counts.bed.merge.USAsom",
                                            "multipleU1234noPEInoH.genome","multipleU123noPEInoH.genome","multipleU12noPEInoH.genome","multipleP12noUSAnoH",
                                            "uniqP1noUSAnoH.genome","uniqP2noUSAnoH.genome","uniqP3noUSAnoH.genome",
                                            "uniqU1noPEInoH.genome","uniqU2noPEInoH.genome","uniqU3noPEInoH.genome","uniqU4noPEInoH.genome","uniqU5noPEInoH.genome"),]
exposures_genome_plot$sample <- factor(c("PEI healthy","USA healthy","All healthy clams",
                                            "All cancer, any healthy","All cancer, no healthy","PEI sublineage","USA sublineage","PEI sublineage, nonLOH","USA sublineage, nonLOH",
                                            #"Maine sublineage","Maine sub-sublineage","Maine sub-sub-sublineage","PEI sub-sublineage",
                                            "USA 1&2&3&4 only","USA 1&2&3 only","USA 1&2 only","PEI 1&2 only",
                                            "Only PEI1","Only PEI2","Only PEI3",
                                            "Only USA1","Only USA2","Only USA3","Only USA4","Only USA5"))
exposures_genome_plot$sample <- factor(exposures_genome_plot$sample, levels = c("PEI healthy","USA healthy","All healthy clams",
                                            "All cancer, any healthy","All cancer, no healthy","PEI sublineage","USA sublineage","PEI sublineage, nonLOH","USA sublineage, nonLOH",
                                            #"Maine sublineage","Maine sub-sublineage","Maine sub-sub-sublineage","PEI sub-sublineage",
                                            "USA 1&2&3&4 only","USA 1&2&3 only","USA 1&2 only","PEI 1&2 only",
                                            "Only PEI1","Only PEI2","Only PEI3",
                                            "Only USA1","Only USA2","Only USA3","Only USA4","Only USA5"))
exposures_genome_plot$color <- c("Healthy","Healthy","Healthy",
                                            "Cancer","Cancer","PEI","USA","PEI","USA",
                                            "USA","USA","USA","PEI",
                                            "PEI","PEI","PEI",
                                            "USA","USA","USA","USA","USA")    
exposures_genome_plot2 <- exposures_genome_and_LOH[c("H_P.genome","H_U.genome","H_RUP.genome",
                                            "allCanyH.genome","allCnoH.genome","allPEInoUSAnoH.genome","allUSAnoPEInoH.genome","USA_LOH_10_hetero_counts.bed.merge.PEIsom","PEI_LOH_10_hetero_counts.bed.merge.USAsom"),]
exposures_genome_plot2$sample <- factor(c("PEI_H","USA_H","All_H",
                                                "All_C_any_H","All_C_no_H","PEI_only","USA_only","PEI_only_nonLOH","USA_only_nonLOH"))
exposures_genome_plot2$sample <- factor(exposures_genome_plot2$sample, levels = c("PEI_H","USA_H","All_H",
                                                "All_C_any_H","All_C_no_H","PEI_only","USA_only","PEI_only_nonLOH","USA_only_nonLOH"))
exposures_genome_plot2$color <- c("Healthy","Healthy","Healthy",
                                            "Cancer","Cancer","PEI","USA","PEI","USA")    
setwd("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/sigfit")
pdf("Signature_exposures.pdf",width=30, height=10)
for(i in c("sigS","sig1")){ #,"sig5","sig40"
plot1 <- ggplot(exposures_genome_plot) + 
    geom_col(aes(x=sample,y=get(i),fill=color))+
    geom_linerange(aes(x=sample,y=get(i),ymin=get(paste0(i,"_L")),ymax=get(paste0(i,"_U"))))+
    geom_text(aes(x=sample,y=-0.02,label = count))+
    ylab(paste(i, "exposure"))+
    scale_fill_manual(values=c("grey", "black", "red", "blue"))+
    theme_classic() +
        theme(axis.text=element_text(size=16,face="bold"),
        axis.title=element_text(size=16,face="bold"),
        text=element_text(size=18,face="bold")) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ggtitle(paste(i, "exposure"))
    print(plot1)
}
dev.off()
# NEW PLOTS
    setwd("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/sigfit")
    pdf("Signature_exposures_Sand1.pdf",width=8, height=4)
    for(i in c("sigS","sig1")){ #,"sig5","sig40"
    plot1 <- ggplot(exposures_genome_plot) + 
        #geom_col(aes(x=sample,y=get(i),fill=color))+
        #geom_linerange(aes(x=sample,y=get(i),ymin=get(paste0(i,"_L")),ymax=get(paste0(i,"_U"))))+
        geom_pointrange(aes(x=sample,y=get(i),ymin=get(paste0(i,"_L")),ymax=get(paste0(i,"_U")),color=color))+
        #geom_text(aes(x=sample,y=-0.02,label = count))+
        ylab(paste(i, "SNV Fraction"))+
        xlab(NULL)+
        #scale_fill_manual(values=c("grey", "black", "red", "blue"))+
        scale_color_manual(values=c("grey", "black", "red", "blue"))+
        theme_classic() +
            theme(axis.text=element_text(size=8,face="bold"),
            axis.title=element_text(size=8,face="bold"),
            text=element_text(size=8,face="bold")) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))#+
        #ggtitle(paste(i, "SNV Fraction"))
        print(plot1)
    }
    dev.off()
###
pdf("Signature_exposures_condensed.pdf",width=7.5, height=10)
for(i in c("sigS","sig1")){ #,"sig5","sig40"
plot1 <- ggplot(exposures_genome_plot2) + 
    #geom_col(aes(x=sample,y=get(i),fill=color))+
    geom_pointrange(aes(x=sample, y=get(i), ymin=get(paste0(i,"_L")), ymax=get(paste0(i,"_U")), color=color), fatten = 3, size = 3)+ #, size = 1))+
    geom_errorbar(aes(x=sample, y=get(i), ymin=get(paste0(i,"_L")), ymax=get(paste0(i,"_U"))),color="grey20",width=0.5)+
    geom_text(aes(x=sample,y=-0.02,label = count))+
    geom_abline(aes(intercept=mean(c(exposures_genome_plot2["USA_LOH_10_hetero_counts.bed.merge.PEIsom",i],exposures_genome_plot2["PEI_LOH_10_hetero_counts.bed.merge.USAsom",i])), slope=0), color="black", linetype="dashed")+
    geom_abline(aes(intercept=mean(c(exposures_genome_plot2["H_P.genome",i],exposures_genome_plot2["H_U.genome",i])), slope=0), color="black", linetype="dashed")+
    geom_abline(aes(intercept=exposures_genome_plot2["allCnoH.genome",i], slope=0), color="black", linetype="dashed")+
    ylab(paste(i, "exposure"))+
    scale_color_manual(values=c("grey", "black", "red", "blue"))+
    theme_classic() +
        theme(axis.text=element_text(size=24,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        text=element_text(size=24,face="bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none") +
    ggtitle(paste(i, "exposure"))
    print(plot1)
}
dev.off()

# OLD NOTES AND CALCULATIONS

# # estimate % AllCnoH that are somatic
#     avg_1H_sigS <- (exposures_genome_plot2["H_P.genome","sigS"] + exposures_genome_plot2["H_U.genome","sigS"])/2
#     avg_somatic_sigS <- (exposures_genome_plot2["USA_LOH_10_hetero_counts.bed.merge.PEIsom","sigS"] + exposures_genome_plot2["PEI_LOH_10_hetero_counts.bed.merge.USAsom","sigS"])/2
#     allCnoH_sigS <- exposures_genome_plot2["allCnoH.genome","sigS"]
#     allCnoH_somatic_fraction_estimate = (allCnoH_sigS - avg_1H_sigS)/(avg_somatic_sigS - avg_1H_sigS) # 0.06674076
#     (0.05095969 - 0.02068163)/(0.4743483 -  0.02068163) # 0.06674076 same as above
#     (0.05095969 - 0)/(0.4743483 -  0) # 0.107431
#     (0.05095969 - 0)/(0.38093483 -  0) # 0.1337753
#     (0.05159662 - 0)/(0.38093483 -  0) # 0.1354474 # MAX
#     (0.05095969 - 0.02996840)/(0.4743483 -  0.02996840) # 0.04723726
#     (0.05095969 - 0.02996840)/(0.48404515 -  0.02996840) # 0.04622851
#     (0.05023588 - 0.02996840)/(0.48404515 -  0.02996840) # 0.04463448 # MIN
#     allCnoH_sigS_somatic_fraction_estimate <- allCnoH_sigS-avg_1H_sigS # 0.03085008
#     0.05095969-0.02068163 # 0.03085008
#     0.05159662-0 # 0.05159662
#     0.05023588-0.02996840 # 0.02026748



#     allCnoH_sigS_somatic_fraction_estimate* 2596657 # 80107.09 sigS mutations in trunk
#     allCnoH_sigS_somatic_fraction_estimate* 2596657 * 1000000000/full_genome_size # 65873.84 sigS mutations in trunk per Gb
#     allCnoH_sigS_somatic_fraction_estimate* 2596657 * 1000000000/full_genome_size # 65873.84 sigS mutations in trunk per Gb
#     # from other calculation in somatypus_output-pairwise.r: sigS_rate = 430
#     allCnoH_sigS_somatic_fraction_estimate* 2596657 * 1000000000/full_genome_size / 430 #  153.195 years
#     allCnoH_sigS_somatic_fraction_estimate* 2596657 * 1000000000/full_genome_size / 240 # 274.4743 years
#     allCnoH_sigS_somatic_fraction_estimate* 2596657 * 1000000000/full_genome_size / 620 # 106.2481 years
#     (0.05159662-0)* 2596657 * 1000000000/full_genome_size / 240 # 459.057 years
#     (0.05023588-0.02996840)* 2596657 * 1000000000/full_genome_size / 620 # 69.80149 years

# how much of genome are we excluding with LOH calls
    LOH_genome_size <- read.table("/ssd3/Mar_genome_analysis/LOH/july_2021/output/new/BOTH_SUBLINEAGES_LOH_10_counts.bed") %>% mutate(V4 = V3-V2) %>% .$V4 %>%as.numeric() %>% sum()
    full_genome_size <- read.table("/ssd3/Mar_genome_analysis/genomes/Mar.3.4.6.p1_Q30Q30A.fasta.fai") %>% .$V2 %>%as.numeric() %>% sum()
    nonLOH_genome_size = full_genome_size - LOH_genome_size 
    perGbcorrection = nonLOH_genome_size/1000000000
    full_genome_correction <- full_genome_size/nonLOH_genome_size

    full_genome_correction*exposures_genome_plot["USA_LOH_10_hetero_counts.bed.merge.PEIsom","count"]
    full_genome_correction*exposures_genome_plot["PEI_LOH_10_hetero_counts.bed.merge.USAsom","count"]


# Extract sigs for just somatic samples
    # helmsman_genome <- helmsman2[(1:(nrow(helmsman2)/6))*6-1,]
    # keeps_somatic <- c("allPEInoUSAnoH", "allUSAnoPEInoH",
    # "multipleP12noUSAnoH",
    # "multipleU1234noPEInoH","multipleU123noPEInoH", "multipleU12noPEInoH",
    # "uniqP1noUSAnoH", "uniqP2noUSAnoH", "uniqP3noUSAnoH", "uniqU1noPEInoH", "uniqU2noPEInoH", "uniqU3noPEInoH", "uniqU4noPEInoH", "uniqU5noPEInoH")
    # helmsman_somatic <- helmsman_genome[keeps_somatic,]
    # helmsman_genome[keeps_somatic,]
    # mutopps_somatic <- mutational_oppertunities_many[rep(5, nrow(helmsman_somatic)),]
    # helmsman_somatic_extr <- extract_signatures(counts = helmsman_somatic,
    #                             nsignatures = 2:4,
    #                             #opportunities =  mutopps_somatic,
    #                             iter = 1000, 
    #                             seed = 1756)
    # helmsman_genome_extr <- extract_signatures(counts = helmsman_genome,
    #                             nsignatures = 2:4,
    #                             #opportunities =  mutopps_somatic,
    #                             iter = 1000, 
    #                             seed = 1756)
    # helmsman_somatic_extr2 <- fit_extract_signatures(counts = helmsman_somatic,
    #                             signatures = helmsman_genome["allCanyH.genome",]/sum(helmsman_genome["allCanyH.genome",]),
    #                             num_extra_sigs=3,
    #                             #opportunities =  mutopps_somatic,
    #                             iter = 1000, 
    #                             seed = 1756)    

    # helmsman_somatic_extr3 <- extract_signatures(counts = helmsman_somatic,
    #                             nsignatures = 2:4,
    #                             #opportunities =  mutopps_somatic,
    #                             iter = 1000, 
    #                             seed = 1756)
    # # Estimated best number of signatures: 3
    # setwd("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/sigfit/somatic")
    # pdf("goodness_of_fit_somatic_bins_nomutopps.pdf")
    # plot_gof(helmsman_somatic_extr)
    # dev.off()
    # signatures_3_somatic <- retrieve_pars(helmsman_somatic_extr[[3]], par = "signatures")
    # signatures_4_somatic <- retrieve_pars(helmsman_somatic_extr[[4]], par = "signatures")
    # plot_spectrum(signatures_3_somatic, pdf_path = "extract_3sigs_somatic_nomutopps.pdf")
    # plot_spectrum(signatures_4_somatic, pdf_path = "extract_4sigs_somatic_nomutopps.pdf")
    # match_signatures(signatures_3_somatic, cosmic_signatures_v3)
    # # 1 => 4, 2 => 9, 3 => 5, 4 => 6
    # match_signatures(signatures_4_somatic, cosmic_signatures_v3)
    # # 1 => 4, 2 => 12, 3 => 5, 4 => 6
    # signatures_3_somatic <- retrieve_pars(helmsman_somatic_extr2, par = "signatures")
    # plot_spectrum(signatures_3_somatic, pdf_path = "extract_3sigs_somatic_nomutopps2.pdf")
    # match_signatures(signatures_3_somatic, cosmic_signatures_v2)

    # signatures_3_genome <- retrieve_pars(helmsman_genome_extr[[3]], par = "signatures")
    # signatures_4_genome <- retrieve_pars(helmsman_genome_extr[[4]], par = "signatures")
    # plot_spectrum(signatures_3_genome, pdf_path = "extract_3sigs_genome_nomutopps.pdf")
    # plot_spectrum(signatures_4_genome, pdf_path = "extract_4sigs_genome_nomutopps.pdf")    
    # match_signatures(signatures_3_genome, cosmic_signatures_v3)
    # match_signatures(signatures_4_genome, cosmic_signatures_v3)

# CDS vs non-CDS
    genome_all <-output_4sigs[(1:(nrow(output_4sigs)/6)*6-1),c(2,6,10,13)]
    cds_all <- output_4sigs[(1:(nrow(output_4sigs)/6)*6-5),c(2,6,10,13)]
    # colnames(cds_all) <- paste0(colnames(cds_all),"_CDS")
    # colnames(genome_all) <- paste0(colnames(genome_all),"_genome")
    bins_all <- rbind(genome_all,cds_all)
    bins_all$rowname <- c(rownames(genome_all),rownames(cds_all))
    bins_all <- separate(bins_all, rowname, c("name","region"), sep = "\\.")
    bins_all <- bins_all[c("H_P.genome","H_U.genome","H_RUP.genome",
                                                "allCanyH.genome","allCnoH.genome","allPEInoUSAnoH.genome","allUSAnoPEInoH.genome",
                                                "multipleU1234noPEInoH.genome","multipleU123noPEInoH.genome","multipleU12noPEInoH.genome","multipleP12noUSAnoH.genome",
                                                "uniqP1noUSAnoH.genome","uniqP2noUSAnoH.genome","uniqP3noUSAnoH.genome",
                                                "uniqU1noPEInoH.genome","uniqU2noPEInoH.genome","uniqU3noPEInoH.genome","uniqU4noPEInoH.genome","uniqU5noPEInoH.genome",
                            "H_P.CDS","H_U.CDS","H_RUP.CDS",
                                                "allCanyH.CDS","allCnoH.CDS","allPEInoUSAnoH.CDS","allUSAnoPEInoH.CDS",
                                                "multipleU1234noPEInoH.CDS","multipleU123noPEInoH.CDS","multipleU12noPEInoH.CDS","multipleP12noUSAnoH.CDS",
                                                "uniqP1noUSAnoH.CDS","uniqP2noUSAnoH.CDS","uniqP3noUSAnoH.CDS",
                                                "uniqU1noPEInoH.CDS","uniqU2noPEInoH.CDS","uniqU3noPEInoH.CDS","uniqU4noPEInoH.CDS","uniqU5noPEInoH.CDS"),]
    bins_all$names <- factor(c("PEI healthy","USA healthy","All healthy clams",
                                                "All cancer, any healthy","All cancer, no healthy","PEI sublineage","USA sublineage",
                                                #"Maine sublineage","Maine sub-sublineage","Maine sub-sub-sublineage","PEI sub-sublineage",
                                                "USA 1&2&3&4 only","USA 1&2&3 only","USA 1&2 only","PEI 1&2 only",
                                                "Only PEI1","Only PEI2","Only PEI3",
                                                "Only USA1","Only USA2","Only USA3","Only USA4","Only USA5",
                                "PEI healthy","USA healthy","All healthy clams",
                                                "All cancer, any healthy","All cancer, no healthy","PEI sublineage","USA sublineage",
                                                #"Maine sublineage","Maine sub-sublineage","Maine sub-sub-sublineage","PEI sub-sublineage",
                                                "USA 1&2&3&4 only","USA 1&2&3 only","USA 1&2 only","PEI 1&2 only",
                                                "Only PEI1","Only PEI2","Only PEI3",
                                                "Only USA1","Only USA2","Only USA3","Only USA4","Only USA5"))
    bins_all$names <- factor(bins_all$names, levels = c("PEI healthy","USA healthy","All healthy clams",
                                                "All cancer, any healthy","All cancer, no healthy","PEI sublineage","USA sublineage",
                                                #"Maine sublineage","Maine sub-sublineage","Maine sub-sub-sublineage","PEI sub-sublineage",
                                                "USA 1&2&3&4 only","USA 1&2&3 only","USA 1&2 only","PEI 1&2 only",
                                                "Only PEI1","Only PEI2","Only PEI3",
                                                "Only USA1","Only USA2","Only USA3","Only USA4","Only USA5"))
    #bins_all$region <- factor(bins_all$region, levels = c(CDS, genome))
    bins_all$color <- c("Healthy","Healthy","Healthy",
                                                "Cancer","Cancer","PEI","USA",
                                                "USA","USA","USA","PEI",
                                                "PEI","PEI","PEI",
                                                "USA","USA","USA","USA","USA",
                        "Healthy","Healthy","Healthy",
                                                "Cancer","Cancer","PEI","USA",
                                                "USA","USA","USA","PEI",
                                                "PEI","PEI","PEI",
                                                "USA","USA","USA","USA","USA")    
    bins_all

    setwd("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/sigfit")
    pdf("Signature1_CDSvsGENOME.pdf",width=8, height=4)
    ggplot(bins_all) + 
        #geom_col(aes(x=sample,y=get(i),fill=color))+
        geom_pointrange(aes(x=names,y=sig1,ymin=sig1l,ymax=sig1u,color=color,shape= region), position=position_dodge(width=0.25))+
        #geom_text(aes(x=sample,y=-0.02,label = count))+
        ylab("Fraction of mutations attributed to sig1")+
        scale_color_manual(values=c("grey", "black", "red", "blue"))+
        theme_classic() +
            theme(axis.text=element_text(size=8,face="bold"),
            axis.title=element_text(size=8,face="bold"),
            text=element_text(size=8,face="bold")) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        ggtitle("Sig1: CDS vs full genome")
    dev.off()

# Count bias for single nucleotide mutations
    # sigS
    signatures_4$mean[4,1:16] %>% sum() # 0.2566099
    signatures_4$mean[4,17:32] %>% sum() # 0.04221173
    signatures_4$mean[4,33:48] %>% sum() # 0.2282434
    signatures_4$mean[4,49:64] %>% sum() # 0.09237359
    signatures_4$mean[4,65:80] %>% sum() # 0.1325403
    signatures_4$mean[4,81:96] %>% sum() # 0.2480211

    #Sig5
    signatures_4$mean[1,1:16] %>% sum()
    signatures_4$mean[1,17:32] %>% sum()
    signatures_4$mean[1,33:48] %>% sum()
    signatures_4$mean[1,49:64] %>% sum()
    signatures_4$mean[1,65:80] %>% sum()
    signatures_4$mean[1,81:96] %>% sum()
    #Sig1
    signatures_4$mean[2,1:16] %>% sum()
    signatures_4$mean[2,17:32] %>% sum()
    signatures_4$mean[2,33:48] %>% sum()
    signatures_4$mean[2,49:64] %>% sum()
    signatures_4$mean[2,65:80] %>% sum()
    signatures_4$mean[2,81:96] %>% sum()
    #Sig40
    signatures_4$mean[3,1:16] %>% sum()
    signatures_4$mean[3,17:32] %>% sum()
    signatures_4$mean[3,33:48] %>% sum()
    signatures_4$mean[3,49:64] %>% sum()
    signatures_4$mean[3,65:80] %>% sum()
    signatures_4$mean[3,81:96] %>% sum()

    healthy_freq <- helmsman_freq["allCanyH",]
    healthy_freq[1,1:16] %>% sum() # 0.162631
    healthy_freq[1,17:32] %>% sum() # 0.07369336
    healthy_freq[1,33:48] %>% sum() #  0.3410991
    healthy_freq[1,49:64] %>% sum() # 0.1523371
    healthy_freq[1,65:80] %>% sum() # 0.182695
    healthy_freq[1,81:96] %>% sum() # 0.08754439


# "humanize" Sig9 to compare to COSMIC - didn't look any better
    # mutational_oppertunities_many <- readRDS("/ssd3/Mar_genome_analysis/mut_sig/trinuc_freq/mutational_oppertunities_all.rds")
    # mutopps_genome <- mutational_oppertunities_many[5,]
    # mutopps_sigs <- mutational_oppertunities_many[rep(5, nrow(signatures_4$mean)),]
    # mutopps_sigs


    # signatures_4 <- readRDS(file = "/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/sigfit/signatures_4_extracted_full.rds")

    # # humanized_signatures1 <- convert_signatures(signatures_4, 
    # #                                        opportunities_from = mutopps_sigs,
    # #                                        opportunities_to = "human-genome")

    # #correct:
    # humanized_signatures <- convert_signatures(signatures_4, 
    #                                     opportunities_to = "human-genome")
    # unnormalized_signatures <- convert_signatures(signatures_4, 
    #                                     opportunities_to = mutopps_sigs)

    # plot_spectrum(signatures_4, pdf_path = "norm_to_clam_signatures.pdf")
    # plot_spectrum(unnormalized_signatures, pdf_path = "unnormalized_signatures.pdf")
    # # plot_spectrum(humanized_signatures1, pdf_path = "humanized_signatures_wrong.pdf")
    # plot_spectrum(humanized_signatures, pdf_path = "humanized_signatures.pdf")
    # humanized_signatures2

    # match_signatures(signatures_4, cosmic_signatures_v3)  
    # # 1 => 5, 2 => 1, 3 => 45, 4 => 12  
    #     cosine(as.numeric(signatures_4$mean[1,]), as.numeric(cosmic_signatures_v3[5,])) # 0.8378334
    #     cosine(as.numeric(signatures_4$mean[2,]), as.numeric(cosmic_signatures_v3[1,])) # 0.8115095
    #     cosine(as.numeric(signatures_4$mean[3,]), as.numeric(cosmic_signatures_v3[45,])) # 0.7519232
    #     cosine(as.numeric(signatures_4$mean[4,]), as.numeric(cosmic_signatures_v3[12,])) # 0.6984928
    # match_signatures(humanized_signatures, cosmic_signatures_v3)
    # # 1 => 3, 2 => 5, 3 => 45, 4 => 12
    #     cosine(as.numeric(humanized_signatures[1,]), as.numeric(cosmic_signatures_v3[3,])) # 0.7576426
    #     cosine(as.numeric(humanized_signatures[2,]), as.numeric(cosmic_signatures_v3[5,])) # 0.853395
    #     cosine(as.numeric(humanized_signatures[3,]), as.numeric(cosmic_signatures_v3[45,])) # 0.8746741
    #     cosine(as.numeric(humanized_signatures[4,]), as.numeric(cosmic_signatures_v3[12,]))  # 0.7695406
    # match_signatures(unnormalized_signatures, cosmic_signatures_v3)    
    # # 1 => 46, 2 => 5, 3 => 45, 4 => 12
    #     cosine(as.numeric(unnormalized_signatures[1,]), as.numeric(cosmic_signatures_v3[46,])) # 0.7562813
    #     cosine(as.numeric(unnormalized_signatures[2,]), as.numeric(cosmic_signatures_v3[5,])) # 0.7693429
    #     cosine(as.numeric(unnormalized_signatures[3,]), as.numeric(cosmic_signatures_v3[45,])) # 0.8102292
    #     cosine(as.numeric(unnormalized_signatures[4,]), as.numeric(cosmic_signatures_v3[12,])) # 0.748946

            