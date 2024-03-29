# original file here: C:\Users\shart\Metzger Lab Dropbox\Sam_data\Mya_genome\SNPs\trunk_estimates3.R

library(tidyverse)
library(sigfit)
library(dndscv)
setwd("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/pairwise/noLOH")

############################################################### SigS ############################################################
# Load mutation counts: rare clam alleles and allCnoH
    helmsman_output_file <- "/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/pairwise/noLOH/helmsman_all_output/subtype_count_matrix.txt"
    helmsman <- read.table(helmsman_output_file,
                                    sep="\t",
                    skip = 1)
    rownames(helmsman) <- helmsman$V1
    helmsman <- select(helmsman, -V1)
    colnames(helmsman) <- colnames(cosmic_signatures_v2)
    helmsman <- helmsman[order(row.names(helmsman)), ]


# Fit to signatures
    helmsman_run <- helmsman
    mutational_oppertunities_many <- readRDS("/ssd3/Mar_genome_analysis/mut_sig/trinuc_freq/mutational_oppertunities_all.rds")
    mutopps <- mutational_oppertunities_many[rep(5, nrow(helmsman_run)),]
    signatures_4 <- readRDS("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/sigfit/signatures_4_extracted_full.rds")

    helmsman_run_fit <- fit_signatures(counts = helmsman_run,  
                                    signatures = signatures_4,
                                    opportunities =  mutopps,
                                    iter = 2000, 
                                    warmup = 1000, 
                                    chains = 1, 
                                    seed = 1756)
    helmsman_run_exposures <- retrieve_pars(helmsman_run_fit, par = "exposures")
    # exp_names <- tibble::rownames_to_column(helmsman_run_exposures$mean, var = "rowname") %>%
    #     separate(rowname, c("region",NA,NA,"sample"), sep = "\\.") %>%
    #     select(sample,region)
    exp_mean <- helmsman_run_exposures$mean
    exp_upper <- helmsman_run_exposures$upper_95%>%
        rename(sigS_U=sigS,sig1_U=sig1,sig5_U=sig5,sig40_U=sig40)
    exp_lower <- helmsman_run_exposures$lower_95 %>%
        rename(sigS_L=sigS,sig1_L=sig1,sig5_L=sig5,sig40_L=sig40)
    mutation_counts <- rowSums(helmsman_run) %>% as.data.frame() 
    colnames(mutation_counts) <- "count"
    exposures_run_total <- cbind(mutation_counts,exp_mean,exp_upper,exp_lower)
    exposures_run_total

############################################################### dnds ############################################################
# No longer used
    # bins_to_keep <- c("H_any", "H_USA", "H_PEI", "H_REF", "H_any_hm", "H_USA_hm", "H_PEI_hm", "H_REF_hm", "H_any_ht", "H_USA_ht", "H_PEI_ht", "H_REF_ht",
    #                   "allCanyH", "allCnoH", "allCnoH_hm", "allCnoH_ht", "anyUSAnoPEInoH", "anyPEInoUSAnoH", "somatic", "allUSAnoPEInoH", "allPEInoUSAnoH", "sublineages")
    # # bins_to_keep <- c("allCnoH", "allCnoH_hm", "allCnoH_ht","H_any", "H_USA", "H_PEI", "H_REF", "H_any_ht", "H_USA_ht", "H_PEI_ht", "allUSAnoPEInoH", "allPEInoUSAnoH", "sublineages")

    # count=0
    # for(i in bins_to_keep){
    #     print(i)
    #     dndsout <- dndscv(read.table(paste0(i,".dnds")),
    #                 refdb="/ssd3/Mar_genome_analysis/dnds/dndscv/Mar.3.4.6.p1_snap02_refCDS.rda",
    #                 cv=NULL,
    #                 max_coding_muts_per_sample = 1000000000,
    #                 max_muts_per_gene_per_sample = 1000000000)
    #     #saveRDS(dndsout, paste0(i,".dnds.rds"))   
    #     dndsout$globaldnds[1,] %>% print()
    # # save to file
    #     if(count == 0){
    #         dnds_summary <- data.frame(dndsout$globaldnds[1,], stringsAsFactors = FALSE) %>%
    #         mutate(name = i)
    #     }
    #     if(count > 0){
    #         dnds_value <- data.frame(dndsout$globaldnds[1,], stringsAsFactors = FALSE) %>%
    #             mutate(name = i)
    #         dnds_summary <- rbind(dnds_summary,dnds_value)
    #     }
    #     count=count+1
    #     print("     done", quote = FALSE)	
    # }
    
    # rownames(dnds_summary) <- dnds_summary$name
    # dnds_summary
    # saveRDS(dnds_summary, "/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/bins/dnds/outputs/noLOH.dnds.summary.rds")


# Plots
    # dnds_to_plot <- dnds_summary[c("H_USA", "H_PEI","H_REF","allCanyH", "allCnoH_hm", "allCnoH_ht", "allUSAnoPEInoH", "allPEInoUSAnoH"),] # "H_any",
    # dnds_to_plot$name <- factor(c("H_USA", "H_PEI","H_REF","allCanyH", "allCnoH_hm", "allCnoH_ht", "USA", "PEI"),
    #                             levels = c("H_USA", "H_PEI","H_REF","allCanyH", "allCnoH_hm", "allCnoH_ht", "USA", "PEI"))
    # dnds_to_plot$color <- factor(c("Healthy","Healthy","Healthy","Cancer","Cancer","Cancer","USA","PEI"))
    # dnds_to_plot
    sigS_to_plot <- exposures_run_total[c("H_USA", "H_PEI","H_REF","allCanyH", "allCnoH_hm", "allCnoH_ht", "allUSAnoPEInoH", "allPEInoUSAnoH"),] %>%
        select(count, sigS, sigS_U, sigS_L)
    sigS_to_plot$name <- factor(c("H_USA", "H_PEI","H_REF","allCanyH", "allCnoH_hm", "allCnoH_ht", "USA", "PEI"),
                                levels = c("H_USA", "H_PEI","H_REF","allCanyH", "allCnoH_hm", "allCnoH_ht", "USA", "PEI"))
    sigS_to_plot$color <- factor(c("Healthy","Healthy","Healthy","Cancer","Cancer","Cancer","USA","PEI"))
    sigS_to_plot
    # sigS_to_plot <- exposures_run_total[c("H_any", "allCnoH", "allUSAnoPEInoH", "allPEInoUSAnoH"),]
    # sigS_to_plot$name <- factor(c("Healthy", "allCnoH", "USA", "PEI"), levels = c("Healthy", "allCnoH", "USA", "PEI"))
    # sigS_to_plot$color <- factor(c("Healthy","Cancer","USA","PEI"))
    # sigS_to_plot

#pdf("dnds_sigS_main_plot.pdf", width=2, height=2)
pdf("sigS_main_plot.pdf", width=2, height=2.5)
# ggplot(dnds_to_plot, aes(x=as.factor(name),y=mle,ymin=cilow,ymax=cihigh,color=color))+
#     geom_abline(aes(intercept=dnds_summary["sublineages","mle"], slope=0), color="black", linetype="dashed")+
#     geom_abline(aes(intercept=dnds_summary["allCnoH","mle"], slope=0), color="black", linetype="dashed")+
#     geom_abline(aes(intercept=dnds_summary["H_any","mle"], slope=0), color="black", linetype="dashed")+
#     #geom_pointrange(fatten = 3, size = 3)+  
#     geom_point(size = 3)+  #size = 9
#     geom_errorbar(color="grey20",width=0.25)+
#     ylab("dN/dS")+
#     scale_color_manual(values=c("grey", "black", "red", "blue"))+
#     scale_y_continuous(limits = c(0,1.1),breaks=seq(0,1.0,0.2))+
#     theme_classic() +
#         theme(axis.text=element_text(size=8,face="bold"),
#         axis.title=element_text(size=8,face="bold"),
#         text=element_text(size=8,face="bold"),
#         axis.title.x = element_blank(),
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#         legend.position = "none") #+
#     #ggtitle("dN/dS")
ggplot(sigS_to_plot, aes(x=as.factor(name),y=sigS,ymin=sigS_L,ymax=sigS_U,color=color))+
    geom_abline(aes(intercept=exposures_run_total["sublineages","sigS"], slope=0), color="black", linetype="dashed")+
    #geom_abline(aes(intercept=exposures_run_total["allCnoH","sigS"], slope=0), color="black", linetype="dashed")+
    geom_abline(aes(intercept=exposures_run_total["H_any","sigS"], slope=0), color="black", linetype="dashed")+
    geom_pointrange()+  # fatten = 3, size = 2
    #geom_point(size = 3)+  #size = 9
    #geom_errorbar(color="grey20",width=0.25)+
    ylab("SigS Fraction")+
    scale_color_manual(values=c("grey", "black", "red", "blue"))+
    scale_y_continuous(limits = c(0,0.55),breaks=seq(0,0.5,0.1))+
    #ylim(0,0.6) +
    theme_classic() +
        theme(axis.text=element_text(size=8,face="bold"),
        axis.title=element_text(size=8,face="bold"),
        text=element_text(size=8,face="bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")
dev.off()


# Extract important values
        # allCnoH_dnds <- dnds_summary["allCnoH","mle"]
        # allCnoH_dnds_L <- dnds_summary["allCnoH","cilow"]
        # allCnoH_dnds_U <- dnds_summary["allCnoH","cihigh"] 
        # allCnoH_dnds_Lsd <- allCnoH_dnds - allCnoH_dnds_L
        # allCnoH_dnds_Usd <- allCnoH_dnds_U - allCnoH_dnds
        # allCnoH_dnds_sd <- mean(allCnoH_dnds_Lsd,allCnoH_dnds_Usd)

        # H_dnds <- dnds_summary["H_any","mle"]
        # H_dnds_L <- dnds_summary["H_any","cilow"]
        # H_dnds_U <- dnds_summary["H_any","cihigh"] 
        # H_dnds_Lsd <- H_dnds - H_dnds_L
        # H_dnds_Usd <- H_dnds_U - H_dnds
        # H_dnds_sd <- mean(H_dnds_Lsd,H_dnds_Usd)

        # somatic_dnds <- dnds_summary["sublineages","mle"]
        # somatic_dnds_L <- dnds_summary["sublineages","cilow"]
        # somatic_dnds_U <- dnds_summary["sublineages","cihigh"] 
        # somatic_dnds_Lsd <- somatic_dnds - somatic_dnds_L
        # somatic_dnds_Usd <- somatic_dnds_U - somatic_dnds
        # somatic_dnds_sd <- mean(somatic_dnds_Lsd,somatic_dnds_Usd)

        allCnoH_sigS <- exposures_run_total["allCnoH_ht","sigS"]
        allCnoH_sigS_L <- exposures_run_total["allCnoH_ht","sigS_L"]
        allCnoH_sigS_U <- exposures_run_total["allCnoH_ht","sigS_U"] 
        allCnoH_sigS_Lsd <- allCnoH_sigS - allCnoH_sigS_L
        allCnoH_sigS_Usd <- allCnoH_sigS_U - allCnoH_sigS
        allCnoH_sigS_sd <- mean(allCnoH_sigS_Lsd,allCnoH_sigS_Usd)
        allCnoH_count <- exposures_run_total["allCnoH_ht","count"]
        # allCnoH_sigS <- exposures_run_total["allCnoH","sigS"]
        # allCnoH_sigS_L <- exposures_run_total["allCnoH","sigS_L"]
        # allCnoH_sigS_U <- exposures_run_total["allCnoH","sigS_U"] 
        # allCnoH_sigS_Lsd <- allCnoH_sigS - allCnoH_sigS_L
        # allCnoH_sigS_Usd <- allCnoH_sigS_U - allCnoH_sigS
        # allCnoH_sigS_sd <- mean(allCnoH_sigS_Lsd,allCnoH_sigS_Usd)
        # allCnoH_count <- exposures_run_total["allCnoH","count"]
        H_sigS <- exposures_run_total["H_any","sigS"]
        H_sigS_L <- exposures_run_total["H_any","sigS_L"]
        H_sigS_U <- exposures_run_total["H_any","sigS_U"] 
        H_sigS_Lsd <- H_sigS - H_sigS_L
        H_sigS_Usd <- H_sigS_U - H_sigS
        H_sig_sd <- mean(H_sigS_Lsd,H_sigS_Usd)

        somatic_sigS <- exposures_run_total["sublineages","sigS"]
        somatic_sigS_L <- exposures_run_total["sublineages","sigS_L"]
        somatic_sigS_U <- exposures_run_total["sublineages","sigS_U"] 
        somatic_sigS_Lsd <- somatic_sigS - somatic_sigS_L
        somatic_sigS_Usd <- somatic_sigS_U - somatic_sigS
        somatic_sigS_sd <- mean(somatic_sigS_Lsd,somatic_sigS_Usd)

        

############################################################### calculations ############################################################
#   View all data
    # dnds_summary
    exposures_run_total

# rate stats
setwd("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/pairwise/withLOH")
    sigS_regression <- readRDS("sigS_regression_noLOH_USAonly.rds") 
    sig55_regression <- readRDS("sig55_regression_noLOH_USAonly.rds") 
    summary(sigS_regression)
    sigS_rate <- summary(sigS_regression)$coefficients[2,1] # 495.6361
    sigS_rate_sd <- (confint(sigS_regression)[2,2]-sigS_rate)/2  #  311.7626
    sigS_mutation_count <- summary(sigS_regression)$coefficients[1,1] # 156057
    sig55_rate <- summary(sig55_regression)$coefficients[2,1] # 226.4911
    sig55_rate_sd <- (confint(sig55_regression)[2,2]-sig55_rate)/2  # 68.15301
    sig55_mutation_count <- summary(sig55_regression)$coefficients[1,1] # 72280.41



    postdiv_age_Sest <- sigS_mutation_count/sigS_rate # 314.8621 YEARS FROM DIVERGENCE TO PRESENT ESTIMATE FROM SIGS
    postdiv_age_Sest_L <- sigS_mutation_count/(sigS_rate+sigS_rate_sd*2) # 139.441
    postdiv_age_Sest_U <- sigS_mutation_count/(sigS_rate-sigS_rate_sd*2) # -1220.252
        postdiv_age_Sest_U <- Inf
    postdiv_age_55est <- sig55_mutation_count/sig55_rate # 319.1314 YEARS FROM DIVERGENCE TO PRESENT ESTIMATE FROM SIGS
    postdiv_age_55est_L <- sig55_mutation_count/(sig55_rate+sig55_rate_sd*2) # 199.2309
    postdiv_age_55est_U <- sig55_mutation_count/(sig55_rate-sig55_rate_sd*2) # 801.4676

# Genome size corrections: how much of genome are we excluding with LOH calls
    LOH_genome_size <- read.table("/ssd3/Mar_genome_analysis/LOH/july_2021/output/new/BOTH_SUBLINEAGES_LOH_10_counts.bed") %>% mutate(V4 = V3-V2) %>% .$V4 %>%as.numeric() %>% sum()
    full_genome_size <- read.table("/ssd3/Mar_genome_analysis/genomes/Mar.3.4.6.p1_Q30Q30A.fasta.fai") %>% .$V2 %>%as.numeric() %>% sum() # 1216068291
    full_genome_Gb <- full_genome_size/1000000000 # 1.216068
    nonLOH_genome_size = full_genome_size - LOH_genome_size 
    nonLOH_genome_Gb <- nonLOH_genome_size/1000000000 # 0.978504
    # Rate in mu/Gb/yr
        sigS_rate/nonLOH_genome_Gb # 506.5244
        2*sigS_rate_sd/nonLOH_genome_Gb # 637.223

# estimate age 
    allCnoH_somatic_sigS_est <- allCnoH_sigS - H_sigS # 0.03112695
    allCnoH_somatic_sigS_count_est <- allCnoH_somatic_sigS_est * allCnoH_count # 53350.19
    allCnoH_somatic_sigS_age_est <- allCnoH_somatic_sigS_count_est / sigS_rate # 107.6398 YEARS FROM ORIGIN TO DIVERGENCE ESTIMATE FROM SIGS
    
    allCnoH_somatic_sigS_age_est_L <- allCnoH_somatic_sigS_count_est / (sigS_rate+sigS_rate_sd*2) # 47.6698
    allCnoH_somatic_sigS_age_est_U <- allCnoH_somatic_sigS_count_est / (sigS_rate-sigS_rate_sd*2) # -417.1597
        allCnoH_somatic_sigS_age_est_U <- Inf
# Total ages
    total_age_est <- postdiv_age_Sest + allCnoH_somatic_sigS_age_est # 422.5019
    total_age_est_L <- postdiv_age_Sest_L + allCnoH_somatic_sigS_age_est_L # 187.1108
    total_age_est_U <- postdiv_age_Sest_U + allCnoH_somatic_sigS_age_est_U # -1637.412

# Plot age
    age_summary <- data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c("name", "age", "lower","upper"))))
    age_summary[1,] <- c("Origin", total_age_est, total_age_est_L, total_age_est_U)
    age_summary[2,] <- c("MRCA_S", postdiv_age_Sest, postdiv_age_Sest_L, postdiv_age_Sest_U)
    age_summary[3,] <- c("MRCA_55", postdiv_age_55est, postdiv_age_55est_L, postdiv_age_55est_U)
    # age_summary[,2:4] <- as.numeric(age_summary[1:2,2:4])age_summary[2,4]
    age_summary
    pdf("cancer_age_estimate.pdf", width=2, height=2) # 
    ggplot(age_summary, aes(x=as.factor(name),y=as.numeric(age),ymin=as.numeric(lower),ymax=as.numeric(upper)))+
        geom_abline(aes(intercept=-50, slope=0), color="black", linetype="dashed")+
        geom_pointrange()+  # fatten = 3, size = 2
        #geom_point(size = 3)+  #size = 9
        #geom_errorbar(color="grey20",width=0.25)+
        ylab("Years before present")+
        #scale_y_continuous(limits = c(900,0),breaks=c(800,700,600,500,400,300,200,100,0))+ #as.numeric(age_summary[1,"age"]),as.numeric(age_summary[2,"age"]))
        #ylim(1000,0) +
        theme_classic() +
            theme(axis.text=element_text(size=8,face="bold"),
            axis.title=element_text(size=8,face="bold"),
            text=element_text(size=8,face="bold")) +
        scale_y_reverse(limits = c(900,0)) +
        coord_flip()
    dev.off()

# # age estimate error propagation. NOTE - rate error is almost all the error
#     # upper = high allCnoH_sigS, low H_sigS, low sigS_rate
#     allCnoH_somatic_sigS_exp_est_U <- allCnoH_sigS_U - H_sigS_L #  0.03221753
#     allCnoH_somatic_sigS_count_est_U <- allCnoH_somatic_sigS_exp_est_U * allCnoH_count # 55219.39
#     allCnoH_somatic_sigS_age_est_U <- allCnoH_somatic_sigS_count_est_U / (sigS_rate-sigS_rate_sd*2) # 231.3356 MAX WITH ERROR PROPAGATION
#     # lower = low allCnoH_sigS, high H_sigS, high sigS_rate
#     allCnoH_somatic_sigS_exp_est_L <- allCnoH_sigS_L - H_sigS_U # 0.03005806
#     allCnoH_somatic_sigS_count_est_L <- allCnoH_somatic_sigS_exp_est_L * allCnoH_count # 51518.16
#     allCnoH_somatic_sigS_age_est_L <- allCnoH_somatic_sigS_count_est_L / (sigS_rate+sigS_rate_sd*2) # 84.80453 MIN WITH ERROR PROPAGATION


# mutation count estimate
    allCnoH_sigS_somatic_fraction_estimate = (allCnoH_sigS - H_sigS)/(somatic_sigS - H_sigS) # 0.068126
    allCnoH_sigS_somatic_fraction_estimate * allCnoH_count # 116764.9
    allCnoH_sigS_somatic_fraction_estimate_U = (allCnoH_sigS_U - H_sigS_L)/(somatic_sigS_L - H_sigS_L) #  0.07083365
    allCnoH_sigS_somatic_fraction_estimate_U * allCnoH_count # 121405.7
    allCnoH_sigS_somatic_fraction_estimate_L = (allCnoH_sigS_L - H_sigS_U)/(somatic_sigS_U - H_sigS_U) # 0.06544987
    allCnoH_sigS_somatic_fraction_estimate_L * allCnoH_count # 112178.1

    # allCnoH_dnds_somatic_fraction_estimate = (allCnoH_dnds - H_dnds)/(somatic_dnds - H_dnds) # 0.1301147
    # allCnoH_dnds_somatic_fraction_estimate * allCnoH_count # 289705.8
    # allCnoH_dnds_somatic_fraction_estimate_U = (allCnoH_dnds_U - H_dnds_L)/(somatic_dnds_L - H_dnds_L) # 0.169512
    # allCnoH_dnds_somatic_fraction_estimate_U * allCnoH_count # 377425.7
    # allCnoH_dnds_somatic_fraction_estimate_L = (allCnoH_dnds_L - H_dnds_U)/(somatic_dnds_U - H_dnds_U) # 0.09611125
    # allCnoH_dnds_somatic_fraction_estimate_L * allCnoH_count #  213995.7

# Total counts and per Mb estimate
    total_nonLOH_mu_est <- (allCnoH_sigS_somatic_fraction_estimate * allCnoH_count) + exposures_run_total["sublineages","count"]/2 # 437205.4
    total_nonLOH_mu_est/nonLOH_genome_size*1000000 #  446.81
    total_nonLOH_mu_est <- (allCnoH_sigS_somatic_fraction_estimate * allCnoH_count) + exposures_run_total["allPEInoUSAnoH","count"] # 431986.9
    total_nonLOH_mu_est/nonLOH_genome_size*1000000 #  441.4769
    total_nonLOH_mu_est <- (allCnoH_sigS_somatic_fraction_estimate * allCnoH_count) + exposures_run_total["allUSAnoPEInoH","count"] # 442423.9
    total_nonLOH_mu_est/nonLOH_genome_size*1000000 #  452.1432