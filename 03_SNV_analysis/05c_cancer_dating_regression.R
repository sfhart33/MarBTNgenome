# Original file here: C:\Users\shart\Metzger Lab Dropbox\Sam_data\Mya_genome\SNPs\somatypus_output_dating2withLOH.r

# Back in R: sigfit
    library(tidyverse)
    library(ape)
    library(sigfit)
    library(lsa)
    data(cosmic_signatures_v2)
    data(cosmic_signatures_v3)
    library(ape)
    library(geiger)
    library(nlme)
    library(phytools)
    library(gridExtra)

##############################
# Data WITH and without LOH
for(run in c("noLOH")){ #"withLOH",
    # Load data
        setwd(paste0("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/pairwise/",run))
        helmsman_output_file <- paste0("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/pairwise/",run,"/helmsman_div_output/subtype_count_matrix.txt")
        div_counts <- readRDS(file = "post_div_counts.rds")
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
        mutopps <- mutational_oppertunities_many[5,] # for full genome
        # Latest full extraction finalized
        signatures4 <- readRDS(file = "/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/sigfit/signatures_4_extracted_full.rds")


    # FIT HELMSMAN DATA TO SIGNATURES
        helmsman_div_fit <- fit_signatures(counts = helmsman,  
                                        signatures = signatures4,
                                        opportunities =  mutational_oppertunities_many[rep(5, nrow(helmsman)),],
                                        iter = 2000, 
                                        warmup = 1000, 
                                        chains = 1, 
                                        seed = 1756)
        helmsman_div_exposures <- retrieve_pars(helmsman_div_fit, 
                                par = "exposures")
        #plot_exposures(mcmc_samples = helmsman_div_fit, pdf_path = "post_div_exposures.pdf")
        #plot_reconstruction(mcmc_samples = helmsman_div_fit, pdf_path = "post_div_reconstruction.pdf")


    # extract data needed to do aging analysis in R rather than excel
        exp_names_div <- tibble::rownames_to_column(helmsman_div_exposures$mean, var = "rowname") %>%
            separate(rowname, c("sample","div"), sep = "_") %>%
            select(sample,div)
        exp_mean_div <- helmsman_div_exposures$mean
        exp_upper_div <- helmsman_div_exposures$upper_95 %>%
            rename(sigS_U=sigS,sig1_U=sig1,sig5_U=sig5,sig40_U=sig40)
        exp_lower_div <- helmsman_div_exposures$lower_95 %>%
            rename(sigS_L=sigS,sig1_L=sig1,sig5_L=sig5,sig40_L=sig40)
        exposures_total_div <- cbind(exp_names_div, exp_mean_div,exp_upper_div,exp_lower_div)
        exposures_total_div <- tibble::rownames_to_column(exposures_total_div, var = "rowname") %>%
            separate(rowname, c("sample","div"), sep = "_")
        exposures_total_div[,3:14] <- exposures_total_div[,3:14]*div_counts
        exposures_total_div <- arrange(exposures_total_div,div,sample)
        # sampling dates (number is years before 8/24/21)
        dates <- c(-10.83,-11.28,-11.28,-0.40,-0.98,-0.98,-7.94,-7.37)
        exposures_total_div$dates <- dates # c(dates,dates)
        exposures_total_div
        exposures_total_div$sub <- c("PEI","PEI","PEI","USA","USA","USA","USA","USA")
        exposures_total_div1 <- filter(exposures_total_div,div=="div1") %>%
            mutate(LOH = run,
                   sigALL = sigS + sig1 + sig5 + sig40,
                   sigALL_U = sigS + sig1 + sig5 + sig40,
                   sigALL_L = sigS + sig1 + sig5 + sig40
            )
        assign(paste0("exp_",run), exposures_total_div1)
        #write.table(exposures_total_div,file = "mutation_sigs_count_postdiv.txt", quote=FALSE, sep = '\t')
}

#######################################



# how much of genome are we excluding with LOH calls
    LOH_genome_size <- read.table("/ssd3/Mar_genome_analysis/LOH/july_2021/output/new/BOTH_SUBLINEAGES_LOH_10_counts.bed") %>% mutate(V4 = V3-V2) %>% .$V4 %>%as.numeric() %>% sum()
    full_genome_size <- read.table("/ssd3/Mar_genome_analysis/genomes/Mar.3.4.6.p1_Q30Q30A.fasta.fai") %>% .$V2 %>%as.numeric() %>% sum()
    nonLOH_genome_size = full_genome_size - LOH_genome_size 
    # perGbcorrection = nonLOH_genome_size/1000000000
    full_genome_correction <- full_genome_size/nonLOH_genome_size




# New 2023 revision supp figure
    setwd(paste0("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/pairwise/","withLOH"))
    pdf("new_supp_fig_regression.pdf",width=3,height=3)
      perGbcorrection = nonLOH_genome_size/1000000000
        
          data_colors <- c("red", "blue")
          datasubset <- get(paste0("exp_",run)) %>% 
            filter(sub == "USA")
        for(sig in c("sigS","sig1","sig5","sig40","sigALL")){
          sig_regression <- lm(get(sig) ~ dates, data=datasubset)
          #assign(paste("regression", dataset, run, sig, sep = "_"), regression)
          sig_rsquared <- summary(sig_regression)$r.squared
          sig_p_value <- summary(sig_regression)$coefficients[2,4]
          sig_intercept <- summary(sig_regression)$coefficients[1,1]
          sig_rate <- summary(sig_regression)$coefficients[2,1]
          #sig_rate_sd <- summary(sig_regression)$coefficients[2,2]
          sig_rate_ci <- confint(sig_regression)[2,2]-sig_rate 
          div_estimate <- round(sig_intercept/sig_rate)
          div_estimate_lower95 <- round(sig_intercept/(sig_rate+sig_rate_ci))
          if((-sig_rate)<sig_rate_ci & sig_rate<0){div_estimate_lower95 <- "-Inf"}
          div_estimate_upper95 <- round(sig_intercept/(sig_rate-sig_rate_ci))
          if(sig_rate<sig_rate_ci & sig_rate>0){div_estimate_upper95 <- "Inf"}
          plot1 <- ggplot(datasubset, aes(x = jitter(dates+2021.73), y = get(sig)/perGbcorrection/1000, ymin=get(paste0(sig, "_L"))/perGbcorrection/1000, ymax=get(paste0(sig, "_U"))/perGbcorrection/1000)) + 
            scale_color_manual(values=data_colors)+
            stat_smooth(method = "lm",color="blue",fullrange = T) + 
            geom_pointrange(data = exp_noLOH, aes(x = jitter(dates+2021.73), y = get(sig)/perGbcorrection/1000, ymin=get(paste0(sig, "_L"))/perGbcorrection/1000, ymax=get(paste0(sig, "_U"))/perGbcorrection/1000, color=sub))+ # , size=1, fatten = 5
            scale_x_continuous(breaks=seq(2010,2022,2))+
            xlab("Year sampled")+
            ylab("Mutations per Mb\n(since USA-PEI split)")+   
            scale_fill_manual(values=data_colors)+
            ggtitle(paste0(#"Regression: ", run, " ", sig, "\n",
              "R^2 = ",round(sig_rsquared,3),", p = ",round(sig_p_value,3), "\n",
              #round(sig_rate/perGbcorrection)," +/- ",round(2*sig_rate_sd/perGbcorrection)," mu/Gb/year", "\n",
              round(sig_rate/perGbcorrection)," mu/Gb/year (95%CI: ",round((sig_rate-sig_rate_ci)/perGbcorrection)," - ",round((sig_rate+sig_rate_ci)/perGbcorrection), ")\n",
              "Split: ",div_estimate," years ago (95%CI: ",div_estimate_lower95," - ",div_estimate_upper95,")")
            ) +
            theme_classic() +
            theme(axis.text=element_text(size=8,face="bold"),
                  axis.title=element_text(size=8,face="bold"),
                  #axis.title=element_blank(),
                  text=element_text(size=8,face="bold"),
                  #title=element_text(size=8, face='bold'),
                  plot.title = element_text(size=8, face='bold'),
                  legend.position = "none")
          print(plot1) 
          assign(paste0(sig, "_plot"), plot1)
        }
        # grid.arrange(grobs=gs, ncol=4, 
        #        top="top label", bottom="bottom\nlabel", 
        #        left="left label", right="right label")
        # grid.rect(gp=gpar(fill=NA))
    dev.off()
