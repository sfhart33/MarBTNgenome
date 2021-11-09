library(tidyverse)
library(ape)
library(sigfit)
library(lsa)
data(cosmic_signatures_v2)
data(cosmic_signatures_v3)
    library(geiger)
    library(nlme)
    library(phytools)



# Load data
    setwd("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/pairwise/noLOH")
    helmsman_output_file <- paste0("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/pairwise/noLOH/helmsman_div_output/subtype_count_matrix.txt")
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
    plot_exposures(mcmc_samples = helmsman_div_fit, pdf_path = "post_div_exposures.pdf")
    plot_reconstruction(mcmc_samples = helmsman_div_fit, pdf_path = "post_div_reconstruction.pdf")


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
    exposures_total_div1 <- filter(exposures_total_div,div=="div1")
    # exposures_total_div2 <- filter(exposures_total_div,div=="div2")
    # exposures_total_div2 <- exposures_total_div2[4:8,]
    # exposures_total_div_usa
    write.table(exposures_total_div,file = "noLOH_mutation_sigs_count_postdiv.txt", quote=FALSE, sep = '\t')

# make plots

pdf("signatures_time_regressions_postdiv.pdf")
#for(file in c("div1","div2")){
file="div"
    for(sig in c("sigS","sig1","sig5","sig40")){
        plot1 <- ggplot(get(paste0("exposures_total_",file)), aes(x = dates, y = get(sig), ymin=get(paste0(sig,"_L")), ymax=get(paste0(sig,"_U")))) + 
            geom_pointrange()+
            stat_smooth(method = "lm") + #, fullrange = T) +
            # xlim(-600,0)+
            # ylim(-100000,350000)+
            xlab("Years before present")+
            ylab(paste(sig,"mutations"))+
            theme_classic() +
                theme(axis.text=element_text(size=12,face="bold"),
                axis.title=element_text(size=16,face="bold"),
                text=element_text(size=18,face="bold")) +
            ggtitle(paste("Regression of",sig,"mutations,",file))
        print(plot1)
    }
            plot1 <- ggplot(get(paste0("exposures_total_",file)), aes(x = dates, y = (sigS+sig1+sig5+sig40), ymin=(sigS_L+sig1_L+sig5_L+sig40_L), ymax=(sigS_U+sig1_U+sig5_U+sig40_U))) + 
                geom_pointrange()+
                stat_smooth(method = "lm") + #, fullrange = T) +
                # xlim(-600,0)+
                # ylim(-100000,350000)+
                xlab("Years before present")+
                xlab(paste("All mutations"))+
                theme_classic() +
                    theme(axis.text=element_text(size=12,face="bold"),
                    axis.title=element_text(size=16,face="bold"),
                    text=element_text(size=18,face="bold")) +
                ggtitle(paste("Regression of all mutations",file))
        print(plot1)
#}
dev.off()

sigS_regression <- lm(sigS ~ dates, data=exposures_total_div1)
sig1_regression <- lm(sig1 ~ dates, data=exposures_total_div1)
sig5_regression <- lm(sig5 ~ dates, data=exposures_total_div1)
sig40_regression <- lm(sig40 ~ dates, data=exposures_total_div1)
saveRDS(sigS_regression, "sigS_regression.rds")

summary(sigS_regression)
    # Coefficients:
    #             Estimate Std. Error t value Pr(>|t|)
    # (Intercept) 155902.92     727.82 214.205 6.99e-13 ***
    # dates          420.57      92.88   4.528  0.00398 **

summary(sig1_regression)
summary(sig5_regression)
summary(sig40_regression)

##################################
# NEW incorperate Phylogenetic Generalised Least Squares to make sure effect is not phylogeny-based
    library(ape)
    library(geiger)
    library(nlme)
    library(phytools)
    BTN.tree <- readRDS("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/pairwise/BTN.tree.phy.rds") 
    BTN.tree$edge.length[16] <- 0 # set root length to 0
    BTN.tree <- root(BTN.tree, outgroup = "genome", resolve.root = TRUE)

    # Set names
    rownames(exposures_total_div1) <- paste0("C_",exposures_total_div1$sample)
    name.check(BTN.tree, exposures_total_div1)
    BTN.tree<-drop.tip(BTN.tree, name.check(BTN.tree, exposures_total_div1)$tree_not_data)
    name.check(BTN.tree, exposures_total_div1)

    # PIC method
        sigS_reg <- exposures_total_div1[, "sigS"]
        dates_reg <- exposures_total_div1[, "dates"]
        names(sigS_reg) <- rownames(exposures_total_div1)
        names(dates_reg) <- rownames(exposures_total_div1)
        sigS_Pic2 <- pic(sigS_reg, BTN.tree)
        dates_Pic2 <- pic(dates_reg, BTN.tree)
        picModel2 <- lm(sigS_Pic2 ~ dates_Pic2)
        summary(picModel2)
            # Coefficients:
            #             Estimate Std. Error t value Pr(>|t|)
            # (Intercept)    7.167      2.411   2.973   0.0311 *
            # dates_Pic2   332.907    118.399   2.812   0.0375 *

    pdf("BTN.tree.datingcorrection.pdf")
        plot.phylo(BTN.tree)
    dev.off()

    # PGLS method
        glsModel <- gls(sigS ~ dates, data = exposures_total_div1, method = "ML")
        summary(glsModel) # significant (same as lm regression done previously)
            # Coefficients:
            #                 Value Std.Error   t-value p-value
            # (Intercept) 155902.92  727.8208 214.20509   0.000
            # dates          420.57   92.8832   4.52793   0.004
            # 420.57+92.8832*2 = 606.3364
        pglsModel <- gls(sigS ~ dates, correlation = corBrownian(phy = BTN.tree), data = exposures_total_div1, method = "ML")
        summary(pglsModel) # Not significant
            # Coefficients:
            #                 Value Std.Error  t-value p-value
            # (Intercept) 154610.33  4823.482 32.05368  0.0000
            # dates          256.72   174.964  1.46726  0.1927
            # 256.72+174.964*2 = 606.648
        # pglsModel <- gls(sigS ~ dates, correlation = corPagel(1, phy = BTN.tree, fixed = FALSE), data = exposures_total_div1, method = "ML")
        # summary(pglsModel) # Doesn't work
        # pglsModel <- gls(sigS ~ dates, correlation = corMartins(1, BTN.tree), data = exposures_total_div1, method = "ML")
        # summary(pglsModel) # Same output as without a correction - I think samples are getting scrambled
