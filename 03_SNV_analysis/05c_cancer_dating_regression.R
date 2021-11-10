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
    # (Intercept) 155919.1      722.5 215.818 6.68e-13 ***
    # dates          423.1       92.2   4.589  0.00374 **

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
            # (Intercept)    7.139      2.391   2.986   0.0306 *
            # dates_Pic2   338.472    117.412   2.883   0.0345 *

            # 338.472+117.412*2 =  573.296

    pdf("BTN.tree.datingcorrection.pdf")
        plot.phylo(BTN.tree)
    dev.off()

    # PGLS method
        glsModel <- gls(sigS ~ dates, data = exposures_total_div1, method = "ML")
        summary(glsModel) # significant (same as lm regression done previously)
            # Coefficients:
            #             Value Std.Error   t-value p-value
            # (Intercept) 155919.1  722.4572 215.81781  0.0000
            # dates          423.1   92.1987   4.58895  0.0037

            # 423.1+92.1987*2 = 607.5
        pglsModel <- gls(sigS ~ dates, correlation = corBrownian(phy = BTN.tree), data = exposures_total_div1, method = "ML")
        summary(pglsModel) # Not significant
            # Coefficients:
            #                 Value Std.Error  t-value p-value
            # (Intercept) 154656.32  4797.512 32.23678  0.0000
            # dates          262.47   174.022  1.50827  0.1822

            # 262.47+174.022*2 = 610.5

        # pglsModel <- gls(sigS ~ dates, correlation = corPagel(1, phy = BTN.tree, fixed = FALSE), data = exposures_total_div1, method = "ML")
        # summary(pglsModel) # Doesn't work
        # pglsModel <- gls(sigS ~ dates, correlation = corMartins(1, BTN.tree), data = exposures_total_div1, method = "ML")
        # summary(pglsModel) # Same output as without a correction - I think samples are getting scrambled

################################# Old notes

# how much of genome are we excluding with LOH calls
    LOH_genome_size <- read.table("/ssd3/Mar_genome_analysis/LOH/july_2021/output/new/BOTH_SUBLINEAGES_LOH_10_counts.bed") %>% mutate(V4 = V3-V2) %>% .$V4 %>%as.numeric() %>% sum()
    full_genome_size <- read.table("/ssd3/Mar_genome_analysis/genomes/Mar.3.4.6.p1_Q30Q30A.fasta.fai") %>% .$V2 %>%as.numeric() %>% sum()
    nonLOH_genome_size = full_genome_size - LOH_genome_size 
    perGbcorrection = nonLOH_genome_size/1000000000
    full_genome_correction <- full_genome_size/nonLOH_genome_size

# More final sigS regression plot
    exposures_total_div1$sub <- c("PEI","PEI","PEI","USA","USA","USA","USA","USA")
    sigS_rsquared <- summary(sigS_regression)$r.squared
    sigS_p_value <- summary(sigS_regression)$coefficients[2,4]
    sigS_intercept <- summary(sigS_regression)$coefficients[1,1]
    sigS_rate <- summary(sigS_regression)$coefficients[2,1]
    sigS_rate_sd <- summary(sigS_regression)$coefficients[2,2]
    div_estimate <- sigS_intercept/sigS_rate
    div_estimate_lower95 <- sigS_intercept/(sigS_rate+2*sigS_rate_sd)
    div_estimate_upper95 <- sigS_intercept/(sigS_rate-2*sigS_rate_sd)
    origin_estimate_dnds <- (allCnoH_sigS_somatic_count_estimate+sigS_intercept)/sigS_rate
    origin_estimate_dnds <- (allCnoH_sigS_somatic_count_estimate+sigS_intercept)/sigS_rate
    print(paste("Estimate of PEI-USA split:",round(div_estimate),"years ago (95% CI:",round(div_estimate_lower95),"-",round(div_estimate_upper95),")"))
    # print(paste("Estimate of origin:",round(origin_estimate),"years ago"))


    pdf("sigS_time_regression3.pdf",width=3,height=3)
    #dates changed to year
    ggplot(exposures_total_div1, aes(x = jitter(dates+2021.73), y = sigS/perGbcorrection/1000, ymin=sigS_L/perGbcorrection/1000, ymax=sigS_U/perGbcorrection/1000)) + 
        scale_color_manual(values=c("red", "blue"))+
        stat_smooth(method = "lm",color="black",fullrange = T) + 
        geom_pointrange(aes(color=sub))+ # , size=1, fatten = 5
        #geom_text(aes(x=2017.5, y=157900/1000, label = paste(round(sigS_rate/perGbcorrection),"+/-",round(2*sigS_rate_sd/perGbcorrection),"mu/Gb/year"), angle = 32), size=3)+
        #geom_text(aes(x=2018, y=154000/1000, label = paste("R-squared:",round(sigS_rsquared,3),"\n p-value:",round(sigS_p_value,3))), size=3)+
        #geom_text(aes(x=2014, y=160000/1000, label = paste("Estimate of PEI-USA split: \n",round(div_estimate),"years ago \n (95% CI:",round(div_estimate_lower95),"-",round(div_estimate_upper95),"years)")), size=3)+
        #geom_text(aes(x=2014, y=160000/1000,
        #    label = paste(round(sigS_rate/perGbcorrection),"+/-",round(2*sigS_rate_sd/perGbcorrection),"mu/Gb/year","\n","R-squared:",round(sigS_rsquared,3),"\n", "p-value:",round(sigS_p_value,3))), size=3)+
        #xlim(2010,2022)+
        scale_x_continuous(breaks=seq(2010,2022,2))+
        scale_y_continuous(breaks=seq(152,162,2))+
        xlab("Year sampled")+
        ylab("SigS mutations per Mb\n(since USA-PEI split)")+   
        scale_fill_manual(values=c("red", "blue"))+
        theme_classic() +
            theme(axis.text=element_text(size=8,face="bold"),
            axis.title=element_text(size=8,face="bold"),
            text=element_text(size=8,face="bold"),
            legend.position = "none") #+
        #ggtitle(paste("Estimate of PEI-USA split:",round(div_estimate),"years ago (95% CI:",round(div_estimate_lower95),"-",round(div_estimate_upper95),"years)"))
    dev.off()

