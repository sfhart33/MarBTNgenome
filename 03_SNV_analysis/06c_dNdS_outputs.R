############################################################################################
# after creating dnds files for all and running dndscv
# /ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/bins/dnds/README
    library(tidyverse)
    library(gridExtra)
    library(dndscv)
    #library(tidyverse)
    setwd("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/bins/dnds/outputs")
    dnds_files <- list.files(pattern = ".dnds.rds")

    count=0
    for(dnds_file in dnds_files){
        dnds_name <- str_split(dnds_file,".dnds")[[1]][1]
        print(paste("Loading",dnds_name), quote = FALSE)
        dnds_output <- readRDS(dnds_file)
        print(dnds_output$globaldnds[1,])
        if(count == 0){
            dnds_summary <- data.frame(dnds_output$globaldnds[1,], stringsAsFactors = FALSE) %>%
            mutate(name = dnds_name)
        }
        if(count > 0){
            dnds_value <- data.frame(dnds_output$globaldnds[1,], stringsAsFactors = FALSE) %>%
                mutate(name = dnds_name)
            dnds_summary <- rbind(dnds_summary,dnds_value)
        }
        count=count+1
        print("     done", quote = FALSE)	
    }

    print(dnds_summary)
    saveRDS(dnds_summary, "dnds_global_summary.rds")
    # dnds_summary <- readRDS("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/bins/dnds/outputs/dnds_global_summary.rds")
    dnds_summary_raw <- dnds_summary
    dnds_healthy <- dnds_summary[5:12,]
    dnds_cancer_real <- dnds_summary[c(1:4,27:36),]

    healthy_dnds <- readRDS("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/bins/dnds/outputs/anyH.dnds.rds")
    healthy_dnds$globaldnds

pdf("global_dnds_raw.pdf", width=30, height=10)
ggplot(dnds_summary, aes(x=as.factor(name),y=mle,ymin=cilow,ymax=cihigh))+
    geom_pointrange()+  
    geom_abline(aes(intercept=1, slope=0), color="black", linetype="dashed")+
    geom_abline(aes(intercept=dnds_summary[29,2], slope=0), color="black", linetype="dashed")+
    xlab("SNV bin")+
    ylab("dNdS (corrected: 1 = neutral)")+
    #ylim(0.4,1.6)+
    theme_classic() +
        theme(axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=16,face="bold"),
        text=element_text(size=18,face="bold")) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ggtitle("dNdS for bins (no correction for LOH regions)")
ggplot(dnds_summary, aes(x=as.factor(name),y=mle,ymin=cilow,ymax=cihigh))+
    geom_pointrange()+  
    geom_abline(aes(intercept=1, slope=0), color="black", linetype="dashed")+
    geom_abline(aes(intercept=dnds_summary[29,2], slope=0), color="black", linetype="dashed")+
    xlab("SNV bin")+
    ylab("dNdS (corrected: 1 = neutral)")+
    ylim(0.4,1.2)+
    theme_classic() +
        theme(axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=16,face="bold"),
        text=element_text(size=18,face="bold")) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ggtitle("dNdS for bins (no correction for LOH regions)")
ggplot(dnds_healthy, aes(x=as.factor(name),y=mle,ymin=cilow,ymax=cihigh))+
    geom_pointrange()+  
    geom_abline(aes(intercept=1, slope=0), color="black", linetype="dashed")+
    geom_abline(aes(intercept=dnds_summary[5,2], slope=0), color="black", linetype="dashed")+
    xlab("SNV bin")+
    ylab("dNdS (corrected: 1 = neutral)")+
    #ylim(0.4,1.6)+
    theme_classic() +
        theme(axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=16,face="bold"),
        text=element_text(size=18,face="bold")) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ggtitle("Healthy SNVs: dNdS for bins (no correction for LOH regions)")
ggplot(dnds_cancer_real, aes(x=as.factor(name),y=mle,ymin=cilow,ymax=cihigh))+
    geom_pointrange()+  
    geom_abline(aes(intercept=1, slope=0), color="black", linetype="dashed")+
    geom_abline(aes(intercept=dnds_summary[29,2], slope=0), color="black", linetype="dashed")+
    xlab("SNV bin")+
    ylab("dNdS (corrected: 1 = neutral)")+
    #ylim(0.4,1.6)+
    theme_classic() +
        theme(axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=16,face="bold"),
        text=element_text(size=18,face="bold")) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ggtitle("Cancer SNVs: dNdS for bins (no correction for LOH regions)")
dev.off()

# looking at somtic bins and incorperating LOH
    setwd("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/bins/dnds/somatic")
    dnds_files <- list.files(pattern = ".dnds.rds")

    count=0
    for(dnds_file in dnds_files){
        dnds_name <- str_split(dnds_file,".dnds")[[1]][1]
        print(paste("Loading",dnds_name), quote = FALSE)
        dnds_output <- readRDS(dnds_file)
        print(dnds_output$globaldnds[1,])
        if(count == 0){
            dnds_summary <- data.frame(dnds_output$globaldnds[1,], stringsAsFactors = FALSE) %>%
            mutate(name = dnds_name)
        }
        if(count > 0){
            dnds_value <- data.frame(dnds_output$globaldnds[1,], stringsAsFactors = FALSE) %>%
                mutate(name = dnds_name)
            dnds_summary <- rbind(dnds_summary,dnds_value)
        }
        count=count+1
        print("     done", quote = FALSE)	
    }
    
    # dnds_summary <- dnds_summary[c(-4,-8,-12,-16),]
    print(dnds_summary)
    saveRDS(dnds_summary, "dnds_global_summary_somatic.rds")
    #readRDS("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/bins/dnds/somatic/dnds_global_summary_somatic.rds")

# Final plot
    usa_somatic <- readRDS("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/bins/dnds/LOH/allUSAnoPEInoH_regular_10_nonLOH.dnds.rds")
    pei_somatic <- readRDS("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/bins/dnds/LOH/allPEInoUSAnoH_regular_10_nonLOH.dnds.rds")
    peiusa_somatic <- readRDS("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/bins/dnds/somatic/sublineages_regular_10_nonLOH_uniqmult.dnds.rds")
    usa_somatic2 <- data.frame(usa_somatic$globaldnds[1,], stringsAsFactors = FALSE) %>%
                mutate(name = "onlyUSAnoLOH")    
    pei_somatic2 <- data.frame(pei_somatic$globaldnds[1,], stringsAsFactors = FALSE) %>%
                mutate(name = "onlyPEInoLOH")
    peiusa_somatic2 <- data.frame(peiusa_somatic$globaldnds[1,], stringsAsFactors = FALSE) %>%
                mutate(name = "USAorPEInoLOH")
    raw_global <- readRDS("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/bins/dnds/outputs/dnds_global_summary.rds")
    raw_global <- rbind(raw_global,usa_somatic2,pei_somatic2,peiusa_somatic2)
    rownames(raw_global) <- raw_global$name
    raw_global_plot <- raw_global[c("H_RUP", "anyH","H_P","H_U",
                                            "allCanyH","allCnoH","allPEInoUSAnoH","allUSAnoPEInoH","onlyPEInoLOH","onlyUSAnoLOH", "USAorPEInoLOH"),]
    raw_global_plot$name <- factor(c("All_H", "Any_H","PEI_H_only","USA_H_only", "All_C_any_H","All_C_no_H","PEI_only","USA_only","PEI_only_nonLOH","USA_only_nonLOH","USA/PEI_nonLOH"),
                                levels = c("All_H", "Any_H","PEI_H_only","USA_H_only", "All_C_any_H","All_C_no_H","PEI_only","USA_only","PEI_only_nonLOH","USA_only_nonLOH","USA/PEI_nonLOH"))
    raw_global_plot$color <- c("Healthy","Healthy","Healthy","Healthy",
                                                "Cancer","Cancer","PEI","USA","PEI","USA", "Cancer")    

pdf("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/bins/dnds/LOH/global_dnds_condensed.pdf", width=7.5, height=10)
ggplot(raw_global_plot, aes(x=as.factor(name),y=mle,ymin=cilow,ymax=cihigh,color=color))+
    #geom_pointrange(fatten = 3, size = 3)+  
    geom_point(size = 9)+  
    geom_errorbar(color="grey20",width=0.5)+
    geom_abline(aes(intercept=1, slope=0), color="black", linetype="dashed")+
    #geom_abline(aes(intercept=raw_global_plot["allCnoH","mle"], slope=0), color="black", linetype="dashed")+
    #geom_abline(aes(intercept=mean(c(raw_global_plot["H_P","mle"],raw_global_plot["H_U","mle"])), slope=0), color="black", linetype="dashed")+
    ylab("dN/dS (corrected: 1 = neutral)")+
    scale_color_manual(values=c("grey", "black", "red", "blue"))+
    theme_classic() +
        theme(axis.text=element_text(size=24,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        text=element_text(size=24,face="bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none") +
    ggtitle("dN/dS")
dev.off()

#####################################
# start looking at individual genes
    setwd("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/bins/dnds/somatic")
    dnds_best_somatic <- readRDS("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/bins/dnds/somatic/sublineages_regular_10_nonLOH_uniqmult.dnds.rds") # readRDS("sublineages_regular_25_nonLOH_uniqmult.dnds.rds")
    usa_somatic <- readRDS("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/bins/dnds/LOH/allUSAnoPEInoH_regular_10_nonLOH.dnds.rds")
    pei_somatic <- readRDS("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/bins/dnds/LOH/allPEInoUSAnoH_regular_10_nonLOH.dnds.rds")
    sel_cv <- dnds_best_somatic$sel_cv
    sel_cv_usa <- usa_somatic$sel_cv
    sel_cv_pei <- pei_somatic$sel_cv
    filter(sel_cv, qmis_cv <0.1) %>% arrange(qmis_cv) # 9
    filter(sel_cv, qmis_cv <0.1, wmis_cv >1) %>% arrange(qmis_cv) # 7 positive selection
    filter(sel_cv, qmis_cv <0.05, wmis_cv >1) %>% arrange(qmis_cv) # 5 positive selection
    filter(sel_cv, qmis_cv <0.1, wmis_cv <1, n_syn >0) %>% arrange(qmis_cv) # 1 negative selection
    #filter(sel_cv, qmis_cv <0.9, wmis_cv <1, n_syn >0) %>% arrange(qmis_cv) %>% head(20)
    poshits <- filter(sel_cv, qmis_cv <0.05, wmis_cv >1) %>% arrange(qmis_cv)
    neghits <- filter(sel_cv, qmis_cv <0.05, wmis_cv <1, n_syn >0) %>% arrange(qmis_cv) # none
    rbind(poshits,neghits) %>% 
        mutate(name = paste0(gene_name,"-mRNA-1")) %>%
        select(name) %>%
        write.table(file = "pos.hits.txt", row.names=FALSE, col.names=FALSE, quote=FALSE, sep = '\t')	# blast these hits
    filter(sel_cv, qmis_cv <0.1, wmis_cv >1) %>% arrange(qmis_cv) %>% select(gene_name,n_syn,n_mis)
    filter(sel_cv_usa, qmis_cv <0.1, wmis_cv >1) %>% arrange(qmis_cv) %>% select(gene_name,n_syn,n_mis)
    filter(sel_cv_pei, qmis_cv <0.1, wmis_cv >1) %>% arrange(qmis_cv) %>% select(gene_name,n_syn,n_mis)
    filter(sel_cv, qmis_cv <0.5, wmis_cv <1, n_syn >0) %>% arrange(qmis_cv) %>% select(gene_name,n_syn,n_mis)
    filter(sel_cv_usa, qmis_cv <0.5, wmis_cv <1, n_syn >0) %>% arrange(qmis_cv) %>% select(gene_name,n_syn,n_mis)
    filter(sel_cv_pei, qmis_cv <0.5, wmis_cv <1, n_syn >0) %>% arrange(qmis_cv) %>% select(gene_name,n_syn,n_mis)
##################

    healthy_cv <- healthy_dnds$sel_cv
    healthy_pos <- filter(healthy_cv, qmis_cv <0.1, wmis_cv >1) %>% arrange(qmis_cv)
    filter(healthy_cv, gene_name == "maker-Mar.3.4.6.p1_scaffold3-snap-gene-789.7") # is under negitive selection
    
    filter(sel_cv, qmis_cv <0.05, wmis_cv >1) %>% arrange(qmis_cv) %>% select(gene_name)
    somatic_pos_H <- filter(healthy_cv, gene_name == "maker-Mar.3.4.6.p1_scaffold15-snap-gene-127.5"|
                        gene_name == "maker-Mar.3.4.6.p1_scaffold9-snap-gene-700.19"|
                        gene_name == "maker-Mar.3.4.6.p1_scaffold11-snap-gene-477.12"|
                        gene_name == "maker-Mar.3.4.6.p1_scaffold4-snap-gene-785.0"|
                        gene_name == "maker-Mar.3.4.6.p1_scaffold6-snap-gene-524.0" #|
                        #gene_name == "maker-Mar.3.4.6.p1_scaffold6-snap-gene-524.3"|
                        #gene_name == "maker-Mar.3.4.6.p1_scaffold6-snap-gene-345.0"
                        ) %>% select(gene_name,n_syn,n_mis,wmis_cv,qmis_cv) 
    #somatic_neg_H <- filter(healthy_cv, gene_name == "maker-Mar.3.4.6.p1_scaffold3-snap-gene-789.7") %>% select(gene_name,n_syn,n_mis,wmis_cv,qmis_cv) 

                        maker-Mar.3.4.6.p1_scaffold3-snap-gene-789.7 
    somatic_pos <- filter(sel_cv, qmis_cv <0.05, wmis_cv >1) %>% arrange(qmis_cv) %>% select(gene_name,n_syn,n_mis,wmis_cv,qmis_cv)
    left_join(somatic_pos,somatic_pos_H,by="gene_name", suffix = c("_C","_H"))
    #     gene_name	n_syn_C	n_mis_C	wmis_cv_C	qmis_cv_C	n_syn_H	n_mis_H	wmis_cv_H	qmis_cv_H	AA		blastp	IGV Notes (All appear to be a single haplotype)
    # 1	maker-Mar.3.4.6.p1_scaffold15-snap-gene-127.5	2	30	35.963162	2.42E-08	0	1	0.05840293	9.22E-07	312	all PEI	no hits	Low USA mapping, almost no SNVs in healthy, LOTS in PEI
    # 2	maker-Mar.3.4.6.p1_scaffold9-snap-gene-700.19	2	12	34.868098	1.64E-04	13	14	0.42658891	4.53E-02	112	11/2 USA	receptor-type tyrosine-protein phosphatase mu-like (Mizuhopecten yessoensis)	Higher CN in USA vs PEI
    # 3	maker-Mar.3.4.6.p1_scaffold11-snap-gene-477.12	1	9	105.558072	1.64E-04	3	3	0.57100099	5.11E-01	80	9/0 USA	no hits	single exon mutation cluster
    # 4	maker-Mar.3.4.6.p1_scaffold4-snap-gene-785.0	4	15	11.320769	9.25E-03	4	0	0	1.06E-05	211	All USA	uncharacterized protein (Crassostrea virginica)	
    # 5	maker-Mar.3.4.6.p1_scaffold6-snap-gene-524.0	17	41	3.491338	2.42E-02	10	35	1.33074384	4.77E-01	447	All USA	ATP-dependent DNA helicase PIF1 (Mytilus galloprovincialis)	All homozygous, both PEI and USA appear to be low CN, especially PEI

healthy_dnds$globaldnds
dnds_best_somatic$globaldnds
usa_somatic$globaldnds
pei_somatic$globaldnds

gene_count <- nrow(healthy_dnds$genemuts) # 23273
filter(sel_cv, n_mis >0) %>% nrow() # 5606
filter(sel_cv, n_mis >0) %>% nrow()/gene_count # 0.24088
filter(sel_cv_usa, n_mis >0) %>% nrow()/gene_count # 0.1361234
filter(sel_cv_pei, n_mis >0) %>% nrow()/gene_count # 0.1375414
