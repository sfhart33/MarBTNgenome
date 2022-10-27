# Original file here: C:\Users\shart\Metzger Lab Dropbox\Sam_data\Mya_genome\steamer\steamer_downstream.r

library(tidyverse) # Need to summarize the data and for %>% function
library(ape) # for pairwise phylogeny

# LOAD INPUTS
    setwd("/ssd3/Mar_genome_analysis/steamer/final_pipeline/merge")
    sample1 <- "MELC-2E11" # Healthy REFERENCE
    sample2 <- "MELC-A9" # Healthy
    sample3 <- "PEI-DF490" # Healthy
    sample4 <- "FFM-19G1" # Cancerous
    sample5 <- "FFM-20B2" # Cancerous
    sample6 <- "FFM-22F10" # Cancerous
    sample7 <- "MELC-A11_S1" # Cancerous
    sample8 <- "NYTC-C9_S2" # Cancerous
    sample9 <- "DF-488" # Cancerous
    sample10 <- "DN-HL03" # Cancerous
    sample11 <- "PEI-DN08_S3" # Cancerous
    for(mult in c("updown","multiple")){
        for(i in 1:11){
            # Steamer insertion sites
                sample <- paste0("/ssd3/Mar_genome_analysis/steamer/final_pipeline/", get(paste0("sample",i)), "_",mult,".bed")
                sample_file <- read.delim(sample) %>%
                    filter( total > 5 | (up > 0 & down>0)) %>%
                    select(contig,start,end, strand, name, total)
                names(sample_file) <- c("contig", "min", "max", "strand", "name", str_replace(get(paste0("sample",i)), "-", "_"))
            # Depth of reads info
                depth <- paste0("/ssd3/Mar_genome_analysis/steamer/final_pipeline/coverage/",get(paste0("sample",i)),"_",mult,"_depth.bed")
                depth_file <- read.delim(depth, header = FALSE, col.names = c("name", paste0(str_replace(get(paste0("sample",i)), "-", "_"),"_cov")))
                col_name <- paste0(str_replace(i, "-", "_"),"cov")
            # Merge the two
                merged <- right_join(sample_file, depth_file, by = "name") #%>%
                    #mutate(!!col_name := get(str_replace(get(paste0("sample",i)), "-", "_")) / get(str_replace(get(paste0("sample",i)), "-", "_cov_")))
                assign(paste0("sample",i,"_premerge"),merged)
            }
        # MERGE INTO ONE FILE
            merge_all <- full_join(sample1_premerge, sample2_premerge, by = c("contig", "min", "max", "strand", "name")) 
            for(i in 3:11){
                merge_all <- full_join(merge_all, get(paste0("sample",i,"_premerge")), by = c("contig", "min", "max", "strand", "name"))
            }
            merge_all[is.na(merge_all)] <- 0
            nrow(merge_all)
            if(mult == "updown"){ merge_all_updown <- merge_all}
    }

    head(merge_all)
    head(merge_all_updown)
    nrow(merge_all)
    nrow(merge_all_updown)


## Phylogeny using correct bootstrapping

samplesALL
sites_binary <- merge_all %>%
    mutate(MELC_2E11 = case_when(MELC_2E11 > 0 ~ 1, MELC_2E11 == 0 ~ 0),
           MELC_A9 = case_when(MELC_A9 > 0 ~ 1, MELC_A9 == 0 ~ 0),
           PEI_DF490 = case_when(PEI_DF490 > 0 ~ 1, PEI_DF490 == 0 ~ 0),
           FFM_19G1 = case_when(FFM_19G1 > 0 ~ 1, FFM_19G1 == 0 ~ 0),
           FFM_20B2 = case_when(FFM_20B2 > 0 ~ 1, FFM_20B2 == 0 ~ 0),
           FFM_22F10 = case_when(FFM_22F10 > 0 ~ 1, FFM_22F10 == 0 ~ 0),
           MELC_A11 = case_when(MELC_A11_S1 > 0 ~ 1, MELC_A11_S1 == 0 ~ 0),
           NYTC_C9 = case_when(NYTC_C9_S2 > 0 ~ 1, NYTC_C9_S2 == 0 ~ 0),
           DF_488 = case_when(DF_488 > 0 ~ 1, DF_488 == 0 ~ 0),
           DN_HL03 = case_when(DN_HL03 > 0 ~ 1, DN_HL03 == 0 ~ 0),
           PEI_DN08 = case_when(PEI_DN08_S3 > 0 ~ 1, PEI_DN08_S3 == 0 ~ 0)
    ) %>%
    select(MELC_2E11,MELC_A9,PEI_DF490,FFM_19G1,FFM_20B2,FFM_22F10,MELC_A11,NYTC_C9,DF_488,DN_HL03,PEI_DN08)

sites_binary_input <- sites_binary[rowSums(sites_binary)>0,] %>% t() %>% as.matrix()
colnames(sites_binary_input)<-NULL



  rootedNJtree <- function(alignment, theoutgroup)
  {
     # this function requires the ape and seqinR packages:
     require("ape")
     require("seqinr")
     # define a function for making a tree:
     makemytree <- function(alignmentmat, outgroup=`theoutgroup`)
     {
        mydist <- dist.gene(alignmentmat)
        mytree <- nj(mydist)
        mytree <- makeLabel(mytree, space="") # get rid of spaces in tip names.
        myrootedtree <- root(mytree, outgroup, r=TRUE)
        return(myrootedtree)
     }
     # infer a tree
     mymat  <- alignment
     mytree <- makemytree(mymat)
     # bootstrap the tree
     myboot <- boot.phylo(mytree, mymat, makemytree)
     # plot the tree:
     plot.phylo(mytree, type="p")   # plot the unrooted phylogenetic tree
     nodelabels(myboot,cex=0.7)    # plot the bootstrap values
     mytree$node.label <- myboot   # make the bootstrap values be the node labels
     return(mytree)
  }                                       

    setwd("/ssd3/Mar_genome_analysis/revision/phylogeny/steamer")
    pdf("nj.steamer.tree.pdf")
    steamer_tree = nj(dist.gene(sites_binary_input)) 
    plot.phylo(steamer_tree, type="u")
    dev.off()

    pdf("nj.steamer.tree.pdf")
    datatree <- rootedNJtree(sites_binary_input, "MELC_A9")
    dev.off()
    saveRDS(datatree, file="nj.steamer.tree.rds")
    write.tree(datatree, file = "nj.steamer.tree")
                                



# Count healthys
    filter(merge_all, MELC_2E11 > 0) %>% nrow()
    filter(merge_all, MELC_A9 > 0) %>% nrow()
    filter(merge_all, PEI_DF490 > 0 ) %>% nrow()
    filter(merge_all, FFM_19G1 > 0) %>% nrow()
    filter(merge_all, FFM_20B2 > 0) %>% nrow()
    filter(merge_all, FFM_22F10 > 0 ) %>% nrow()
    filter(merge_all, MELC_A11_S1 > 0) %>% nrow()
    filter(merge_all, NYTC_C9_S2 > 0) %>% nrow()
    filter(merge_all, DF_488 > 0 ) %>% nrow()
    filter(merge_all, DN_HL03 > 0) %>% nrow()
    filter(merge_all, PEI_DN08_S3 > 0) %>% nrow()
    filter(merge_all_updown, MELC_2E11 > 0) %>% nrow()
    filter(merge_all_updown, MELC_A9 > 0) %>% nrow()
    filter(merge_all_updown, PEI_DF490 > 0 ) %>% nrow()
    filter(merge_all_updown, FFM_19G1 > 0) %>% nrow()
    filter(merge_all_updown, FFM_20B2 > 0) %>% nrow()
    filter(merge_all_updown, FFM_22F10 > 0 ) %>% nrow()
    filter(merge_all_updown, MELC_A11_S1 > 0) %>% nrow()
    filter(merge_all_updown, NYTC_C9_S2 > 0) %>% nrow()
    filter(merge_all_updown, DF_488 > 0 ) %>% nrow()
    filter(merge_all_updown, DN_HL03 > 0) %>% nrow()
    filter(merge_all_updown, PEI_DN08_S3 > 0) %>% nrow()

# other bins
    no_healthy <- filter(merge_all, MELC_2E11 == 0 & MELC_A9 == 0 & PEI_DF490 == 0 )
    all_cancer <- filter(no_healthy,
                            FFM_19G1 > 0 & FFM_20B2 > 0 & MELC_A11_S1 > 0 & FFM_22F10 > 0 & NYTC_C9_S2 > 0 &
                            DF_488 > 0 & DN_HL03 > 0 & PEI_DN08_S3 > 0)%>%
                    mutate(score = (FFM_19G1+FFM_20B2+MELC_A11_S1+FFM_22F10+NYTC_C9_S2+DF_488+DN_HL03+PEI_DN08_S3)/
                                (FFM_19G1_cov+FFM_20B2_cov+MELC_A11_S1_cov+FFM_22F10_cov+
                                NYTC_C9_S2_cov+DF_488_cov+DN_HL03_cov+PEI_DN08_S3))
    all_cancer2 <- filter(no_healthy,
                            (FFM_19G1 > 0 | FFM_20B2 > 0 | MELC_A11_S1 > 0 | FFM_22F10 > 0 | NYTC_C9_S2 > 0) &
                            (DF_488 > 0 | DN_HL03 > 0 | PEI_DN08_S3 > 0))
    any_cancer <- filter(no_healthy,
                            (FFM_19G1 > 0 | FFM_20B2 > 0 | MELC_A11_S1 > 0 | FFM_22F10 > 0 | NYTC_C9_S2 > 0 |
                            DF_488 > 0 | DN_HL03 > 0 | PEI_DN08_S3 > 0))
    no_pei <- filter(no_healthy,
                            DF_488 == 0 & DN_HL03 == 0 & PEI_DN08_S3 == 0)
    no_usa <- filter(no_healthy,
                            FFM_19G1 == 0 & FFM_20B2 == 0 & MELC_A11_S1 == 0 & FFM_22F10 == 0 & NYTC_C9_S2 == 0)
    any_pei <- filter(no_usa, DF_488 > 0 | DN_HL03 > 0 | PEI_DN08_S3 > 0)
    all_pei <- filter(no_usa,
                            DF_488 > 0 & DN_HL03 > 0 & PEI_DN08_S3 > 0)%>%
                    mutate(score = (DF_488+DN_HL03+PEI_DN08_S3)/(DF_488_cov+DN_HL03_cov+PEI_DN08_S3))
    notall_pei <- filter(no_usa,
                            (DF_488 == 0 | DN_HL03 == 0 | PEI_DN08_S3 == 0) & (DF_488 > 0 | DN_HL03 > 0 | PEI_DN08_S3 > 0))                      
    all_usa <- filter(no_pei,
                            FFM_19G1 > 0 & FFM_20B2 > 0 & MELC_A11_S1 > 0 & FFM_22F10 > 0 & NYTC_C9_S2 > 0)%>%
                    mutate(score = (FFM_19G1+FFM_20B2+MELC_A11_S1+FFM_22F10+NYTC_C9_S2)/
                                (FFM_19G1_cov+FFM_20B2_cov+MELC_A11_S1_cov+FFM_22F10_cov+NYTC_C9_S2_cov))
    notall_usa <- filter(no_pei,
                            (FFM_19G1 == 0 | FFM_20B2 == 0 | MELC_A11_S1 == 0 | FFM_22F10 == 0 | NYTC_C9_S2 == 0) &
                            (FFM_19G1 > 0 | FFM_20B2 > 0 | MELC_A11_S1 > 0 | FFM_22F10 > 0 | NYTC_C9_S2 > 0))
    somatic <- rbind(all_cancer,all_usa,all_pei)

    nrow(no_healthy) # 1061
    nrow(any_cancer) # 550
    nrow(all_cancer) # 193
    nrow(all_cancer2) # 238
    nrow(all_usa) # 187
    nrow(all_pei) # 51
    nrow(notall_usa) # 61
    nrow(notall_pei) # 13
    nrow(somatic) # 431

# Compare multiple read support to JUST updown support
    no_healthy_updown <- filter(merge_all_updown, MELC_2E11 == 0 & MELC_A9 == 0 & PEI_DF490 == 0 )
    all_cancer_updown <- filter(merge_all_updown, MELC_2E11 == 0 & MELC_A9 == 0 & PEI_DF490 == 0 &
                            FFM_19G1 > 0 & FFM_20B2 > 0 & MELC_A11_S1 > 0 & FFM_22F10 > 0 & NYTC_C9_S2 > 0 &
                            DF_488 > 0 & DN_HL03 > 0 & PEI_DN08_S3 > 0)
    any_cancer_updown <- filter(no_healthy_updown,
                            (FFM_19G1 > 0 | FFM_20B2 > 0 | MELC_A11_S1 > 0 | FFM_22F10 > 0 | NYTC_C9_S2 > 0 |
                            DF_488 > 0 | DN_HL03 > 0 | PEI_DN08_S3 > 0))
    nrow(no_healthy_updown) # 401
    nrow(all_cancer_updown) # 131
    nrow(any_cancer_updown) # 401

# Merge to find discrepancies between multiple and updown (possibly mediating rearrangement?)
    head(merge_all_updown)
    merge_all_updown2 <- merge_all_updown
    colnames(merge_all_updown2)[6:27] <- paste0(colnames(merge_all_updown2)[6:27], "_updown")
    head(merge_all_updown2)
    merge_mult_updown <- full_join(merge_all, merge_all_updown2, by = c("contig", "min", "max", "strand", "name")) 
    merge_mult_updown[is.na(merge_mult_updown)] <- 0

    all_cancer_mult_updown <- filter(merge_mult_updown, MELC_2E11 == 0 & MELC_A9 == 0 & PEI_DF490 == 0 &
                            FFM_19G1 > 0 & FFM_20B2 > 0 & MELC_A11_S1 > 0 & FFM_22F10 > 0 & NYTC_C9_S2 > 0 &
                            DF_488 > 0 & DN_HL03 > 0 & PEI_DN08_S3 > 0)
    all_cancer_mult_updown_disc <- filter(all_cancer_mult_updown,
                            FFM_19G1_updown == 0 | FFM_20B2_updown == 0 | MELC_A11_S1_updown == 0 | FFM_22F10_updown == 0 | NYTC_C9_S2_updown == 0 |
                            DF_488_updown == 0 | DN_HL03_updown == 0 | PEI_DN08_S3 == 0)
    all_cancer_mult_updown_disc_USA <- filter(all_cancer_mult_updown_disc,
                            FFM_19G1_updown > 0 & FFM_20B2_updown > 0 & MELC_A11_S1_updown > 0 & FFM_22F10_updown > 0 & NYTC_C9_S2_updown > 0 &
                            DF_488_updown == 0 & DN_HL03_updown == 0 & PEI_DN08_S3 == 0)
    all_cancer_mult_updown_disc_PEI <- filter(all_cancer_mult_updown_disc,
                            FFM_19G1_updown == 0 & FFM_20B2_updown == 0 & MELC_A11_S1_updown == 0 & FFM_22F10_updown == 0 & NYTC_C9_S2_updown == 0 &
                            DF_488_updown > 0 & DN_HL03_updown > 0 & PEI_DN08_S3 > 0)
    nrow(all_cancer_mult_updown_disc) # 61
    nrow(all_cancer_mult_updown_disc_USA) # 0
    nrow(all_cancer_mult_updown_disc_PEI) # 2
    # Overall, doesn't look like much evidence for Steamer-mediated rearrangement

# Print bed files of insertion sites
    print_steamer_bed <- function(file1, name){
        if(!"score" %in% colnames(file1)){
                file1 <- mutate(file1, score = 0)
        }
        file1 %>%
            filter(substr(contig, 1, 12) == "Mar.3.4.6.p1") %>%
            arrange(as.character(contig),min) %>%
            select(contig,min,max,name,score,strand) %>%
            write.table(file = name, row.names=FALSE, col.names=FALSE, quote=FALSE, sep = '\t')
    }

    print_steamer_bed(somatic,"allCnoH_allPEI_allUSA.bed")
    print_steamer_bed(all_cancer,"allCnoH.bed")
    print_steamer_bed(any_cancer,"anyCnoH.bed")
    print_steamer_bed(all_usa,"allUSA.bed")
    print_steamer_bed(all_pei,"allPEI.bed")
    print_steamer_bed(merge_all,"all_sites.bed")

# allele frequency plots

all_cancer_freq <- mutate(all_cancer,
                                    usa_freq = (FFM_19G1/FFM_19G1_cov + FFM_20B2/FFM_20B2_cov + FFM_22F10/FFM_22F10_cov + MELC_A11_S1/MELC_A11_S1_cov + NYTC_C9_S2/NYTC_C9_S2_cov)/5,
                                    pei_freq = (DF_488/DF_488_cov + DN_HL03/DN_HL03_cov + PEI_DN08_S3/PEI_DN08_S3_cov)/3)
all_usa_freq <- mutate(all_usa,
                                    usa_freq = (FFM_19G1/FFM_19G1_cov + FFM_20B2/FFM_20B2_cov + FFM_22F10/FFM_22F10_cov + MELC_A11_S1/MELC_A11_S1_cov + NYTC_C9_S2/NYTC_C9_S2_cov)/5,
                                    pei_freq = (DF_488/DF_488_cov + DN_HL03/DN_HL03_cov + PEI_DN08_S3/PEI_DN08_S3_cov)/3)
all_pei_freq <- mutate(all_pei,
                                    usa_freq = (FFM_19G1/FFM_19G1_cov + FFM_20B2/FFM_20B2_cov + FFM_22F10/FFM_22F10_cov + MELC_A11_S1/MELC_A11_S1_cov + NYTC_C9_S2/NYTC_C9_S2_cov)/5,
                                    pei_freq = (DF_488/DF_488_cov + DN_HL03/DN_HL03_cov + PEI_DN08_S3/PEI_DN08_S3_cov)/3)
                                    
                                    # usa_freq = (FFM_19G1 + FFM_20B2 + FFM_22F10 + MELC_A11_S1 + NYTC_C9_S2)/((FFM_19G1_cov + FFM_20B2_cov + FFM_22F10_cov + MELC_A11_S1_cov + NYTC_C9_S2_cov) - (FFM_19G1 + FFM_20B2 + FFM_22F10 + MELC_A11_S1 + NYTC_C9_S2)/2),
                                    # pei_freq = (DF_488 + DN_HL03 + PEI_DN08_S3)/((DF_488_cov + DN_HL03_cov + PEI_DN08_S3_cov)-(DF_488 + DN_HL03 + PEI_DN08_S3)/2))

# pdf("all_cancer_freq5.pdf")
#     ggplot(all_cancer_freq) +
#         geom_freqpoly(aes(usa_freq), binwidth=0.05, color = "blue") +
#         geom_freqpoly(aes(pei_freq), binwidth=0.05, color = "red") +
#             theme_bw() +
#             theme(axis.text=element_text(size=12,face="bold"),
#                     axis.title=element_text(size=14,face="bold"),
#                     text=element_text(size=14,face="bold"),
#                     plot.title = element_text(hjust = 0.5)) +
#             xlim(0,1.5)+
#             xlab("Allele freq (avg for sub-lineage)") +
#             ggtitle("Allele freq of Steamer insertions (allCnoH)")
#     ggplot() +
#         geom_freqpoly(data = all_usa_freq, aes(usa_freq), binwidth=0.05, color = "blue") +
#         geom_freqpoly(data = all_pei_freq, aes(pei_freq), binwidth=0.05, color = "red") +
#             theme_bw() +
#             theme(axis.text=element_text(size=12,face="bold"),
#                     axis.title=element_text(size=14,face="bold"),
#                     text=element_text(size=14,face="bold"),
#                     plot.title = element_text(hjust = 0.5)) +
#             xlim(0,1.5)+
#             xlab("Allele freq (avg for sub-lineage)") +
#             ggtitle("Allele freq of Steamer insertions (sublineage only)")
#     ggplot(all_cancer_freq) +
#         geom_point(aes(usa_freq, pei_freq)) +
#             theme_bw() +
#             theme(axis.text=element_text(size=12,face="bold"),
#                     axis.title=element_text(size=14,face="bold"),
#                     text=element_text(size=14,face="bold"),
#                     plot.title = element_text(hjust = 0.5)) +
#             xlim(0,1.5)+
#             ylim(0,1.5)+
#             ggtitle("Allele freq of Steamer insertions (allCnoH)")
# dev.off()

            summary_freq <- rbind(
                mutate(all_cancer_freq, bin = "allC_PEI", freq = pei_freq, sub = "PEI"),
                mutate(all_cancer_freq, bin = "allC_USA", freq = usa_freq, sub = "USA"),
                mutate(all_usa_freq, bin = "allUSA_noPEI", freq = usa_freq,sub = "USA"),
                mutate(all_pei_freq, bin = "allPEI_noUSA", freq = pei_freq, sub = "PEI")
                )

    pdf("steamer_freq_violin.pdf", width=8, height=4)
        ggplot(summary_freq, aes(x=as.factor(bin), y=freq, color=as.factor(sub))) +
            geom_hline(aes(yintercept=1), color="black", size=0.25)+
            geom_hline(aes(yintercept=3/4), color="black", size=0.25)+
            geom_hline(aes(yintercept=1/4), color="black", size=0.25)+
            geom_hline(aes(yintercept=1/3), color="black", size=0.25)+
            geom_hline(aes(yintercept=2/3), color="black", size=0.25)+
            geom_hline(aes(yintercept=1/2), color="black", size=0.25)+
            geom_hline(aes(yintercept=0), color="black", size=0.25)+
            geom_violin() + #, binwidth = 0.01
            scale_fill_manual(values=c("red", "blue")) +
            # ylim(0,1.25)+
            scale_y_continuous(limits = c(0,1.25),breaks=c(0,1/4,1/3,1/2,2/3,3/4,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1)) +
            ylab("Steamer insertion allele freq estimate")+
            xlab("Bin")+
            theme_classic() +
            theme(axis.text=element_text(size=8,face="bold"),
                axis.title=element_text(size=12,face="bold"),
                text=element_text(size=12,face="bold"),
                legend.title = element_blank()) #+
            #ggtitle("USA_somaticSNV_CN_comparison")
        #print(plot1)
    dev.off()

#### NEW BELOW
# Steamer pileups
# LOAD INPUTS - just those supported by both up and down
    setwd("/ssd3/Mar_genome_analysis/steamer/final_pipeline/merge")
    sample1 <- "MELC-2E11" # Healthy REFERENCE
    sample2 <- "MELC-A9" # Healthy
    sample3 <- "PEI-DF490" # Healthy
    sample4 <- "FFM-19G1" # Cancerous
    sample5 <- "FFM-20B2" # Cancerous
    sample6 <- "FFM-22F10" # Cancerous
    sample7 <- "MELC-A11_S1" # Cancerous
    sample8 <- "NYTC-C9_S2" # Cancerous
    sample9 <- "DF-488" # Cancerous
    sample10 <- "DN-HL03" # Cancerous
    sample11 <- "PEI-DN08_S3" # Cancerous
    for(mult in c("updown","multiple")){
        for(i in 1:11){
            # Steamer insertion sites
                sample <- paste0("/ssd3/Mar_genome_analysis/steamer/final_pipeline/", get(paste0("sample",i)), "_",mult,".bed")
                sample_file <- read.delim(sample) %>%
                    filter((up > 0 & down>0)) %>%
                    select(contig,start,end, strand, name, total)
                names(sample_file) <- c("contig", "min", "max", "strand", "name", str_replace(get(paste0("sample",i)), "-", "_"))
            # Depth of reads info
                depth <- paste0("/ssd3/Mar_genome_analysis/steamer/final_pipeline/coverage/",get(paste0("sample",i)),"_",mult,"_depth.bed")
                depth_file <- read.delim(depth, header = FALSE, col.names = c("name", paste0(str_replace(get(paste0("sample",i)), "-", "_"),"_cov")))
                col_name <- paste0(str_replace(i, "-", "_"),"cov")
            # Merge the two
                merged <- right_join(sample_file, depth_file, by = "name") #%>%
                    #mutate(!!col_name := get(str_replace(get(paste0("sample",i)), "-", "_")) / get(str_replace(get(paste0("sample",i)), "-", "_cov_")))
                assign(paste0("sample",i,"_premerge"),merged)
            }
        # MERGE INTO ONE FILE
            merge_all <- full_join(sample1_premerge, sample2_premerge, by = c("contig", "min", "max", "strand", "name")) 
            for(i in 3:11){
                merge_all <- full_join(merge_all, get(paste0("sample",i,"_premerge")), by = c("contig", "min", "max", "strand", "name"))
            }
            merge_all[is.na(merge_all)] <- 0
            nrow(merge_all)
            if(mult == "updown"){ merge_all_updown <- merge_all}
    }

# other bins
    no_healthy <- filter(merge_all, MELC_2E11 == 0 & MELC_A9 == 0 & PEI_DF490 == 0 )
    all_cancer <- filter(no_healthy,
                            FFM_19G1 > 0 & FFM_20B2 > 0 & MELC_A11_S1 > 0 & FFM_22F10 > 0 & NYTC_C9_S2 > 0 &
                            DF_488 > 0 & DN_HL03 > 0 & PEI_DN08_S3 > 0)%>%
                    mutate(score = (FFM_19G1+FFM_20B2+MELC_A11_S1+FFM_22F10+NYTC_C9_S2+DF_488+DN_HL03+PEI_DN08_S3)/
                                (FFM_19G1_cov+FFM_20B2_cov+MELC_A11_S1_cov+FFM_22F10_cov+
                                NYTC_C9_S2_cov+DF_488_cov+DN_HL03_cov+PEI_DN08_S3))
    all_cancer2 <- filter(no_healthy,
                            (FFM_19G1 > 0 | FFM_20B2 > 0 | MELC_A11_S1 > 0 | FFM_22F10 > 0 | NYTC_C9_S2 > 0) &
                            (DF_488 > 0 | DN_HL03 > 0 | PEI_DN08_S3 > 0))
    any_cancer <- filter(no_healthy,
                            (FFM_19G1 > 0 | FFM_20B2 > 0 | MELC_A11_S1 > 0 | FFM_22F10 > 0 | NYTC_C9_S2 > 0 |
                            DF_488 > 0 | DN_HL03 > 0 | PEI_DN08_S3 > 0))
    no_pei <- filter(no_healthy,
                            DF_488 == 0 & DN_HL03 == 0 & PEI_DN08_S3 == 0)
    no_usa <- filter(no_healthy,
                            FFM_19G1 == 0 & FFM_20B2 == 0 & MELC_A11_S1 == 0 & FFM_22F10 == 0 & NYTC_C9_S2 == 0)

    all_pei <- filter(no_usa,
                            DF_488 > 0 & DN_HL03 > 0 & PEI_DN08_S3 > 0)%>%
                    mutate(score = (DF_488+DN_HL03+PEI_DN08_S3)/(DF_488_cov+DN_HL03_cov+PEI_DN08_S3))                  
    all_usa <- filter(no_pei,
                            FFM_19G1 > 0 & FFM_20B2 > 0 & MELC_A11_S1 > 0 & FFM_22F10 > 0 & NYTC_C9_S2 > 0)%>%
                    mutate(score = (FFM_19G1+FFM_20B2+MELC_A11_S1+FFM_22F10+NYTC_C9_S2)/
                                (FFM_19G1_cov+FFM_20B2_cov+MELC_A11_S1_cov+FFM_22F10_cov+NYTC_C9_S2_cov))


    steamer_ins <- rbind(mutate(all_cancer, bin = "allBTN"),mutate(all_usa, bin = "allUSA"),mutate(all_pei, bin = "allPEI")) %>% 
        separate(col = contig, into = c(NA, "contig"), sep = "scaffold") %>% # contig to scaffold number so that I can sort, note it is still a character until changed
        arrange(contig,min) %>% # sort
        mutate(lead = lead(min) - min,
               lag = min - lag(min)) # find how close to surrounding insertions

    select(steamer_ins, contig, min, lead, lag)
    filter(steamer_ins, contig == 11)

    filter(steamer_ins, lead > 0 & lag > 0 & (lead < 10000 | lag < 10000)) %>%
        select(name, contig, min, lead, lag, bin)
        # Manually found in /ssd3/Mar_genome_analysis/steamer/final_pipeline/genes/allCnoH_allPEI_allUSA.proximity ->
        # /ssd3/Mar_genome_analysis/steamer/final_pipeline/merge/pileups_manual.txt

    filter(steamer_ins, lead > 0 & lag > 0 & (lead < 100000 & lag < 100000)) %>%
        select(name, contig, min, lead, lag, bin)
    
