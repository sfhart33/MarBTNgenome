# From: C:\Users\shart\Metzger Lab Dropbox\Sam_data\Mya_genome\mitochondria\mito_somatypus_analysis.r
    # module load R/3.6.0
    # R
    library(tidyverse)
    library(ape)
    library(Biostrings)
    library(gridExtra)
    library(devtools)
    library(dndscv)

    setwd("/ssd3/Mar_genome_analysis/bwa_mapping/mito/all_samples/somatypus")
    snvs <- read.delim("Somatypus_SNVs_final.counts", header = TRUE)
    head(snvs)

#not in dloop region
snvs2 <- filter(snvs, pos < 12060 | pos > 12971)

# # plot freq for all snvs for each sample
    # pdf("test.pdf")
    # for(i in c(	"Href_f", "Husa_f", "Hpei_f", 
    # 	"Cpei0_f", "Cpei1_f", "Cpei2_f", "Cpei3_f",
    # 	"Cusa0a_f", "Cusa0b_f", "Cusa1_f", "Cusa2_f", "Cusa3_f", "Cusa4_f", "Cusa5_f",
    # 	"Tpei1_f", "Tpei2_f", "Tpei3_f",
    # 	"Tusa1_f", "Tusa2_f", "Tusa3_f", "Tusa4_f", "Tusa5_f")){
    #     plot1 <- ggplot(snvs2, aes(get(i)))+
    #         geom_histogram(binwidth=0.01)+
    #         ggtitle(i)+
    #         theme_classic() +
    #         theme(axis.text=element_text(size=12,face="bold"),
    #                 axis.title=element_text(size=16,face="bold"),
    #                 text=element_text(size=18,face="bold"))
    #     print(plot1)
    # }
    # dev.off()

    # pdf("mito_snp_freq.pdf")
    # for(i in c(	"Href_f", "Husa_f", "Hpei_f", 
    # 	"Cpei0_f", "Cpei1_f", "Cpei2_f", "Cpei3_f",
    # 	"Cusa0a_f", "Cusa0b_f", "Cusa1_f", "Cusa2_f", "Cusa3_f", "Cusa4_f", "Cusa5_f")){ #	"Tpei1_f", "Tpei2_f", "Tpei3_f", "Tusa1_f", "Tusa2_f", "Tusa3_f", "Tusa4_f", "Tusa5_f"
    #     plot1 <- ggplot(snvs2, aes(get(i)))+
    #         geom_histogram(binwidth=0.01)+
    #         ggtitle(i)+
    #         theme_classic() +
    #         theme(axis.text=element_text(size=12,face="bold"),
    #                 axis.title=element_text(size=16,face="bold"),
    #                 text=element_text(size=18,face="bold"))
    #     print(plot1)
    #     plot1 <- ggplot(snvs2, aes(get(i)))+
    #         geom_histogram(binwidth=0.001)+
    #         xlim(0,0.1)+
    #         ggtitle(i)+
    #         theme_classic() +
    #         theme(axis.text=element_text(size=12,face="bold"),
    #                 axis.title=element_text(size=16,face="bold"),
    #                 text=element_text(size=18,face="bold"))
    #     print(plot1)
    #     plot1 <- ggplot(snvs2, aes(get(i)))+
    #         geom_histogram(binwidth=0.001)+
    #         xlim(0.9,1)+
    #         ggtitle(i)+
    #         theme_classic() +
    #         theme(axis.text=element_text(size=12,face="bold"),
    #                 axis.title=element_text(size=16,face="bold"),
    #                 text=element_text(size=18,face="bold"))
    #     print(plot1)
    # }
    # dev.off()

# Filter for snvs not found in healthies
    snvs_noH <- filter(snvs2, Href_f < 0.5, Husa_f < 0.5, Hpei_f < 0.5)
    snvs_inH <- filter(snvs2, Href_f > 0.5 | Husa_f > 0.5 | Hpei_f > 0.5)

# Median allele frequency to look for host contamination
    # for(i in c(	"Href", "Husa", "Hpei", 
    #     "Cpei0", "Cpei1", "Cpei2", "Cpei3",
    #     "Cusa0a", "Cusa0b", "Cusa1", "Cusa2", "Cusa3", "Cusa4", "Cusa5")){
    #     print(i)
    #     #filter(snvs, get(i) > 0.5) %>% pull(get(i)) %>% mean() %>% print()
    #     filter(snvs_noH, get(paste0(i,"_f")) > 0.5) %>% pull(get(paste0(i,"_f"))) %>% median() %>% print()
    # }
    # for(i in c(	"Href", "Husa", "Hpei", 
    #     "Cpei0", "Cpei1", "Cpei2", "Cpei3",
    #     "Cusa0a", "Cusa0b", "Cusa1", "Cusa2", "Cusa3", "Cusa4", "Cusa5")){
    #     print(i)
    #     #filter(snvs, get(i) > 0.5) %>% pull(get(i)) %>% mean() %>% print()
    #     filter(snvs_inH, get(paste0(i,"_f")) > 0.5) %>% pull(get(paste0(i,"_f"))) %>% median() %>% print()
    # }

    # for(i in c(	"Href", "Husa", "Hpei", 
    #     "Cpei0", "Cpei1", "Cpei2", "Cpei3",
    #     "Cusa0a", "Cusa0b", "Cusa1", "Cusa2", "Cusa3", "Cusa4", "Cusa5")){
    #     print(i)
    #     #filter(snvs, get(i) > 0.5) %>% pull(get(i)) %>% mean() %>% print()
    #     filter(snvs, get(paste0(i,"_f")) > 0.5) %>% pull(get(paste0(i,"_a"))) %>% median() %>% print()
    # }
    # for(i in c(	"Href", "Husa", "Hpei", 
    #     "Cpei0", "Cpei1", "Cpei2", "Cpei3",
    #     "Cusa0a", "Cusa0b", "Cusa1", "Cusa2", "Cusa3", "Cusa4", "Cusa5")){
    #     print(i)
    #     #filter(snvs, get(i) > 0.5) %>% pull(get(i)) %>% mean() %>% print()
    #     filter(snvs, get(paste0(i,"_f")) > 0.5) %>% pull(get(paste0(i,"_t"))) %>% median() %>% print()
    # }
    #######
    samples <- NULL
    withH <- NULL
    noH <- NULL
    for(i in c(	"Href", "Husa", "Hpei", 
        "Cpei0", "Cpei1", "Cpei2", "Cpei3",
        "Cusa0a", "Cusa0b", "Cusa1", "Cusa2", "Cusa3", "Cusa4", "Cusa5")){
        print(i)
        #filter(snvs, get(i) > 0.5) %>% pull(get(i)) %>% mean() %>% print()
        noH1 <- filter(snvs_noH, get(paste0(i,"_f")) > 0.5) %>% pull(get(paste0(i,"_f"))) %>% median()
        withH1 <- filter(snvs, get(paste0(i,"_f")) > 0.5) %>% pull(get(paste0(i,"_f"))) %>% median()
        noH <- c(noH, noH1)
        withH <- c(withH, withH1)
        samples <- c(samples , i)
    }
        purity <- data.frame(samples,noH,withH) 
        healthy_corr <- mean(purity[1:3,3])
        purity <- purity %>%
            mutate(withH_purity = 1-withH) %>%
            mutate(noH_purity = 1-noH) %>%
            mutate(withH_purity_corr = healthy_corr-withH) %>%
            print()
    pdf("sample_purity_plots.pdf")
    ggplot(purity, aes(samples,withH_purity)) +
            geom_col()+
            theme_classic() +
                theme(axis.text=element_text(size=12,face="bold"),
                axis.title=element_text(size=16,face="bold"),
                text=element_text(size=18,face="bold")) +
            ggtitle("Mitochondiral SNP contamination estimate")
    ggplot(purity, aes(samples,withH_purity_corr)) +
            geom_col()+
            theme_classic() +
                theme(axis.text=element_text(size=12,face="bold"),
                axis.title=element_text(size=16,face="bold"),
                text=element_text(size=18,face="bold")) +
            ggtitle("Mitochondiral SNP contamination estimate - healthy corrected")
    dev.off()

# comparison for pairwise differences
    # set cutoff parameters
        top=0.6
        bottom=0.1
    # samples
        samplesBTN=c("Cpei0", "Cpei1", "Cpei3", "Cusa0a", "Cusa0b", "Cusa2", "Cusa4", "Cusa5", "Href", "Husa", "Hpei", "genome") # Host cont to exclude: "Cusa1", "Cusa3", "Cpei2"
    # calculate differences
        pairwise=data.frame(A=character(), B=character(), count=numeric(), stringsAsFactors = FALSE)
        for(A in samplesBTN){
            if(A != "genome"){
                altA1=top
                altA2=bottom
            }
            for(B in samplesBTN){
                if(B != "genome"){
                    altB1=top
                    altB2=bottom
                }
                if(A == "genome"){
                    if(B != "genome"){ 
                        count <- filter(snvs2, get(paste0(B,"_f"))>altB1) %>% nrow()
                    }
                }
                if(B == "genome"){
                    if(A != "genome"){
                        count <- filter(snvs2, get(paste0(A,"_f"))>altA1) %>% nrow()
                    }
                }
                if(A != "genome" & B != "genome"){
                    count <- filter(snvs2, (get(paste0(A,"_f"))>altA1 & get(paste0(B,"_f"))<altB2) |
                                                (get(paste0(A,"_f"))<altA2 & get(paste0(B,"_f"))>altB1)
                                                ) %>% nrow()
                }
                if(A == "genome" & B == "genome"){
                    count=0
                }
            pairwise[nrow(pairwise) + 1,] = c(A,B,count)
            #print(paste(A,B,count,"done"))
            }
        }
    # Combine into single table
        #head(pairwise)
        col1 <- filter(pairwise, A =="Cpei0") %>% .$count
        col2 <- filter(pairwise, A =="Cpei1") %>% .$count
        col3 <- filter(pairwise, A =="Cpei3") %>% .$count
        col4 <- filter(pairwise, A =="Cusa0a") %>% .$count
        col5 <- filter(pairwise, A =="Cusa0b") %>% .$count
        col6 <- filter(pairwise, A =="Cusa2") %>% .$count
        col7 <- filter(pairwise, A =="Cusa4") %>% .$count
        col8 <- filter(pairwise, A =="Cusa5") %>% .$count
        col9 <- filter(pairwise, A =="Href") %>% .$count
        col10 <- filter(pairwise, A =="Husa") %>% .$count
        col11 <- filter(pairwise, A =="Hpei") %>% .$count
        col12 <- filter(pairwise, A =="genome") %>% .$count
        #extras
        col13 <- filter(pairwise, A =="Cpei2") %>% .$count
        col14 <- filter(pairwise, A =="Cusa1") %>% .$count
        col15 <- filter(pairwise, A =="Cusa3") %>% .$count
        pairwise_table <- data.frame(Cpei0 = as.numeric(col1), Cpei1 = as.numeric(col2), Cpei3 = as.numeric(col3),
                                    Cusa0a = as.numeric(col4), Cusa0b = as.numeric(col5), Cusa2 = as.numeric(col6), Cusa4 = as.numeric(col7), Cusa5 = as.numeric(col8), 
                                    Href = as.numeric(col9), Husa = as.numeric(col10), Hpei = as.numeric(col11), genome = as.numeric(col12))
    # label with real aliases
        samplesBTN_alias=c("DF-488", "DN-HL03", "PEI-DN08", "FFM-19G1", "FFM-20B2", "FFM-22F10", "MELC-A11", "NYTC-C9", "Href", "Husa", "Hpei", "genome") # Host cont to exclude: "Cusa1", "Cusa3", "Cpei2"
        rownames(pairwise_table) <- samplesBTN_alias
        colnames(pairwise_table) <- samplesBTN_alias
        pairwise_table
    # print tree to plot with FigTree
        njtree <- root(nj(as.dist(pairwise_table)), outgroup = "genome")
        write.tree(njtree, file = "BTN_mt.tree")
        njtree2 <- root(nj(as.dist(pairwise_table)), outgroup = "Href")
        write.tree(njtree2, file = "BTN_mt2.tree")


# MUTATIONAL BIASES AND CODING MUTATIONS


    #cutoff parameter
        cut=0.5

    snvs_H <- snvs2 %>% filter(Href_f > cut | Husa_f > cut | Hpei_f > cut)
        nrow(snvs_H) # 39
    snvs_noH <- snvs2 %>% filter(Href_f < cut & Husa_f < cut & Hpei_f < cut)
    snvs_any <- snvs2 %>% filter(Href_f > cut | Husa_f > cut | Hpei_f > cut |
                                Cpei0_f > cut | Cpei1_f > cut | Cpei3_f > cut | 
                                Cusa0a_f > cut | Cusa0b_f > cut | Cusa2_f > cut | Cusa4_f > cut | Cusa5_f > cut )
        nrow(snvs_any) # 102
    snvs_anyC <- snvs_noH %>% filter(Cpei0_f > cut | Cpei1_f > cut | Cpei3_f > cut | 
                                        Cusa0a_f > cut | Cusa0b_f > cut | Cusa2_f > cut | Cusa4_f > cut | Cusa5_f > cut )
        nrow(snvs_anyC) # 63
    snvs_allC <- snvs_noH %>% filter(Cpei0_f > cut & Cpei1_f > cut & Cpei3_f > cut & 
                                        Cusa0a_f > cut & Cusa0b_f > cut & Cusa2_f > cut & Cusa4_f > cut & Cusa5_f > cut )
        nrow(snvs_allC) # 13
    snvs_someC <- snvs_anyC %>% filter(!(Cpei0_f > cut & Cpei1_f > cut & Cpei3_f > cut & 
                                        Cusa0a_f > cut & Cusa0b_f > cut & Cusa2_f > cut & Cusa4_f > cut & Cusa5_f > cut ))
        nrow(snvs_someC) # 50    
    snvs_allpei <- snvs_noH %>% filter((Cpei0_f > cut & Cpei1_f > cut & Cpei3_f > cut) & 
                                        !(Cusa0a_f > cut | Cusa0b_f > cut | Cusa2_f > cut | Cusa4_f > cut | Cusa5_f > cut))
        nrow(snvs_allpei) # 26
    snvs_allusa <- snvs_noH %>% filter(!(Cpei0_f > cut | Cpei1_f > cut | Cpei3_f > cut) &
                                        (Cusa0a_f > cut & Cusa0b_f > cut & Cusa2_f > cut & Cusa4_f > cut & Cusa5_f > cut))
        nrow(snvs_allusa) # 21
    snvs_subsublineage <- snvs_noH %>% filter((!(Cpei0_f > cut & Cpei1_f > cut & Cpei3_f > cut) &
                                                (Cpei0_f > cut | Cpei1_f > cut | Cpei3_f > cut)) |
                                                (!(Cusa0a_f > cut & Cusa0b_f > cut & Cusa2_f > cut & Cusa4_f > cut & Cusa5_f > cut) &
                                                (Cusa0a_f > cut | Cusa0b_f > cut | Cusa2_f > cut | Cusa4_f > cut | Cusa5_f > cut)))
    nrow(snvs_subsublineage) # 3 

# count nuc freq
    mt <- readDNAStringSet("/ssd3/Mar_genome_analysis/genomes/mito/mt_genome.fasta")
    mt_nucfreq <- oligonucleotideFrequency(mt, as.prob=TRUE, width=1)         #   A C G T # 0.2856188 0.1172898 0.2256645 0.371427
    mt_masked <- readDNAStringSet("/ssd3/Mar_genome_analysis/genomes/mito/mt_genome_maskedloop.fasta")    
        oligonucleotideFrequency(mt_masked, as.prob=TRUE, width=1)  # 0.2791899 0.1177576 0.2270619 0.3759906

# count mutation freq
    count_nucs <- function(data, name){
        data2 <- data.frame(matrix(ncol=5,nrow=0, dimnames=list(NULL, c("sample","ref", "alt", "count", "expected"))))
        i = 0
        for(ref_nuc in c("A","T","G","C")){
            for(alt_nuc in c("A","T","G","C")){
                if(ref_nuc != alt_nuc){
                    i=i+1
                    nuc_count <- filter(data, ref == ref_nuc, alt == alt_nuc) %>% nrow()
                    pos_count <- nrow(data)*mt_nucfreq[1,ref_nuc]/3
                    data2[i,] <- c(name, ref_nuc, alt_nuc, nuc_count, pos_count)
                }
            }
        }
        return(data2)
    }

# plot mutation freq
    plot_mt_muts <- function(data, title){
        obs <- select(data, -expected)
        obs$data <- rep("obs", nrow(data))
        exp <- data %>%
            mutate(count = expected) %>%
            select(-expected)
        exp$data <- rep("exp", nrow(data))
        data2 <- rbind(obs,exp) %>%
            mutate(name = as.factor(paste0(ref,">",alt,"_",data))) %>%
            mutate(count = as.numeric(count)) %>%
            arrange(name)
        #print(data2)
        ggplot(data2,aes(x=name,y=count,fill=data)) +
            geom_col()+
            scale_fill_manual(values=c("grey","black"))+
            xlab("Mutation in + stand direction") +
            theme_classic() +
            theme(axis.text=element_text(size=12,face="bold"),
                axis.title=element_text(size=16,face="bold"),
                text=element_text(size=18,face="bold"),
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
            ggtitle(title)
    }

    s <- list()
    s[[1]] <- plot_mt_muts(count_nucs(snvs_H, "snvs_H"), "Healthy SNVs (39)")
    s[[2]] <- plot_mt_muts(count_nucs(snvs_allC, "snvs_allC"), "All Cancer (13)")
    s[[3]] <- plot_mt_muts(count_nucs(snvs_someC, "snvs_someC"), "Likely somatic(50)")
    s[[4]] <- plot_mt_muts(count_nucs(snvs_allusa, "snvs_allusa"), "All USA (21)")
    s[[5]] <- plot_mt_muts(count_nucs(snvs_allpei, "snvs_allpei"), "All PEI (26)")


    pdf("mut_bias_plot.pdf", height=20,width=10)
        #do.call(grid.arrange,s)
        grid.arrange(grobs = s, ncol=1)
    dev.off()

# Count up biases for written
    total_bias <- count_nucs(snvs_any, "total")
        as.numeric(total_bias$count) %>% sum() # 102
    somatic_bias <- count_nucs(snvs_someC, "somatic")
        as.numeric(somatic_bias$count) %>% sum() # 50
        filter(somatic_bias, (ref == "C" & alt == "T") | (ref == "G" & alt == "A")) %>%
            .$count %>%
            as.numeric %>%
            sum() # 41
        as.numeric(somatic_bias$count[])
    H_bias <- count_nucs(snvs_H, "snvs_H")
        as.numeric(H_bias$count) %>% sum() # 39



# Previous trinucleotide analysis, switch to just single nucleotide
        # helmsmanfile <- function(input, name){
        #     input %>% select(chr, pos, ref, alt) %>%
        #         mutate(id = name) %>%
        #         write.table(file = paste0("helmsman.txt"), append = TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep = '\t')
        # }
        # helmsmanfile(snvs_H, "anyH")
        # helmsmanfile(snvs_allC, "allC")
        # helmsmanfile(snvs_allpei, "allPEInoUSA")
        # helmsmanfile(snvs_allusa, "allUSAnoPEI")
        # helmsmanfile(snvs_somatic, "somatic")
#         mt_masked <- readDNAStringSet("/ssd3/Mar_genome_analysis/genomes/mito/mt_genome_maskedloop.fasta")
#         mt_masked2 <- trinucleotideFrequency(mt_masked, as.prob=TRUE)
#         mut <- as.data.frame(mt_masked2)
#         mut_C <- mutate(mut,
#             "ACA>" = (ACA + TGT)/3, 
#             "ACC>" = (ACC + GGT)/3, 
#             "ACG>" = (ACG + CGT)/3, 
#             "ACT>" = (ACT + AGT)/3,
#             "CCA>" = (CCA + TGG)/3, 
#             "CCC>" = (CCC + GGG)/3, 
#             "CCG>" = (CCG + CGG)/3, 
#             "CCT>" = (CCT + AGG)/3,
#             "GCA>" = (GCA + TGC)/3, 
#             "GCC>" = (GCC + GGC)/3, 
#             "GCG>" = (GCG + CGC)/3, 
#             "GCT>" = (GCT + AGC)/3,
#             "TCA>" = (TCA + TGA)/3, 
#             "TCC>" = (TCC + GGA)/3, 
#             "TCG>" = (TCG + CGA)/3, 
#             "TCT>" = (TCT + AGA)/3)
#         mut_T <- mutate(mut,
#             "ATA>" = (ATA + TAT)/3, 
#             "ATC>" = (ATC + GAT)/3, 
#             "ATG>" = (ATG + CAT)/3, 
#             "ATT>" = (ATT + AAT)/3,
#             "CTA>" = (CTA + TAG)/3, 
#             "CTC>" = (CTC + GAG)/3, 
#             "CTG>" = (CTG + CAG)/3, 
#             "CTT>" = (CTT + AAG)/3,
#             "GTA>" = (GTA + TAC)/3, 
#             "GTC>" = (GTC + GAC)/3, 
#             "GTG>" = (GTG + CAC)/3, 
#             "GTT>" = (GTT + AAC)/3,
#             "TTA>" = (TTA + TAA)/3, 
#             "TTC>" = (TTC + GAA)/3, 
#             "TTG>" = (TTG + CAA)/3, 
#             "TTT>" = (TTT + AAA)/3) 
#         mut_opp <- cbind(mut_C[,65:80],mut_C[,65:80],mut_C[,65:80],mut_T[,65:80],mut_T[,65:80],mut_T[,65:80])
#         mut_opp
#         library(sigfit)
#         data("cosmic_signatures_v2")
#         colnames(mut_opp) <- colnames(cosmic_signatures_v2)
#         rownames(mut_opp) <-  c("mt_genome_dloop_masked")
#         mut_opp
#         rowSums(mut_opp) # all equal 1
#         plot_spectrum(mut_opp, pdf_path = "mt_dloopmasked_mutational_oppertunities.pdf")
#         saveRDS(mut_opp, file = "mutational_oppertunities_mt_dloopmasked.rds")
#         write.table(mut_opp, file = "mutational_oppertunities_mt_dloopmasked.txt")
#     # run helmsman outside R
#         OUTPUT=/ssd3/Mar_genome_analysis/bwa_mapping/mito/all_samples/somatypus/helmsman
#         SNPS=/ssd3/Mar_genome_analysis/bwa_mapping/mito/all_samples/somatypus/helmsman.txt
#         GENOME=/ssd3/Mar_genome_analysis/genomes/mito/mt_genome.fasta
#         cd /ssd3/Mar_genome_analysis/mut_sig/helmsman
#         module load conda/
#         source activate helmsman
#         cd /ssd3/Mar_genome_analysis/bwa_mapping/mito/all_samples/somatypus
#         python /ssd3/Mar_genome_analysis/mut_sig/helmsman/helmsman.py --verbose --decomp nmf --mode txt --input $SNPS --fastafile $GENOME --projectdir $OUTPUT --length 3
#         python /ssd3/Mar_genome_analysis/mut_sig/helmsman/helmsman.py --verbose --decomp nmf --mode txt --input $SNPS --fastafile $GENOME --projectdir $OUTPUT/length1 --length 1
#     # back in R
#         helmsman_length1_file <- "/ssd3/Mar_genome_analysis/bwa_mapping/mito/all_samples/somatypus/helmsman/length1/subtype_count_matrix.txt"
#         helmsman_length1 <- read.table(helmsman_output_file,sep="\t",skip = 1)

#         helmsman_output_file <- "/ssd3/Mar_genome_analysis/bwa_mapping/mito/all_samples/somatypus/helmsman/subtype_count_matrix.txt"
#         helmsman <- read.table(helmsman_output_file,
#                                         sep="\t",
#                         skip = 1)
#         rownames(helmsman) <- helmsman$V1
#         helmsman <- select(helmsman, -V1)
#         colnames(helmsman) <- colnames(cosmic_signatures_v2)
#         helmsman <- helmsman[order(row.names(helmsman)), ]
#         head(helmsman)
#         plot_spectrum(helmsman, pdf_path = "mt_snps_spectrum.pdf")


# dnds and coding changes
# downloaded from published mya genome, manipulate outside R
    # head /ssd3/Mar_genome_analysis/genomes/mito/mt_genome.gff3

    # awk ' $3 ~ "CDS" {
    #     split($9, names, ";" );
    #     split(names[6], cds, "=" );
    #     split($1, genome, "." );
    #     length_cds = $5 - $4 + 1;
    #     if($7 ~ "+"){strand = 1};
    #     if($7 ~ "-"){strand = -1};
    #     print cds[2], cds[2], cds[2], genome[1], $4, $5, 1, length_cds, length_cds, strand} '  OFS='\t'  /ssd3/Mar_genome_analysis/genomes/mito/mt_genome.gff3 > /ssd3/Mar_genome_analysis/genomes/mito/mt_genome.dnds.input

# back in R
    library(devtools)
    library(dndscv)

# Build reference for dndscv
    buildref(cdsfile="/ssd3/Mar_genome_analysis/genomes/mito/mt_genome.dnds.input",
            genomefile="/ssd3/Mar_genome_analysis/genomes/mito/mt_genome.fasta",
            outfile = "/ssd3/Mar_genome_analysis/genomes/mito//mt_genome_refCDS.rda",
            numcode=5)

# (Old notes)
    # all_snvs <- read.delim("helmsman.txt", head=FALSE) %>%
    #     select(V5,V1,V2,V3,V4) %>%
    #     write.table(file = paste0("helmsman.txt"), append = TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep = '\t')

# Create dndscv file for each bin
    dndscvfile <- function(input, name){
        input %>% 
        mutate(id = name) %>% 
        select(id, chr, pos, ref, alt) %>%
            write.table(file = paste0(name,".dnds.txt"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep = '\t')
    }
    dndscvfile(snvs_H, "anyH")
    dndscvfile(snvs_allC, "allC")
    dndscvfile(snvs_allpei, "allPEInoUSA")
    dndscvfile(snvs_allusa, "allUSAnoPEI")
    dndscvfile(snvs_someC, "somatic")

# run dnds on various bins
    SH_dndscv <- function(input){
    mutations_list <- read.table(input,
                                    sep="\t")
    dndsout <- dndscv(mutations_list,
                            refdb="/ssd3/Mar_genome_analysis/genomes/mito//mt_genome_refCDS.rda",
                            cv=NULL,
                            outp = 1,
                            max_coding_muts_per_sample = 1000000000,
                            max_muts_per_gene_per_sample = 1000000000)
    return(dndsout)
    }
    anyH_dnds <- SH_dndscv("anyH.dnds.txt")
    allC_dnds <- SH_dndscv("allC.dnds.txt")
    PEI_dnds <- SH_dndscv("allPEInoUSA.dnds.txt")
    USA_dnds <- SH_dndscv("allUSAnoPEI.dnds.txt")
    somatic_dnds <- SH_dndscv("somatic.dnds.txt")

    comparison_dnds <- rbind(anyH_dnds$globaldnds[1,], allC_dnds$globaldnds[1,], somatic_dnds$globaldnds[1,])
    comparison_dnds$name <- as.factor(c("Healthy", "AllC_noH", "Somatic"))
    comparison_dnds$name <- ordered(comparison_dnds$name, levels = c("Healthy", "AllC_noH", "Somatic"))
# Plot
    pdf("mt_dnds.pdf")
    ggplot(comparison_dnds) +
        #geom_col(aes(name, mle))+
        geom_pointrange(aes(x=name, y=mle, ymin=cilow, ymax=cihigh))+
            #scale_fill_manual(values=c("grey","black"))+
            xlab("Sample bin") +
            ylab("dN/dS (neutral = 1)") +
            theme_classic() +
            theme(axis.text=element_text(size=12,face="bold"),
                axis.title=element_text(size=16,face="bold"),
                text=element_text(size=18,face="bold") #,
                #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
                ) +
            ggtitle("Mitochondrial SNVs dN/dS")
    dev.off()



# (Old notes)
    # visualize in igv instead
    # all_snvs <- read.delim("helmsman.txt", head=FALSE)
    # all_snvs2 <- all_snvs %>% mutate(chr = paste0(V1,".1"), start=V2-1, end =V2) %>%
    #     select(chr,start,end, V5) %>%
    #     arrange(start)

    # all_snvs3 <- all_snvs %>% filter(V5 != "somatic") %>%
    #     mutate(chr = paste0(V1,".1"), start=V2-1, end =V2) %>%
    #     select(chr,start,end, V5, V3, V4) %>%
    #     arrange(start)

    # snvs_somatic_bed
    # write.table(all_snvs3, file = "mt_snvs2.bed", quote=FALSE, row.names=FALSE, col.names=FALSE)

    # write.table(somatic_dnds$annotmuts, file = "somatic_annotmuts", quote=FALSE, row.names=FALSE, col.names=TRUE,sep="\t")
    # write.table(anyH_dnds$annotmuts, file = "anyH_annotmuts", quote=FALSE, row.names=FALSE, col.names=TRUE,sep="\t")

