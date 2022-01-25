# From original file: C:\Users\shart\Metzger Lab Dropbox\Sam_data\Mya_genome\SNPs\somatypus_output_pairwise_noLOH.r

library(tidyverse)
library(ape)
library(sigfit)
library(lsa)
library(gridExtra)
data(cosmic_signatures_v2)
data(cosmic_signatures_v3)

setwd("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run")


# Load data if running on all SNVs
    # snvs.nr <- read.table("Somatypus_SNVs_final_NR.txt", header=T, check.names=F)
    # snvs.nv <- read.table("Somatypus_SNVs_final_NV.txt", header=T, check.names=F)
    # snvs.metadata <- read.table("Somatypus_SNVs_final_Metadata.txt", header=T)
    # setwd("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/pairwise")  
# Load data if running on SNVs outside possible LOH regions
    snvs.nr <- read.table("Somatypus_SNVs_final_noLOH2_NR.txt", header=T, check.names=F)
    snvs.nv <- read.table("Somatypus_SNVs_final_noLOH2_NV.txt", header=T, check.names=F)
    snvs.metadata <- read.table("Somatypus_SNVs_final_noLOH2_Metadata.txt", header=T)
    setwd("/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/pairwise/noLOH")

# Define sample names and combine into one table
    samples <- c("H_ref","H_usa","H_pei",
    "C_pei1","C_pei2","C_pei4","C_pei3", # C_pei4 has high host contamination
    "C_usa1","C_usa2","C_usa6","C_usa3","C_usa7","C_usa4","C_usa5", # C_usa6 & 7 have high host contamination
    "C_pei2T","C_pei4T","C_pei3T","C_usa6T","C_usa3T","C_usa7T","C_usa4T","C_usa5T")
    samplesNR <- paste0(samples, "_NR")
    samplesNV <- paste0(samples, "_NV")
    colnames(snvs.nr) <- samplesNR
    colnames(snvs.nv) <- samplesNV
    colnames(indels.nr) <- samplesNR
    colnames(indels.nv) <- samplesNV
    snvs.metadata2 <- snvs.metadata %>%
        select(CHROM,POS,REF,ALT)
    snvs <- cbind(snvs.metadata2,snvs.nv,snvs.nr)

# Calculate average read depth across the genome
    for(sample in samples){
        depth <- mean(snvs[,paste0(sample,"_NR")])
        assign(paste0(sample,"_avg"), depth)
    }
    for(sample in samples){
        print(sample)
        print(get(paste0(sample,"_avg")))
    }
    # average depths
        # [1] "H_ref"
        # [1] 59.91782
        # [1] "H_usa"
        # [1] 60.0784
        # [1] "H_pei"
        # [1] 60.3944
        # [1] "C_pei1"
        # [1] 94.64185
        # [1] "C_pei2"
        # [1] 89.04225
        # [1] "C_pei4"
        # [1] 56.09117
        # [1] "C_pei3"
        # [1] 65.3303
        # [1] "C_usa1"
        # [1] 92.81066
        # [1] "C_usa2"
        # [1] 92.78821
        # [1] "C_usa6"
        # [1] 68.59992
        # [1] "C_usa3"
        # [1] 63.84493
        # [1] "C_usa7"
        # [1] 58.19498
        # [1] "C_usa4"
        # [1] 63.49766
        # [1] "C_usa5"
        # [1] 59.25352
        # [1] "C_pei2T"
        # [1] 5.811725
        # [1] "C_pei4T"
        # [1] 18.54592
        # [1] "C_pei3T"
        # [1] 9.330324
        # [1] "C_usa6T"
        # [1] 31.91118
        # [1] "C_usa3T"
        # [1] 35.81058
        # [1] "C_usa7T"
        # [1] 33.20985
        # [1] "C_usa4T"
        # [1] 29.47542
        # [1] "C_usa5T"
        # [1] 31.76414


#Set thresholds
    min_cutoff=3
    stringent_depth=6 # was 8
    weak_depth=16

# Not found in any healthy, found in at least one cancer above stringent threshold
    realSNVs <- filter(snvs, !(H_ref_NV > min_cutoff |
                        H_usa_NV > min_cutoff |
                        H_pei_NV > min_cutoff) &
                        (C_pei1_NV > C_pei1_avg/stringent_depth |
                        C_pei2_NV > C_pei2_avg/stringent_depth |
                        C_pei3_NV > C_pei3_avg/stringent_depth |
                        C_usa1_NV > C_usa1_avg/stringent_depth |
                        C_usa2_NV > C_usa2_avg/stringent_depth |
                        C_usa3_NV > C_usa3_avg/stringent_depth |
                        C_usa4_NV > C_usa4_avg/stringent_depth |
                        C_usa5_NV > C_usa5_avg/stringent_depth))

# Function to write helmsman output
    # helmsmanfile <- function(input, sample1, sample2){
    #     name <- paste0(sample1,".",sample2)
    #     input %>% select(CHROM, POS, REF, ALT) %>%
    #         mutate(ID = name) %>%
    #         write.table(file = paste0("helmsman.txt"), append = TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep = '\t')
    #         #write.table(file = paste0("helmsman.txt"), append = TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep = '\t')
    # }
 

# Run pairwise analysis
    # samplesALL=c("H_ref","H_usa","H_pei","C_pei1","C_pei2","C_pei3","C_usa1","C_usa2","C_usa3","C_usa4","C_usa5","genome")
    samplesBTN=c("C_pei1","C_pei2","C_pei3","C_usa1","C_usa2","C_usa3","C_usa4","C_usa5","genome")
    pairwise=data.frame(A=character(), B=character(), count=numeric(), stringsAsFactors = FALSE)
    # if(file.exists("helmsman.txt")) {
    #     file.remove("helmsman.txt")
    # }
    for(A in samplesBTN){
        print(paste(A,"started"))
        if(A != "genome"){
            # altA1= get(paste0(A,"_avg"))/stringent_depth # high threshold
            altA1= get(paste0(A,"_avg"))/weak_depth # high threshold
            altA2= get(paste0(A,"_avg"))/weak_depth # low threshold
        }
        for(B in samplesBTN){
            if(B != "genome"){
                # altB1= get(paste0(B,"_avg"))/stringent_depth # high threshold
                altB1= get(paste0(B,"_avg"))/weak_depth # high threshold
                altB2= get(paste0(B,"_avg"))/weak_depth # low threshold
            }
            if(A == "genome"){
                if(B != "genome"){ 
                    count <- filter(realSNVs, get(paste0(B,"_NV"))>altB1) %>% nrow()
                            # filter(realSNVs, get(paste0(B,"_NV"))>altB1) %>% helmsmanfile(A,B)
                }
            if(B == "genome"){ 
                    count <- 0
                }
            }
            if(B == "genome"){
                if(A != "genome"){
                    count <- filter(realSNVs, get(paste0(A,"_NV"))>altA1) %>% nrow()
                             #filter(realSNVs, get(paste0(A,"_NV"))>altA1) %>% helmsmanfile(A,B)
                }
            }
            if(A != "genome" & B != "genome"){
                count <- filter(realSNVs, (get(paste0(A,"_NV"))>altA1 & get(paste0(B,"_NV"))<altB2) |
                                            (get(paste0(A,"_NV"))<altA2 & get(paste0(B,"_NV"))>altB1)
                                            ) %>% nrow()
                        #  filter(realSNVs, (get(paste0(A,"_NV"))>altA1 & get(paste0(B,"_NV"))<altB2) |
                        #                     (get(paste0(A,"_NV"))<altA2 & get(paste0(B,"_NV"))>altB1)
                        #                     ) %>% helmsmanfile(A,B)
            }
        pairwise[nrow(pairwise) + 1,] = c(A,B,count)
        print(paste("   ",B,count,"done"))
        }
    }
    #head(pairwise)
    col1 <- filter(pairwise, A =="C_pei1") %>% .$count
    col2 <- filter(pairwise, A =="C_pei2") %>% .$count
    col3 <- filter(pairwise, A =="C_pei3") %>% .$count
    # col4 <- filter(pairwise, A =="C_pei4") %>% .$count
    col5 <- filter(pairwise, A =="C_usa1") %>% .$count
    col6 <- filter(pairwise, A =="C_usa2") %>% .$count
    col7 <- filter(pairwise, A =="C_usa3") %>% .$count
    col8 <- filter(pairwise, A =="C_usa4") %>% .$count
    col9 <- filter(pairwise, A =="C_usa5") %>% .$count
    # col10 <- filter(pairwise, A =="C_usa6") %>% .$count
    # col11 <- filter(pairwise, A =="C_usa7") %>% .$count
    col12 <- filter(pairwise, A =="genome") %>% .$count
    pairwise_table <- data.frame(Cpei1 = as.numeric(col1), Cpei2 = as.numeric(col2), Cpei3 = as.numeric(col3), # Cpei4 = as.numeric(col4),
                                Cusa1 = as.numeric(col5), Cusa2 = as.numeric(col6), Cusa3 = as.numeric(col7),
                                Cusa4 = as.numeric(col8), Cusa5 = as.numeric(col9), genome = as.numeric(col12)) # , Cusa6 = as.numeric(col10), Cusa7 = as.numeric(col11)
    rownames(pairwise_table) <- samplesBTN
    pairwise_table_nogenome <- pairwise_table[1:8,1:8]
    pairwise_table

# Print output
    #pdf("BTN.tree.pdf")
    pdf("BTN.noLOH.tree.pdf")
        njtree <- root(nj(as.dist(pairwise_table)), outgroup = "genome", resolve.root = TRUE)
        plot.phylo(root(njtree, outgroup = "genome"))
    dev.off()
    #write.tree(njtree, file = "BTN.tree")
    #njtree_nogenome <- nj(as.dist(pairwise_table_nogenome))
    saveRDS(pairwise_table, file="BTN.pairwise_table.phy.rds")
    saveRDS(njtree, file="BTN.tree.phy.rds")
    #saveRDS(njtree_nogenome , file="BTNnogenome.tree.phy.rds")
    write.tree(njtree, file = "BTN.noLOH2.tree")

# NEW bootstrapping phylogeny
    #tree_build <- function(x){nj(as.dist(x))}
    tree_build <- function(x) root(nj(as.dist(x)), outgroup = "genome", resolve.root = TRUE)
    for(i in 1:10){
        boot.phylo(tree_build(pairwise_table), pairwise_table, tree_build, B=1000) %>% print()
        #Example output
        # Calculating bootstrap values... done.
        # [1] 1000   31   49   94   51   69  103 1000
    }


# plot in a grid 
    # p <- list()
    # i=0
    # for(A in samplesBTN){
    #     for(B in samplesBTN){
    #         if(A != B & A != "genome" & B != "genome" ){
    #             i=i+1
    #         # assign thresholds
    #             altA1= get(paste0(A,"_avg"))/stringent_depth # high threshold
    #             altB2= get(paste0(B,"_avg"))/weak_depth # low threshold
    #             plot1 <- filter(realSNVs, (get(paste0(A,"_NV"))>altA1 & get(paste0(B,"_NV"))<altB2))
    #         # plot
    #             plotalt(plot1, A, B) %>% print()
    #             p[[i]] <- plotalt(plot1, A, B)
    #             print(paste(A,"not",B,"done"))
    #         }
    #     }
    # }

    # pdf("pairwise_histograms.pdf", height=20, width=20)
    #     #grid.arrange(p[1:4], ncol=2) #length(samplesBTN))
    #     do.call(grid.arrange,p)
    # dev.off()


#################################################################################
# From original file: C:\Users\shart\Metzger Lab Dropbox\Sam_data\Mya_genome\SNPs\somatypus_output_initial_counts.r

#BINNING - same as for signature extraction
    # Is SNV found in any cancer
        anyC <- filter(snvs, C_pei1_NV > C_pei1_avg/stringent_depth |
                            C_pei2_NV > C_pei2_avg/stringent_depth |
                            C_pei3_NV > C_pei3_avg/stringent_depth |
                            C_usa1_NV > C_usa1_avg/stringent_depth |
                            C_usa2_NV > C_usa2_avg/stringent_depth |
                            C_usa3_NV > C_usa3_avg/stringent_depth |
                            C_usa4_NV > C_usa4_avg/stringent_depth |
                            C_usa5_NV > C_usa5_avg/stringent_depth)
        noC <- filter(snvs, !(C_pei1_NV > C_pei1_avg/stringent_depth |
                            C_pei2_NV > C_pei2_avg/stringent_depth |
                            C_pei3_NV > C_pei3_avg/stringent_depth |
                            C_usa1_NV > C_usa1_avg/stringent_depth |
                            C_usa2_NV > C_usa2_avg/stringent_depth |
                            C_usa3_NV > C_usa3_avg/stringent_depth |
                            C_usa4_NV > C_usa4_avg/stringent_depth |
                            C_usa5_NV > C_usa5_avg/stringent_depth))

    # Is SNV found in any healthy clam    
        anyHanyC <- filter(anyC, H_ref_NV > min_cutoff |
                            H_usa_NV > min_cutoff |
                            H_pei_NV > min_cutoff)
        noHanyC <- filter(anyC, !(H_ref_NV > min_cutoff |
                            H_usa_NV > min_cutoff |
                            H_pei_NV > min_cutoff))
        anyHnoC <- filter(noC, H_ref_NV > min_cutoff |
                            H_usa_NV > min_cutoff |
                            H_pei_NV > min_cutoff)
        noHnoC <- filter(noC, !(H_ref_NV > min_cutoff |
                            H_usa_NV > min_cutoff |
                            H_pei_NV > min_cutoff))
    # which healthy clam is it found in
        anyH <- rbind(anyHnoC, anyHanyC)  %>%
            arrange(CHROM,POS)
        H_R <- filter(noC, (H_ref_NV > min_cutoff) &
                            !(H_usa_NV > min_cutoff) &
                            !(H_pei_NV > min_cutoff))
        H_RU <- filter(noC, (H_ref_NV > min_cutoff) &
                            (H_usa_NV > min_cutoff) &
                            !(H_pei_NV > min_cutoff))
        H_RP <- filter(noC, (H_ref_NV > min_cutoff) &
                            !(H_usa_NV > min_cutoff) &
                            (H_pei_NV > min_cutoff))
        H_RUP <- filter(noC, (H_ref_NV > min_cutoff) &
                            (H_usa_NV > min_cutoff) &
                            (H_pei_NV > min_cutoff))
        H_U <- filter(noC, !(H_ref_NV > min_cutoff) &
                            (H_usa_NV > min_cutoff) &
                            !(H_pei_NV > min_cutoff))
        H_UP <- filter(noC, !(H_ref_NV > min_cutoff) &
                            (H_usa_NV > min_cutoff) &
                            (H_pei_NV > min_cutoff))      
        H_P <- filter(noC, !(H_ref_NV > min_cutoff) &
                            !(H_usa_NV > min_cutoff) &
                            (H_pei_NV > min_cutoff))                                            
    # Likely founder germline
        allCanyH <- filter(anyHanyC, C_pei1_NV > C_pei1_avg/weak_depth &
                            C_pei2_NV > C_pei2_avg/weak_depth &
                            C_pei3_NV > C_pei3_avg/weak_depth &
                            C_usa1_NV > C_usa1_avg/weak_depth &
                            C_usa2_NV > C_usa2_avg/weak_depth &
                            C_usa3_NV > C_usa3_avg/weak_depth &
                            C_usa4_NV > C_usa4_avg/weak_depth &
                            C_usa5_NV > C_usa5_avg/weak_depth)
    # Likely founder germline in LOH region or mismap
        subsetCanyH <- filter(anyHanyC, !(C_pei1_NV > C_pei1_avg/weak_depth &
                            C_pei2_NV > C_pei2_avg/weak_depth &
                            C_pei3_NV > C_pei3_avg/weak_depth &
                            C_usa1_NV > C_usa1_avg/weak_depth &
                            C_usa2_NV > C_usa2_avg/weak_depth &
                            C_usa3_NV > C_usa3_avg/weak_depth &
                            C_usa4_NV > C_usa4_avg/weak_depth &
                            C_usa5_NV > C_usa5_avg/weak_depth))
    # Founder germline OR early somatic mutation
        allCnoH <- filter(noHanyC, C_pei1_NV > C_pei1_avg/weak_depth &
                            C_pei2_NV > C_pei2_avg/weak_depth &
                            C_pei3_NV > C_pei3_avg/weak_depth &
                            C_usa1_NV > C_usa1_avg/weak_depth &
                            C_usa2_NV > C_usa2_avg/weak_depth &
                            C_usa3_NV > C_usa3_avg/weak_depth &
                            C_usa4_NV > C_usa4_avg/weak_depth &
                            C_usa5_NV > C_usa5_avg/weak_depth)
    # Somtic mutations OR mismaps, maybe some founder germline that dropped out of some samples
        subsetCnoH <- filter(noHanyC, !(C_pei1_NV > C_pei1_avg/weak_depth &
                            C_pei2_NV > C_pei2_avg/weak_depth &
                            C_pei3_NV > C_pei3_avg/weak_depth &
                            C_usa1_NV > C_usa1_avg/weak_depth &
                            C_usa2_NV > C_usa2_avg/weak_depth &
                            C_usa3_NV > C_usa3_avg/weak_depth &
                            C_usa4_NV > C_usa4_avg/weak_depth &
                            C_usa5_NV > C_usa5_avg/weak_depth))
    # Sublineage specific
        anyUSAnoPEInoH <- filter(subsetCnoH, !(C_pei1_NV > C_pei1_avg/weak_depth |
                            C_pei2_NV > C_pei2_avg/weak_depth |
                            C_pei3_NV > C_pei3_avg/weak_depth) &
                            (C_usa1_NV > C_usa1_avg/weak_depth |
                            C_usa2_NV > C_usa2_avg/weak_depth |
                            C_usa3_NV > C_usa3_avg/weak_depth |
                            C_usa4_NV > C_usa4_avg/weak_depth |
                            C_usa5_NV > C_usa5_avg/weak_depth))
        anyPEInoUSAnoH <- filter(subsetCnoH, (C_pei1_NV > C_pei1_avg/weak_depth |
                            C_pei2_NV > C_pei2_avg/weak_depth |
                            C_pei3_NV > C_pei3_avg/weak_depth) &
                            !(C_usa1_NV > C_usa1_avg/weak_depth |
                            C_usa2_NV > C_usa2_avg/weak_depth |
                            C_usa3_NV > C_usa3_avg/weak_depth |
                            C_usa4_NV > C_usa4_avg/weak_depth |
                            C_usa5_NV > C_usa5_avg/weak_depth))        
        allUSAnoPEInoH <- filter(subsetCnoH, !(C_pei1_NV > C_pei1_avg/weak_depth |
                            C_pei2_NV > C_pei2_avg/weak_depth |
                            C_pei3_NV > C_pei3_avg/weak_depth) &
                            (C_usa1_NV > C_usa1_avg/weak_depth &
                            C_usa2_NV > C_usa2_avg/weak_depth &
                            C_usa3_NV > C_usa3_avg/weak_depth &
                            C_usa4_NV > C_usa4_avg/weak_depth &
                            C_usa5_NV > C_usa5_avg/weak_depth))
        allPEInoUSAnoH <- filter(subsetCnoH, (C_pei1_NV > C_pei1_avg/weak_depth &
                            C_pei2_NV > C_pei2_avg/weak_depth &
                            C_pei3_NV > C_pei3_avg/weak_depth) &
                            !(C_usa1_NV > C_usa1_avg/weak_depth |
                            C_usa2_NV > C_usa2_avg/weak_depth |
                            C_usa3_NV > C_usa3_avg/weak_depth |
                            C_usa4_NV > C_usa4_avg/weak_depth |
                            C_usa5_NV > C_usa5_avg/weak_depth))

nrow(anyC)
nrow(noC)
nrow(anyHanyC)
nrow(noHanyC)
nrow(anyHnoC)
nrow(noHnoC)
nrow(anyH)
nrow(allCanyH)
nrow(allCnoH)
nrow(anyUSAnoPEInoH)
nrow(anyPEInoUSAnoH)
nrow(allUSAnoPEInoH)
nrow(allPEInoUSAnoH)
    # > nrow(anyC)
    # [1] 10,435,090
    # > nrow(noC)
    # [1] 16,855,311
    # > nrow(anyHanyC)
    # [1] 7,297,123
    # > nrow(noHanyC)
    # [1] 3,137,967
    # > nrow(anyHnoC)
    # [1] 9,002,766
    # > nrow(noHnoC)
    # [1] 7852545
    # > nrow(anyH)
    # [1] 16,299,889
    # > nrow(allCanyH)
    # [1] 6,858,790
    # > nrow(allCnoH)
    # [1] 2,264,329
    # > nrow(anyUSAnoPEInoH)
    # [1] 369,653
    # > nrow(anyPEInoUSAnoH)
    # [1] 330,589
    # > nrow(allUSAnoPEInoH)
    # [1] 331,167
    # > nrow(allPEInoUSAnoH)
    # [1] 320,715
