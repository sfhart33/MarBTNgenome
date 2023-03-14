library(tidyverse)

setwd("/ssd3/Mar_genome_analysis/steamer/final_pipeline/genes")

prox_files <- list.files(pattern = "*.proximity") %>% print()
# [1] "all_sites.proximity"             "allCnoH_allPEI_allUSA.proximity"
# [3] "allCnoH.proximity"               "allPEI.proximity"
# [5] "allUSA.proximity"                "anyCnoH.proximity"

for(file in prox_files) {
    assign(file, read.table(file ,sep="\t", col.names=c("insertion","gene","distance")))
}


pdf("proximity_histograms.pdf")
for(file in prox_files){
    # file <- "all_sites.proximity"
    # for(widths in c(50,100,200,500)){    
    plot1 <- ggplot(get(file))+
                geom_histogram(aes(distance), binwidth = 250, color = "black") +  #, bins=100
                xlim(-15000,15000)+
                theme_bw() +
                theme(axis.text=element_text(size=12,face="bold"),
                        axis.title=element_text(size=14,face="bold"),
                        text=element_text(size=14,face="bold")) +
                xlab("distance from nearest gene (bins=250bp)") +
                ggtitle(file)
    print(plot1)
}
dev.off()

MELC_2E11.subset.proximity <- read.table("MELC-2E11.subset.proximity", sep="\t", col.names=c("insertion","gene","distance"))
pdf("proximity_histograms_testfile2.pdf")   
    plot1 <- ggplot(MELC_2E11.subset.proximity)+
                geom_histogram(aes(distance), binwidth = 250, color = "black") +  #, bins=100
                xlim(-15000,15000)+
                theme_bw() +
                theme(axis.text=element_text(size=12,face="bold"),
                        axis.title=element_text(size=14,face="bold"),
                        text=element_text(size=14,face="bold")) +
                xlab("distance from nearest gene (bins=250bp)") +
                ggtitle("MELC-2E11.subset.proximity")
    print(plot1)
dev.off()


filter(all_sites.proximity, distance > -2000, distance < 0) %>% nrow() # 118
filter(all_sites.proximity, distance > 0, distance < 2000) %>% nrow() # 30

filter(anyCnoH.proximity, distance > -2000, distance < 0) %>% nrow() # 115
filter(anyCnoH.proximity, distance > 0, distance < 2000) %>% nrow() # 29

filter(allCnoH.proximity, distance > -2000, distance < 0) %>% nrow() # 40
filter(allCnoH.proximity, distance > 0, distance < 2000) %>% nrow() # 11
filter(allUSA.proximity, distance > -2000, distance < 0) %>% nrow() # 45
filter(allUSA.proximity, distance > 0, distance < 2000) %>% nrow() # 7
filter(allPEI.proximity, distance > -2000, distance < 0) %>% nrow() # 11
filter(allPEI.proximity, distance > 0, distance < 2000) %>% nrow() # 2

regions <- read.table("regions_intersect.txt" ,sep="\t",row.names = 1,header=T) %>%
    t() %>%
    as.data.frame() %>%
    mutate(sublin = allUSA+allPEI) %>%
    as.matrix()
rownames(regions) <- read.table("regions_intersect.txt" ,sep="\t",row.names = 1,header=T) %>% t() %>% rownames()

regions2 <- regions[,c(1:6,8)]/regions[,7] # normalize to genome size
regions3 <- regions2/regions2["genome",]  # normalize to genome insertion rate

regions2["genes",]/regions2["genome",]   
regions2["cds",]/regions2["genome",]        
regions2["upstream1000",]/regions2["genome",]    
regions2["upstream2000",]/regions2["genome",]    
regions2["fiveprime",]/regions2["genome",]  
genome_correction <- regions2["genome",]  

regions_tests <- read.table("regions_intersect_tests.txt" ,sep="\t",row.names = 1, header=T) %>% t()
regions_tests2 <- regions_tests[,1]/regions_tests[,2]
mapping_correction <- regions_tests2/regions_tests2["genome"]

# Make bar charts comparing regions by steamer insertions (normalize by both bp and read mapping)
    regions2
    genome_correction
    mapping_correction
    corrected_output <- as.data.frame((t(t(as.matrix(regions2))/genome_correction)))/mapping_correction
    corrected_output_plot <- corrected_output[c(2,6,4,7),c(4,7)]
    # corrected_output_plot %>% select(allCnoH) %>% mutate(sample = "allCnoH") %>% add_rownames(var = "region")
    steamer_upstream_plot <- data.frame(
                                        samples = as.factor(c("allCnoH","allCnoH","allCnoH","allCnoH","sublineages","sublineages","sublineages","sublineages")),
                                        regions = factor(c("genes","CDS","upstream","5'UTR","genes","CDS","upstream","5'UTR"), levels=c("genes","CDS","upstream","5'UTR")),
                                        values = c(corrected_output_plot[,"allCnoH"],corrected_output_plot[,"sublin"]))
    pdf("steamer_upstream_barchart.pdf")
        ggplot(steamer_upstream_plot) +
            geom_bar(aes(fill=regions, x=samples, y=values), position="dodge", stat="identity") +
                theme_bw() +
                theme(axis.text=element_text(size=20),
                        axis.title=element_text(size=20),
                        text=element_text(size=20),
                        plot.title = element_text(hjust = 0.5)) +
                xlab("Sample bin") +
                ylab("Normalized Steamer insertions \nCorrected for read mapping rates")# +
                #ggtitle("Steamer insertions \nNormalized to full genome \nCorrected for read mapping")
    dev.off()

# # NOT USED - IN RELATION TO START CODON RATHER THAN GENE
# # get start codons from dndscv, correct for strand and 0-based bed format
#     load(file = "/ssd3/Mar_genome_analysis/dnds/dndscv/Mar.3.4.6.p1_snap02_refCDS.rda")
#     gene_count <- nrow(RefCDS)
#     starts_bed <- data.frame(matrix(ncol=6,nrow=0, dimnames=list(NULL, c("chr", "start", "end","gene","length","strand"))))
#     for(i in 1:gene_count){
#         if(RefCDS[[i]]$strand==1){
#             coding=RefCDS[[i]]$intervals_cds[1,1]
#             start=coding-1
#             end=coding+2
#             starts_bed[i,]<- c(RefCDS[[i]]$chr, start, end, RefCDS[[i]]$gene_name, RefCDS[[i]]$CDS_length,"+")
#         }
#         if(RefCDS[[i]]$strand==-1){
#             coding=RefCDS[[i]]$intervals_cds[nrow(RefCDS[[i]]$intervals_cds),2]
#             start=coding-3
#             end=coding
#             starts_bed[i,] <- c(RefCDS[[i]]$chr, start, end, RefCDS[[i]]$gene_name, RefCDS[[i]]$CDS_length,"-")
#         }
#     }
#     head(starts_bed)
#     write.table(starts_bed, "start_codons.bed", sep="\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
#         atg_bed <- data.frame(matrix(ncol=6,nrow=0, dimnames=list(NULL, c("chr", "start", "end","gene","length","strand"))))
#         count=0
#     for(i in 1:gene_count){
#         if(RefCDS[[i]]$strand==1 & as.character(RefCDS[[i]]$seq_cds[1:3]) == "ATG"){
#             count=count+1
#             coding=RefCDS[[i]]$intervals_cds[1,1]
#             start=coding-1
#             end=coding+2
#             atg_bed[count,]<- c(RefCDS[[i]]$chr, start, end, RefCDS[[i]]$gene_name, RefCDS[[i]]$CDS_length,"+")
#         }
#         if(RefCDS[[i]]$strand==-1 & as.character(RefCDS[[i]]$seq_cds[1:3]) == "ATG"){
#             count=count+1
#             coding=RefCDS[[i]]$intervals_cds[nrow(RefCDS[[i]]$intervals_cds),2]
#             start=coding-3
#             end=coding
#             atg_bed[count,] <- c(RefCDS[[i]]$chr, start, end, RefCDS[[i]]$gene_name, RefCDS[[i]]$CDS_length,"-")
#         }
#     }
#     head(atg_bed)
#     write.table(atg_bed, "start_codons_atg.bed", sep="\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
#     nrow(atg_bed)/nrow(starts_bed) # 0.8699781
# # after calculating distance to steamer insertions:
#     prox_files2 <- list.files(pattern = ".start_codons" ) %>% print()
#     for(file in prox_files2) {
#         assign(file, read.table(file ,sep="\t", col.names=c("insertion","gene","distance")))
#     }
#     pdf("proximity_histograms_start_codons_15kb.pdf")
#     for(file in prox_files2){
#         # file <- "all_sites.proximity"
#         # for(widths in c(50,100,200,500)){    
#         plot1 <- ggplot(get(file))+
#                     geom_histogram(aes(as.numeric(distance)), binwidth = 250, color = "black") +  #, bins=100
#                     xlim(-15000,15000)+
#                     theme_bw() +
#                     theme(axis.text=element_text(size=12,face="bold"),
#                             axis.title=element_text(size=14,face="bold"),
#                             text=element_text(size=14,face="bold")) +
#                     xlab("distance from nearest srat codon (bins=250bp)") +
#                     ggtitle(file)
#         print(plot1)
#     }
#     dev.off()
#     pdf("proximity_histograms_start_codons_50kb.pdf")
#     for(file in prox_files2){
#         # file <- "all_sites.proximity"
#         # for(widths in c(50,100,200,500)){    
#         plot1 <- ggplot(get(file))+
#                     geom_histogram(aes(as.numeric(distance)), binwidth = 1000, color = "black") +  #, bins=100
#                     xlim(-50000,50000)+
#                     theme_bw() +
#                     theme(axis.text=element_text(size=12,face="bold"),
#                             axis.title=element_text(size=14,face="bold"),
#                             text=element_text(size=14,face="bold")) +
#                     xlab("distance from nearest srat codon (bins=250bp)") +
#                     ggtitle(file)
#         print(plot1)
#     }
#     dev.off()
#     pdf("proximity_histograms_start_codons_1kb.pdf")
#     for(file in prox_files2){
#         # file <- "all_sites.proximity"
#         # for(widths in c(50,100,200,500)){    
#         plot1 <- ggplot(get(file))+
#                     geom_histogram(aes(as.numeric(distance)), binwidth = 10, color = "black") +  #, bins=100
#                     xlim(-1000,1000)+
#                     theme_bw() +
#                     theme(axis.text=element_text(size=12,face="bold"),
#                             axis.title=element_text(size=14,face="bold"),
#                             text=element_text(size=14,face="bold")) +
#                     xlab("distance from nearest srat codon (bins=250bp)") +
#                     ggtitle(file)
#         print(plot1)
#     }
#     dev.off()
#     pdf("proximity_histograms_start_codons_100bp.pdf")
#     for(file in prox_files2){
#         # file <- "all_sites.proximity"
#         # for(widths in c(50,100,200,500)){    
#         plot1 <- ggplot(get(file))+
#                     geom_histogram(aes(as.numeric(distance)), binwidth = 1, color = "black") +  #, bins=100
#                     xlim(-100,100)+
#                     theme_bw() +
#                     theme(axis.text=element_text(size=12,face="bold"),
#                             axis.title=element_text(size=14,face="bold"),
#                             text=element_text(size=14,face="bold")) +
#                     xlab("distance from nearest srat codon (bins=250bp)") +
#                     ggtitle(file)
#         print(plot1)
#     }
#     dev.off()
    