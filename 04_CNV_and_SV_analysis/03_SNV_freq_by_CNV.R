
# Original file here: C:\Users\shart\Metzger Lab Dropbox\Sam_data\Mya_genome\CNV_final\CNV_calling_pipeline.r

# Load packages
    library(cn.mops)
    library(tidyverse)
    #library(vcfR)
    library(mixtools)
    #library(karyoploteR)
    #library(fastseg)
    library(zoo)
    library(bedr)

    # load regions
    load_cn_snvs <- function(region,snvs,loh,method,cn){
        table_load <- read.table(paste0("/ssd3/Mar_genome_analysis/CNV_calling/FINAL/SNVs/",region,"_",snvs,loh,method,"_CN",cn,".bed"), header=F, check.names=F)
        colnames(table_load) <- c("chr","start","end","pei1_NV","pei2_NV","pei3_NV","usa1_NV","usa2_NV","usa3_NV","usa4_NV","usa5_NV","pei1_NR","pei2_NR","pei3_NR","usa1_NR","usa2_NR","usa3_NR","usa4_NR","usa5_NR")
        table_output <- mutate(table_load,
                    pei1_f = pei1_NV/pei1_NR,
                    pei2_f = pei2_NV/pei2_NR,
                    pei3_f = pei3_NV/pei3_NR,
                    usa1_f = usa1_NV/usa1_NR,
                    usa2_f = usa2_NV/usa2_NR,
                    usa3_f = usa3_NV/usa3_NR,
                    usa4_f = usa4_NV/usa4_NR,
                    usa5_f = usa5_NV/usa5_NR) %>%
                mutate(usa_f = (usa1_f+usa2_f+usa3_f+usa4_f+usa5_f)/5,
                       pei_f = (pei1_f+pei2_f+pei3_f)/3)
            return(table_output)
            #assign(paste0(region,"_",snvs,method,"_CN",cn,".bed"), table_output)
    }

    for(snv_file in c("allUSAnoPEInoH", "allCanyH", "allCnoH")){
                    print(snv_file)
        for(loh_file in c("","_LOH","_noLOH")){
                    print(paste(" ",loh_file))
            for(method_file in c("","_noH0-1")){
                    print(paste("  ",method_file))
                for(cn_file in c(0,1,2,3,4,5,6,"7plus")){
                #for(cn_file in c(0,1)){
                    print(paste("   ",cn_file))
                    assign(paste0("USA_",snv_file,loh_file,method_file,"_CN",cn_file), load_cn_snvs("USA", snv_file, loh_file, method_file, cn_file))
                }
            }
        }
    }
    for(snv_file in c("allPEInoUSAnoH", "allCanyH", "allCnoH")){
                    print(snv_file)
        for(loh_file in c("","_LOH","_noLOH")){
                    print(paste(" ",loh_file))
            for(method_file in c("","_noH0-1")){
                    print(paste("  ",method_file))
                for(cn_file in c(0,1,2,3,4,5,6,"7plus")){
                #for(cn_file in c(0,1)){
                    print(paste("   ",cn_file))
                    assign(paste0("PEI_",snv_file,loh_file,method_file,"_CN",cn_file), load_cn_snvs("PEI", snv_file, loh_file, method_file, cn_file))
                }
            }
        }
    }
    
# plot histograms (eventually overlay)
    pdf("CN_regions_overlay.pdf")
        for(snv_file in c("allUSAnoPEInoH", "allCanyH", "allCnoH")){
            print(snv_file)
            plot1 <- ggplot() +
                    geom_freqpoly(data=get(paste0("USA_",snv_file,"_CN2")), aes(usa_f), color="black", binwidth = 0.01) +
                    geom_freqpoly(data=get(paste0("USA_",snv_file,"_CN3")), aes(usa_f), color="black", binwidth = 0.01) +
                    geom_freqpoly(data=get(paste0("USA_",snv_file,"_CN4")), aes(usa_f), color="black", binwidth = 0.01) +
                    geom_freqpoly(data=get(paste0("USA_",snv_file,"_CN5")), aes(usa_f), color="black", binwidth = 0.01) +
                    geom_freqpoly(data=get(paste0("USA_",snv_file,"_CN6")), aes(usa_f), color="black", binwidth = 0.01) +
                    xlim(0,1.01)+
                    xlab("Allele freq")+
                    theme_classic() +
                    theme(axis.text=element_text(size=12,face="bold"),
                        axis.title=element_text(size=12,face="bold"),
                        text=element_text(size=12,face="bold"),
                        legend.title = element_blank())+
                    ggtitle(paste0("USA_",snv_file,"_CN_comparison"))
            print(plot1)
        }
        for(snv_file in c("allPEInoUSAnoH", "allCanyH", "allCnoH")){
            print(snv_file)
            plot1 <- ggplot() +
                    geom_freqpoly(data=get(paste0("PEI_",snv_file,"_CN2")), aes(pei_f), color="black", binwidth = 0.01) +
                    geom_freqpoly(data=get(paste0("PEI_",snv_file,"_CN3")), aes(pei_f), color="black", binwidth = 0.01) +
                    geom_freqpoly(data=get(paste0("PEI_",snv_file,"_CN4")), aes(pei_f), color="black", binwidth = 0.01) +
                    geom_freqpoly(data=get(paste0("PEI_",snv_file,"_CN5")), aes(pei_f), color="black", binwidth = 0.01) +
                    geom_freqpoly(data=get(paste0("PEI_",snv_file,"_CN6")), aes(pei_f), color="black", binwidth = 0.01) +
                    xlim(0,1.01)+
                    xlab("Allele freq")+
                    theme_classic() +
                    theme(axis.text=element_text(size=12,face="bold"),
                        axis.title=element_text(size=12,face="bold"),
                        text=element_text(size=12,face="bold"),
                        legend.title = element_blank())+
                    ggtitle(paste0("PEI_",snv_file,"_CN_comparison"))
            print(plot1)
        }
    dev.off()

