# Original file here: C:\Users\shart\Metzger Lab Dropbox\Sam_data\Mya_genome\SVs\delly_downstream.R

library(tidyverse)
setwd("/ssd3/Mar_genome_analysis/delly/final_run")

# load data
    #files <- list.files(pattern="*.table")
    files <- list.files(pattern="*.delly.table")[c(1,2,3,4,5,7,8,9,11,13,14)]

    count=0
    for(file in files){
        count=count+1
        print(count)
        if(count==1){
            mutations_list <- read.table(file, header = TRUE)
        }
        else{
            mutations_list2 <- read.table(file, header = TRUE)
            mutations_list <- full_join(mutations_list,mutations_list2)
        }
    }
    head(mutations_list)
    mutations_list2 <- mutations_list
    mutations_list[is.na(mutations_list)] <- 0
    nrow(mutations_list) # 237030
    cancer_list <- filter(mutations_list, MELC.A9_NV == 0, MELC.2E11_NV == 0, PEI.DF490_NV == 0)
    nrow(cancer_list) # 85505
    colnames(cancer_list)
    samples_list <- c("MELC.2E11","MELC.A9","PEI.DF490","PEI.DF488","PEI.DN03","PEI.DN08_S3","FFM.22F10","FFM.19G1","FFM.20B2","MELC.A11_S1","NYTC.C9_S2")
    samples_cancer <- c("PEI.DF488","PEI.DN03","PEI.DN08_S3","FFM.22F10","FFM.19G1","FFM.20B2","MELC.A11_S1","NYTC.C9_S2")

    samples_NR <- paste0(samples_list,"_NR")
    mean_depth <- colMeans(mutations_list2[,paste0(samples_list,"_NR")],  na.rm = TRUE) %>% print()

    for(sample in samples_list){
        #print(sample)
        depth <- mean_depth[paste0(sample,"_NR")]
        assign(paste0(sample,"_avg"), as.numeric(depth))
    }
    for(sample in samples_list){
        print(sample)
        print(get(paste0(sample,"_avg")))
    }

#plot allele counts
pdf("SV_allele_counts.pdf")
    for(i in samples_list){
    test <- filter(mutations_list,get(paste0(i,"_NV"))>0) #%>% #nrow() %>% print()
        plot1 <- ggplot(test, aes(get(paste0(i,"_NV")))) +
            geom_histogram(binwidth=1)+
            xlim(0,250)+
            xlab("Allele count")+
            theme_classic() +
                theme(axis.text=element_text(size=12,face="bold"),
                axis.title=element_text(size=16,face="bold"),
                text=element_text(size=18,face="bold")) +
            ggtitle(i)
        print(plot1)
    }
    for(i in samples_list){
        plot1 <- ggplot(test, aes(get(paste0(i,"_NV"))/get(paste0(i,"_NR")))) +
            geom_histogram(binwidth=0.01)+
            xlim(0,1.01)+
            xlab("Allele count")+
            theme_classic() +
                theme(axis.text=element_text(size=12,face="bold"),
                axis.title=element_text(size=16,face="bold"),
                text=element_text(size=18,face="bold")) +
            ggtitle(i)
        print(plot1)
    }
    dev.off()

# SV sizes histograms (all samples)
    all_dup <- filter(mutations_list2, TYPE == "DUP")
    all_del <- filter(mutations_list2, TYPE == "DEL")
    all_ins <- filter(mutations_list2, TYPE == "INS")
    all_inv <- filter(mutations_list2, TYPE == "INV")
    all_bnd <- filter(mutations_list2, TYPE == "BND")
    mutations_list_grouped <- group_by(mutations_list2,TYPE)
    median(mutations_list_grouped$END - mutations_list_grouped$POS, na.rm=TRUE)
    pdf("SV_sizes.pdf")
    for(i in samples_list){
    test <- filter(all_dup,get(paste0(i,"_NV"))>0) #%>% #nrow() %>% print()
        plot1 <- ggplot(test, aes(log10(END - POS))) +
            geom_histogram(bins=1000)+
            xlim(0,8)+
            xlab("Duplication size (log10)")+
            theme_classic() +
                theme(axis.text=element_text(size=12,face="bold"),
                axis.title=element_text(size=16,face="bold"),
                text=element_text(size=18,face="bold")) +
            ggtitle(i)
        print(plot1)
    }
    for(i in samples_list){
    test <- filter(all_del,get(paste0(i,"_NV"))>0) #%>% #nrow() %>% print()
        plot1 <- ggplot(test, aes(log10(END - POS))) +
            geom_histogram(bins=1000)+
            xlim(0,8)+
            xlab("Deletion size (log10)")+
            theme_classic() +
                theme(axis.text=element_text(size=12,face="bold"),
                axis.title=element_text(size=16,face="bold"),
                text=element_text(size=18,face="bold")) +
            ggtitle(i)
        print(plot1)
    }
    for(i in samples_list){
    test <- filter(all_inv,get(paste0(i,"_NV"))>0) #%>% #nrow() %>% print()
        plot1 <- ggplot(test, aes(log10(END - POS))) +
            geom_histogram(bins=1000)+
            xlim(0,8)+
            xlab("Inversion size (log10)")+
            theme_classic() +
                theme(axis.text=element_text(size=12,face="bold"),
                        axis.title=element_text(size=16,face="bold"),
                text=element_text(size=18,face="bold")) +
            ggtitle(i)
        print(plot1)
    }
    for(i in samples_list){
    test <- filter(all_ins,get(paste0(i,"_NV"))>0) #%>% #nrow() %>% print()
        plot1 <- ggplot(test, aes(INS_LEN)) +
            geom_histogram(bins=150)+
            xlim(0,150)+
            xlab("Insertion size")+
            theme_classic() +
                theme(axis.text=element_text(size=12,face="bold"),
                axis.title=element_text(size=16,face="bold"),
                text=element_text(size=18,face="bold")) +
            ggtitle(i)
        print(plot1)
    }
    dev.off()


# SV sizes histograms (just cancer)
    all_dup <- filter(cancer_list, TYPE == "DUP")
    all_del <- filter(cancer_list, TYPE == "DEL")
    all_ins <- filter(cancer_list, TYPE == "INS")
    all_inv <- filter(cancer_list, TYPE == "INV")
    all_bnd <- filter(cancer_list, TYPE == "BND")
    pdf("SV_sizes_cancer_only.pdf")
    for(i in samples_list){
    test <- filter(all_dup,get(paste0(i,"_NV"))>0) #%>% #nrow() %>% print()
        plot1 <- ggplot(test, aes(log10(END - POS))) +
            geom_histogram(bins=1000)+
            xlim(0,8)+
            xlab("Duplication size (log10)")+
            theme_classic() +
                theme(axis.text=element_text(size=12,face="bold"),
                axis.title=element_text(size=16,face="bold"),
                text=element_text(size=18,face="bold")) +
            ggtitle(i)
        print(plot1)
    }
    for(i in samples_list){
    test <- filter(all_del,get(paste0(i,"_NV"))>0) #%>% #nrow() %>% print()
        plot1 <- ggplot(test, aes(log10(END - POS))) +
            geom_histogram(bins=1000)+
            xlim(0,8)+
            xlab("Deletion size (log10)")+
            theme_classic() +
                theme(axis.text=element_text(size=12,face="bold"),
                axis.title=element_text(size=16,face="bold"),
                text=element_text(size=18,face="bold")) +
            ggtitle(i)
        print(plot1)
    }
    for(i in samples_list){
    test <- filter(all_inv,get(paste0(i,"_NV"))>0) #%>% #nrow() %>% print()
        plot1 <- ggplot(test, aes(log10(END - POS))) +
            geom_histogram(bins=1000)+
            xlim(0,8)+
            xlab("Inversion size (log10)")+
            theme_classic() +
                theme(axis.text=element_text(size=12,face="bold"),
                axis.title=element_text(size=16,face="bold"),
                text=element_text(size=18,face="bold")) +
            ggtitle(i)
        print(plot1)
    }
    for(i in samples_list){
    test <- filter(all_ins,get(paste0(i,"_NV"))>0) #%>% #nrow() %>% print()
        plot1 <- ggplot(test, aes(INS_LEN)) +
            geom_histogram(bins=150)+
            xlim(0,150)+
            xlab("Insertion size")+
            theme_classic() +
                theme(axis.text=element_text(size=12,face="bold"),
                axis.title=element_text(size=16,face="bold"),
                text=element_text(size=18,face="bold")) +
            ggtitle(i)
        print(plot1)
    }
    dev.off()

# count numbers of each that fit criteria
    threshold = 16
    samples_usa <- c("FFM.22F10","FFM.19G1","FFM.20B2","MELC.A11_S1","NYTC.C9_S2")
    samples_pei <- c("PEI.DF488","PEI.DN03","PEI.DN08_S3")

    noUSA_SV <- filter(cancer_list, #PEI.DF488_NV > PEI.DF488_avg/threshold &
                                    #PEI.DN03_NV > PEI.DN03_avg/threshold &
                                    #PEI.DN08_S3_NV > PEI.DN08_S3_avg/threshold &
                                    FFM.22F10_NV == 0 &
                                    FFM.19G1_NV == 0  &
                                    FFM.20B2_NV == 0  &
                                    MELC.A11_S1_NV == 0  &
                                    NYTC.C9_S2_NV == 0 )
    noUSA_SV2 <- filter(cancer_list, #PEI.DF488_NV > PEI.DF488_avg/threshold &
                                    #PEI.DN03_NV > PEI.DN03_avg/threshold &
                                    #PEI.DN08_S3_NV > PEI.DN08_S3_avg/threshold &
                                    #FFM.22F10_NV == 0 &
                                    FFM.19G1_NV == 0  &
                                    FFM.20B2_NV == 0  &
                                    #MELC.A11_S1_NV == 0  &
                                    NYTC.C9_S2_NV == 0 )
    noPEI_SV <- filter(cancer_list, PEI.DF488_NV == 0 &
                                    PEI.DN03_NV == 0  &
                                    PEI.DN08_S3_NV == 0  #&
                                    #FFM.22F10_NV > FFM.22F10_avg/threshold &
                                    #FFM.19G1_NV > FFM.19G1_avg/threshold &
                                    #FFM.20B2_NV > FFM.20B2_avg/threshold &
                                    #MELC.A11_S1_NV > MELC.A11_S1_avg/threshold &
                                    #NYTC.C9_S2_NV > NYTC.C9_S2_avg/threshold
                                    )
    allUSAnoPEI_SV <- filter(noPEI_SV, 
                                    FFM.22F10_NV > 0 &
                                    FFM.19G1_NV > 0 &
                                    FFM.20B2_NV > 0 &
                                    MELC.A11_S1_NV > 0 &
                                    NYTC.C9_S2_NV > 0 )
    allPEInoUSA_SV <- filter(noUSA_SV, PEI.DF488_NV > 0 &
                                    PEI.DN03_NV > 0 &
                                    PEI.DN08_S3_NV > 0)
    saveRDS(noUSA_SV, "SVs_noUSA.rds")
    saveRDS(noPEI_SV, "SVs_noPEI.rds")

    # counting SVs                  
        for(type in c("DEL","INS","DUP","INV","BND")){
            print(type)
            for(i in samples_cancer){
                print(paste0("   ",i))
                #filter(cancer_list, TYPE == type, get(paste0(i,"_NV")) > get(paste0(i,"_avg"))/threshold) %>% nrow() %>% print()
                #filter(cancer_list, TYPE == type, get(paste0(i,"_NV")) > 5) %>% nrow() %>% print()
                filter(cancer_list, TYPE == type, get(paste0(i,"_NV")) > 1) %>% nrow() %>% print()
            }
        }

# nonreferce SVs
    samples_nonref_list <- c("MELC.A9","PEI.DF490","PEI.DF488","PEI.DN03","PEI.DN08_S3","FFM.22F10","FFM.19G1","FFM.20B2","MELC.A11_S1","NYTC.C9_S2")
    SVs_nonref <- filter(mutations_list, MELC.2E11_NV == 0)
    saveRDS(SVs_nonref, "SVs_nonref.rds")
# Overlaid histograms of SV sizes
    pdf("SV_sizes_nonref_overlay.pdf")
    for(type in c("DEL","DUP","INV")){
        test <- filter(SVs_nonref, TYPE == type) #%>% #nrow() %>% print()
        plot1 <- ggplot(test) +
            geom_freqpoly(data = filter(test,MELC.A9_NV>0), aes(log10(END - POS)), binwidth=0.2, color = "green")+
            geom_freqpoly(data = filter(test,PEI.DF490_NV>0), aes(log10(END - POS)), binwidth=0.2, color = "green")+
            geom_freqpoly(data = filter(test,PEI.DF488_NV>0), aes(log10(END - POS)), binwidth=0.2, color = "red")+
            geom_freqpoly(data = filter(test,PEI.DN03_NV>0), aes(log10(END - POS)), binwidth=0.2, color = "red")+
            geom_freqpoly(data = filter(test,PEI.DN08_S3_NV>0), aes(log10(END - POS)), binwidth=0.2, color = "red")+
            geom_freqpoly(data = filter(test,FFM.19G1_NV>0), aes(log10(END - POS)), binwidth=0.2, color = "blue")+
            geom_freqpoly(data = filter(test,FFM.20B2_NV>0), aes(log10(END - POS)), binwidth=0.2, color = "blue")+
            geom_freqpoly(data = filter(test,FFM.22F10_NV>0), aes(log10(END - POS)), binwidth=0.2, color = "blue")+
            geom_freqpoly(data = filter(test,MELC.A11_S1_NV>0), aes(log10(END - POS)), binwidth=0.2, color = "blue")+
            geom_freqpoly(data = filter(test,NYTC.C9_S2_NV>0), aes(log10(END - POS)), binwidth=0.2, color = "blue")+
            xlim(1,8)+
            xlab(paste(type," size (log10)"))+
            theme_classic() +
                theme(axis.text=element_text(size=12,face="bold"),
                axis.title=element_text(size=16,face="bold"),
                text=element_text(size=18,face="bold")) +
            ggtitle(paste(type, "nonref"))
        print(plot1)
        usa_test <- filter(noPEI_SV, TYPE == type)
        pei_test <- filter(noUSA_SV, TYPE == type)
        plot1 <- ggplot(test) +
            geom_freqpoly(data = filter(pei_test,PEI.DF488_NV>0), aes(log10(END - POS)), binwidth=0.2, color = "red")+
            geom_freqpoly(data = filter(pei_test,PEI.DN03_NV>0), aes(log10(END - POS)), binwidth=0.2, color = "red")+
            geom_freqpoly(data = filter(pei_test,PEI.DN08_S3_NV>0), aes(log10(END - POS)), binwidth=0.2, color = "red")+
            geom_freqpoly(data = filter(usa_test,FFM.19G1_NV>0), aes(log10(END - POS)), binwidth=0.2, color = "blue")+
            geom_freqpoly(data = filter(usa_test,FFM.20B2_NV>0), aes(log10(END - POS)), binwidth=0.2, color = "blue")+
            geom_freqpoly(data = filter(usa_test,FFM.22F10_NV>0), aes(log10(END - POS)), binwidth=0.2, color = "blue")+
            geom_freqpoly(data = filter(usa_test,MELC.A11_S1_NV>0), aes(log10(END - POS)), binwidth=0.2, color = "blue")+
            geom_freqpoly(data = filter(usa_test,NYTC.C9_S2_NV>0), aes(log10(END - POS)), binwidth=0.2, color = "blue")+
            xlim(1,8)+
            xlab(paste(type," size (log10)"))+
            theme_classic() +
                theme(axis.text=element_text(size=12,face="bold"),
                axis.title=element_text(size=16,face="bold"),
                text=element_text(size=18,face="bold")) +
            ggtitle(paste(type, "somatic (not in healthy or other sublineage)"))
        print(plot1)
    }
    dev.off()
    SVs_nonref_dup_sizes <- filter(SVs_nonref, TYPE == "DUP") %>% mutate(SIZE = (END - POS), SIZE_LOG = log10(END - POS))
    wilcox.test(filter(SVs_nonref_dup_sizes, MELC.A9_NV>0 | PEI.DF490_NV>0)[,"SIZE"],
        filter(SVs_nonref_dup_sizes, PEI.DF488_NV>0 | PEI.DN03_NV>0 | PEI.DN08_S3_NV>0 | FFM.19G1_NV>0 | FFM.20B2_NV>0 | FFM.22F10_NV>0 | MELC.A11_S1_NV>0 | NYTC.C9_S2_NV>0)[,"SIZE"])
        # W = 11383807, p-value = 1.572e-11
        # alternative hypothesis: true location shift is not equal to 0

# size summary stats
    #nonref_SV_size <- data.frame(matrix(ncol=3,nrow=0, dimnames=list(NULL, c("sample", "sv", "size"))))
    threshold = 8
    count=0
    for(type in c("DEL","DUP","INV")){
        print(type)
        for(i in samples_nonref_list){
            counting <- filter(SVs_nonref, TYPE == type, get(paste0(i,"_NV")) > get(paste0(i,"_avg"))/threshold) %>%
                mutate(sample = i, size = END-POS) %>%
                select(sample, TYPE, size)
            print(sum(counting$size))
        if(count==0){
            nonref_SV_size <- counting
            count= count+1
        }
        else{
            nonref_SV_size <- rbind(nonref_SV_size,counting)
        }    
        }
    }
    nonref_SV_size$sample <- as.factor(nonref_SV_size$sample)
    pdf("size.boxplot2.pdf")
        ggplot(nonref_SV_size, aes(x=as.factor(paste0(TYPE,"_",sample)), y=log10(size))) +
            geom_boxplot(outlier.colour="black", outlier.size=1) +
            theme_classic() +
                theme(axis.text=element_text(size=12,face="bold"),
                axis.title=element_text(size=12,face="bold"),
                text=element_text(size=12,face="bold"),
                axis.title.x = element_blank(),
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                legend.position = "none")
    dev.off()
 #somatic size summary plots
    count=0
    for(type in c("DEL","DUP","INV")){
        print(type)
        for(i in samples_pei){
            counting <- filter(noUSA_SV, TYPE == type, get(paste0(i,"_NV")) > get(paste0(i,"_avg"))/threshold) %>%
                mutate(sample = i, size = END-POS, region="PEI") %>%
                select(sample, TYPE, size, region)
            print(sum(counting$size))
        if(count==0){
            somatic_SV_size <- counting
            count= count+1
        }
        else{
            somatic_SV_size <- rbind(somatic_SV_size,counting)
        }   
        }        
        for(i in samples_usa){
            counting <- filter(noPEI_SV, TYPE == type, get(paste0(i,"_NV")) > get(paste0(i,"_avg"))/threshold) %>%
                mutate(sample = i, size = END-POS, region="USA") %>%
                select(sample, TYPE, size, region)
            print(sum(counting$size))
            somatic_SV_size <- rbind(somatic_SV_size,counting)
        }   
    }
    pdf("size.boxplot_somatic.pdf")
        ggplot(somatic_SV_size, aes(x=as.factor(paste0(TYPE,"_",sample)), y=log10(size),color=region)) +
            geom_boxplot(outlier.colour="black", outlier.size=1) +
            theme_classic() +
                theme(axis.text=element_text(size=12,face="bold"),
                axis.title=element_text(size=12,face="bold"),
                text=element_text(size=12,face="bold"),
                axis.title.x = element_blank(),
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                legend.position = "none")
    dev.off()

# threshold test (DOESN'T MATTER MUCH)
    for(type in c("DEL","INS","DUP","INV","BND")){
        print(type)
        for(i in samples_nonref_list){
            counting <- filter(SVs_nonref, TYPE == type, get(paste0(i,"_NV")) > 1) %>% nrow()
            print(paste0("   ",i))
            print(paste0("          ",counting))
            for(threshold in c(4,6,8,10,12,14,16)){     
            print(paste0("      ",threshold))
            counting <- filter(SVs_nonref, TYPE == type, get(paste0(i,"_NV")) > get(paste0(i,"_avg"))/threshold) %>% nrow()
            print(paste0("          ",counting))
            }

        }
    }

############################################ SUMMARY STATS AND PLOTS #################################################################

# count somatic SVs
    somatic_SV_counts <- data.frame(matrix(ncol=6,nrow=0, dimnames=list(NULL, c("sample","region", "sv", "count","sub_shared","mean_depth"))))
    threshold = 8
    count = 0
    for(type in c("DEL","INS","DUP","INV","BND")){
        print(type)
        for(i in samples_pei){
            count = count + 1
            counting <- filter(noUSA_SV, TYPE == type, get(paste0(i,"_NV")) > get(paste0(i,"_avg"))/threshold) %>% nrow()
            #counting <- filter(noUSA_SV, TYPE == type, get(paste0(i,"_NV")) > 1) %>% nrow()
            counting_shared <- filter(allPEInoUSA_SV, TYPE == type) %>% nrow()
            print(paste0("   ",i))
            print(paste0("          ",counting))
            somatic_SV_counts[count,] <- c(i,"PEI",type,counting,counting_shared,get(paste0(i,"_avg")))
        }
        for(i in samples_usa){
            count = count + 1
            counting <- filter(noPEI_SV, TYPE == type, get(paste0(i,"_NV")) > get(paste0(i,"_avg"))/threshold) %>% nrow()
            #counting <- filter(noPEI_SV, TYPE == type, get(paste0(i,"_NV")) > 1) %>% nrow()
            counting_shared <- filter(allUSAnoPEI_SV, TYPE == type) %>% nrow()
            print(paste0("   ",i))
            print(paste0("          ",counting))
            somatic_SV_counts[count,] <- c(i,"USA",type,counting,counting_shared,get(paste0(i,"_avg")))
        }
    }

    somatic_SV_counts$sample <- as.factor(somatic_SV_counts$sample) 
    somatic_SV_counts$region <- as.factor(somatic_SV_counts$region) 
    somatic_SV_counts$count <- as.integer(somatic_SV_counts$count) 
    somatic_SV_counts$sub_shared <- as.numeric(somatic_SV_counts$sub_shared) 
    somatic_SV_counts$mean_depth <- as.numeric(somatic_SV_counts$mean_depth) 
    somatic_SV_counts_unchanged <- somatic_SV_counts
    somatic_SV_counts[1:8,"count"] <- somatic_SV_counts[1:8,"count"]/10
    somatic_SV_counts[1:8,"sv"] <- paste0(somatic_SV_counts[1:8,"sv"],"/10")
    somatic_SV_counts$sv <- as.factor(somatic_SV_counts$sv) 
    somatic_SV_counts

    somatic_SV_summary <- data.frame(matrix(ncol=5,nrow=0, dimnames=list(NULL, c("region", "sv", "mean","sd","sub_shared"))))
    count = 0
    for(type in c("DEL/10","INS","DUP","INV","BND")){
        for(sub in c("USA","PEI")){
            count = count + 1
            subset <- filter(somatic_SV_counts, sv==type, region==sub)
            subset_mean <- mean(subset$count)
            subset_sd <- sd(subset$count)
            subset_sem <- sd(subset$count)
            subset_sub_shared <- subset[1,"sub_shared"]
            somatic_SV_summary[count,] <- c(sub,type,subset_mean,subset_sd,subset_sub_shared)
        }
    }
    somatic_SV_summary$region <- as.factor(somatic_SV_summary$region) 
    somatic_SV_summary$sv <- as.factor(somatic_SV_summary$sv) 
    somatic_SV_summary$mean <- as.numeric(somatic_SV_summary$mean) 
    somatic_SV_summary$sd <- as.numeric(somatic_SV_summary$sd) 
    somatic_SV_summary$sub_shared <- as.numeric(somatic_SV_summary$sub_shared) 
    somatic_SV_summary

# count SVs (non-ref, including healthies)
    samples_nonref_list <- c("MELC.A9","PEI.DF490","PEI.DF488","PEI.DN03","PEI.DN08_S3","FFM.22F10","FFM.19G1","FFM.20B2","MELC.A11_S1","NYTC.C9_S2")
    SVs_nonref <- filter(mutations_list, MELC.2E11_NV == 0)
    SV_counts <- data.frame(matrix(ncol=5,nrow=0, dimnames=list(NULL, c("sample","region", "sv", "count","mean_depth"))))
    threshold = 8
    count = 0
    for(type in c("DEL","INS","DUP","INV","BND")){
        print(type)
        for(i in samples_nonref_list){
            count = count + 1
            counting <- filter(SVs_nonref, TYPE == type, get(paste0(i,"_NV")) > get(paste0(i,"_avg"))/threshold) %>% nrow()
            # counting <- filter(SVs_nonref, TYPE == type, get(paste0(i,"_NV")) > 1) %>% nrow()
            print(paste0("   ",i))
            print(paste0("          ",counting))
            SV_counts[count,] <- c(i,NA,type,counting,get(paste0(i,"_avg")))
        }
    }
    SV_counts$region <- rep(c("H","H","PEI","PEI","PEI","USA","USA","USA","USA","USA"),5)
    SV_counts$sample <- as.factor(SV_counts$sample) 
    SV_counts$region <- as.factor(SV_counts$region) 
    SV_counts$count <- as.integer(SV_counts$count) 
    SV_counts$mean_depth <- as.numeric(SV_counts$mean_depth) 
    SV_counts$sv <- as.factor(SV_counts$sv) 
    SV_counts


    SV_summary <- data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c("region", "sv", "mean","sd"))))
    count = 0
    for(type in c("DEL","INS","DUP","INV","BND")){
        for(sub in c("H","USA","PEI")){
            count = count + 1
            subset <- filter(SV_counts, sv==type, region==sub)
            subset_mean <- mean(subset$count)
            subset_sd <- sd(subset$count)
            subset_sem <- sd(subset$count)
            SV_summary[count,] <- c(sub,type,subset_mean,subset_sd)
        }
    }
    SV_summary$region <- as.factor(SV_summary$region) 
    SV_summary$sv <- as.factor(SV_summary$sv) 
    SV_summary$mean <- as.numeric(SV_summary$mean) 
    SV_summary$sd <- as.numeric(SV_summary$sd) 
    SV_summary

# Normalize to healthy
    SV_counts_norm <- SV_counts
    SV_counts_norm[1:10,4] <- SV_counts[1:10,4]/SV_summary[1,3]
    SV_counts_norm[11:20,4] <- SV_counts[11:20,4]/SV_summary[4,3]
    SV_counts_norm[21:30,4] <- SV_counts[21:30,4]/SV_summary[7,3]
    SV_counts_norm[31:40,4] <- SV_counts[31:40,4]/SV_summary[10,3]
    SV_counts_norm[41:50,4] <- SV_counts[41:50,4]/SV_summary[13,3]
    SV_summary_norm <- SV_summary
    SV_summary_norm[1:3,3:4] <- SV_summary[1:3,3:4]/SV_summary[1,3]
    SV_summary_norm[4:6,3:4] <- SV_summary[4:6,3:4]/SV_summary[4,3]
    SV_summary_norm[7:9,3:4] <- SV_summary[7:9,3:4]/SV_summary[7,3]
    SV_summary_norm[10:12,3:4] <- SV_summary[10:12,3:4]/SV_summary[10,3]
    SV_summary_norm[13:15,3:4] <- SV_summary[13:15,3:4]/SV_summary[13,3]
    SV_summary_norm[c(1,4,7,10,13),"sd"] <- c(NA,NA,NA,NA,NA)
# Plotting
    pdf("SV_summary_plot.pdf",width=9, height=7)
    ggplot(somatic_SV_summary, aes(x=sv,y=mean,fill=region)) +
    geom_bar(position="dodge", stat="identity")+
    scale_fill_manual(values=c("red","blue"))+ # c("grey","red", "blue")
    geom_point(data=somatic_SV_counts,
                aes(x=sv,y=count,fill=region,position=region),
                position=position_dodge(.9), fill="black", size=2)+
    geom_errorbar(aes(ymin=(mean-sd), ymax=(mean+sd)), width=.2,position=position_dodge(.9), color="grey46")+
    xlab("SV type")+
    ylab("Number of SVs")+
    theme_classic() +
    theme(axis.text=element_text(size=12,face="bold"),
            axis.title=element_text(size=16,face="bold"),
            text=element_text(size=18,face="bold"),
            legend.title = element_blank())+
        ggtitle("SV counts per sample \n (all SVs not found in healthys or other sub-lineage)")

    ggplot(SV_summary_norm, aes(x=sv,y=mean,fill=region)) +
    geom_bar(position="dodge", stat="identity")+
    scale_fill_manual(values=c("grey","red", "blue"))+ # 
    geom_point(data=SV_counts_norm,
                aes(x=sv,y=count,fill=region,position=region),
                position=position_dodge(.9), fill="black", size=2)+
    geom_errorbar(aes(ymin=(mean-sd), ymax=(mean+sd)), width=.2,position=position_dodge(.9), color="grey46")+
    geom_text(data=filter(SV_summary,region=="H"),aes(x=sv,y=-0.05,label=paste0("x",round(mean))))+
    xlab("SV type")+
    ylab("Number of SVs \n (normalized to non-ref healthys)")+
    theme_classic() +
    theme(axis.text=element_text(size=12,face="bold"),
            axis.title=element_text(size=16,face="bold"),
            text=element_text(size=18,face="bold"),
            legend.title = element_blank())+
        ggtitle("SV counts per sample (all SVs not found in reference)")
    dev.off()
# t-tests       
    for(type in c("DEL","INS","DUP","INV","BND")){
        print(type)
        ttest_output <- t.test(filter(SV_counts,sv==type,region=="H")[,"count"],filter(SV_counts,sv==type,(region=="USA"|region=="PEI"))[,"count"])
        print(paste0("   Nonref: H vs C comparison p.value: ",ttest_output$p.value))

        ttest_output <- t.test(filter(SV_counts,sv==type,region=="H")[,"count"],filter(SV_counts,sv==type,region=="USA")[,"count"])
        print(paste0("   Nonref: H vs USA comparison p.value: ",ttest_output$p.value))
        ttest_output <- t.test(filter(SV_counts,sv==type,region=="H")[,"count"],filter(SV_counts,sv==type,region=="PEI")[,"count"])
        print(paste0("   Nonref: H vs PEI comparison p.value: ",ttest_output$p.value))
        ttest_output <- t.test(filter(SV_counts,sv==type,region=="PEI")[,"count"],filter(SV_counts,sv==type,region=="USA")[,"count"])
        print(paste0("   Nonref: PEI vs USA comparison p.value: ",ttest_output$p.value))

        ttest_output <- t.test(filter(somatic_SV_counts_unchanged,sv==type,region=="PEI")[,"count"],filter(somatic_SV_counts_unchanged,sv==type,region=="USA")[,"count"])
        print(paste0("   Somatic: PEI vs USA comparison p.value: ",ttest_output$p.value))
    }
        # [1] "DEL"
        # [1] "   Nonref: H vs C comparison p.value: 0.0106873986162626"
        # [1] "   Nonref: H vs USA comparison p.value: 0.072858563101637"
        # [1] "   Nonref: H vs PEI comparison p.value: 0.10965162856565"
        # [1] "   Nonref: PEI vs USA comparison p.value: 0.556494340623865"
        # [1] "   Somatic: PEI vs USA comparison p.value: 0.890658549847798"
        # [1] "INS"
        # [1] "   Nonref: H vs C comparison p.value: 0.269707387805138"
        # [1] "   Nonref: H vs USA comparison p.value: 0.458846355265254"
        # [1] "   Nonref: H vs PEI comparison p.value: 0.525266140001964"
        # [1] "   Nonref: PEI vs USA comparison p.value: 0.937219857083149"
        # [1] "   Somatic: PEI vs USA comparison p.value: 0.199591237033406"
        # [1] "DUP"
        # [1] "   Nonref: H vs C comparison p.value: 2.78902034330921e-05"
        # [1] "   Nonref: H vs USA comparison p.value: 0.000116511458554244"
        # [1] "   Nonref: H vs PEI comparison p.value: 0.0139109756234978"
        # [1] "   Nonref: PEI vs USA comparison p.value: 0.0121126563005277"
        # [1] "   Somatic: PEI vs USA comparison p.value: 9.67573919614661e-05"
        # [1] "INV"
        # [1] "   Nonref: H vs C comparison p.value: 1.54028716948826e-07"
        # [1] "   Nonref: H vs USA comparison p.value: 9.70342606582611e-05"
        # [1] "   Nonref: H vs PEI comparison p.value: 0.00493067891099698"
        # [1] "   Nonref: PEI vs USA comparison p.value: 0.323523985667811"
        # [1] "   Somatic: PEI vs USA comparison p.value: 0.736253998530645"
        # [1] "BND"
        # [1] "   Nonref: H vs C comparison p.value: 0.000265126870665696"
        # [1] "   Nonref: H vs USA comparison p.value: 0.00552498442235707"
        # [1] "   Nonref: H vs PEI comparison p.value: 0.059357295456871"
        # [1] "   Nonref: PEI vs USA comparison p.value: 0.374813110050265"
        # [1] "   Somatic: PEI vs USA comparison p.value: 0.0530534894023736"
