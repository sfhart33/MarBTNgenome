# originally here: C:\Users\shart\Metzger Lab Dropbox\Sam_data\Mya_genome\repeats\repeat_library_plotting.r

library(tidyverse)
setwd("/ssd3/Mar_genome_analysis/repeat_lib")
# merge into one table of depths
    for(i in c("MELC-A11_S1", "PEI-DN08_S3", "NYTC-C9_S2", "MELC-2E11", "MELC-A9", "PEI-DF490", "DF-488", "DN-HL03", "FFM-19G1", "FFM-20B2", "FFM-22F10")){
        table <- read.delim(paste("repeatLIB", i, "depth", sep="."), header = TRUE)
        if(i == "MELC-A11_S1"){coverage <- table}
        if(i != "MELC-A11_S1"){coverage <- full_join(coverage, table)}
    }
    head(coverage)
    tail(coverage)
    write.table(coverage, "all_samples_repeatLB.coverage", quote = FALSE, row.names=FALSE)


# transform table and calculate average coverages
    coverage2 <- coverage %>%
        separate(col = repeat_name, sep = "#", into = c("name","library","type"), remove = TRUE) %>%
        separate(col = type, sep = "/", into = c("class","subclass"), fill= "right", remove = TRUE) %>%            
        mutate(
            USA_avg = (FFM.19G1 + FFM.20B2 + FFM.22F10 + MELC.A11_S1 + NYTC.C9_S2)/5,
            PEI_avg = (PEI.DN08_S3 + DF.488 + DN.HL03)/3,
            C_avg = (FFM.19G1 + FFM.20B2 + FFM.22F10 + MELC.A11_S1 + NYTC.C9_S2 + PEI.DN08_S3 + DF.488 + DN.HL03)/8,
            H_avg = (MELC.2E11 + MELC.A9 + PEI.DF490)/3,
            C_H_ttest = 1,
            C_Hmax_ttest = 1,
            P_U_ttest = 1,
            C_Hmax_ratio = 1)
    for(i in 1:nrow(coverage2)){
        coverage2[i,"C_H_ttest"]<- t.test(coverage2[i,c("FFM.19G1","FFM.20B2","FFM.22F10","MELC.A11_S1","NYTC.C9_S2","PEI.DN08_S3","DF.488","DN.HL03")],coverage2[i,c("MELC.2E11","MELC.A9","PEI.DF490")])$p.value
        coverage2[i,"C_Hmax_ttest"]<- t.test(coverage2[i,c("FFM.19G1","FFM.20B2","FFM.22F10","MELC.A11_S1","NYTC.C9_S2","PEI.DN08_S3","DF.488","DN.HL03")],mu=max(coverage2[i,c("MELC.2E11","MELC.A9","PEI.DF490")]))$p.value
        coverage2[i,"P_U_ttest"]<- t.test(coverage2[i,c("FFM.19G1","FFM.20B2","FFM.22F10","MELC.A11_S1","NYTC.C9_S2")],coverage2[i,c("PEI.DN08_S3","DF.488","DN.HL03")])$p.value
        coverage2[i,"C_Hmax_ratio"]<- sum(coverage2[i,c("FFM.19G1","FFM.20B2","FFM.22F10","MELC.A11_S1","NYTC.C9_S2","PEI.DN08_S3","DF.488","DN.HL03")])/8/max(coverage2[i,c("MELC.2E11","MELC.A9","PEI.DF490")])
    }
    head(coverage2)
    coverage2 <- coverage2 %>%
        mutate(
            C_H_ratio = C_avg/H_avg,
            P_U_ratio = PEI_avg/USA_avg,
            )
    coverage <- coverage2       
    head(coverage)
    
    
# Which elements are present in all samples but significantly more in cancer?
    all_samples <- filter(coverage, FFM.19G1>1, FFM.20B2>1, FFM.22F10>1, MELC.A11_S1>1, NYTC.C9_S2>1, PEI.DN08_S3>1, DF.488>1, DN.HL03>1, MELC.2E11>1, MELC.A9>1, PEI.DF490>1)
    all_samples$C_H_ttest_corr <- p.adjust(all_samples$C_H_ttest,method ="bonferroni",n=nrow(all_samples))
    all_samples$C_Hmax_ttest_corr <- p.adjust(all_samples$C_Hmax_ttest,method ="bonferroni",n=nrow(all_samples))
    all_samples$P_U_ttest_corr <- p.adjust(all_samples$P_U_ttest,method ="bonferroni",n=nrow(all_samples))
    all_samples$C_H_ttest_fdr <- p.adjust(all_samples$C_H_ttest,method ="fdr",n=nrow(all_samples))
    all_samples$P_U_ttest_fdr <- p.adjust(all_samples$P_U_ttest,method ="fdr",n=nrow(all_samples))

    for(i in c(0.001,0.01,0.05,0.1,0.2,0.5)){
        print(paste0("p<",i))
        filter(all_samples, C_H_ttest_corr < i,C_H_ratio > 5) %>% nrow() %>% print()
        filter(all_samples, C_Hmax_ttest_corr < i,C_Hmax_ratio > 5) %>% nrow() %>% print()
    }
    significant <- filter(all_samples, C_H_ttest_corr < 0.05)
        nrow(significant) # 136
    significant5x <- filter(significant, C_H_ratio > 5)
        nrow(significant5x) # 45
    # significant_fdr <- filter(all_samples, C_H_ttest_fdr < 0.05)
    # nrow(significant_fdr)
    # significant5x_fdr <- filter(significant_fdr, C_H_ratio > 5)
    # nrow(significant5x_fdr)
    # significant5x_fdrDNA <- filter(significant5x_fdr, class == "DNA")
    # nrow(significant5x_fdrDNA)
    # significant10x_fdr <- filter(significant_fdr, C_H_ratio > 10)
    # nrow(significant10x_fdr)
    # significant10x_fdrDNA <- filter(significant10x_fdr, class == "DNA")
    # nrow(significant10x_fdrDNA)
    steamer <- all_samples %>% filter(name== "NEW_CONTIG_MERGE_326") %>% print()
    significant5xDNA <- filter(significant5x, class == "DNA")
# Enrichment of DNA elements
        nrow(significant5xDNA) # 8
        nrow(significant5x) # 45
        
        DNA <- filter(all_samples, class == "DNA")
        nrow(DNA) # 171
        nrow(all_samples) # 4471
        #chisqared in excel: 1.06e-6
    # Holds up using other filtering methods
        filter(all_samples, C_H_ttest_corr < 0.05,C_H_ratio > 10, class == "DNA") %>% nrow() # 4
        filter(all_samples, C_H_ttest_corr < 0.05,C_H_ratio > 10) %>% nrow() # 23
        filter(all_samples, C_H_ttest_fdr < 0.05,C_H_ratio > 10, class == "DNA") %>% nrow() # 8
        filter(all_samples, C_H_ttest_fdr < 0.05,C_H_ratio > 10) %>% nrow() # 38    
        filter(all_samples, C_H_ttest_fdr < 0.05,C_H_ratio > 5, class == "DNA") %>% nrow() # 15
        filter(all_samples, C_H_ttest_fdr < 0.05,C_H_ratio > 5) %>% nrow() # 81

# PEI vs USA
    sublineage_diff <- filter(all_samples, P_U_ttest_corr < 0.05)
    filter(sublineage_diff, P_U_ratio < 1) %>% nrow() %>% print() # 36
    filter(sublineage_diff, P_U_ratio > 1) %>% nrow() %>% print() # 20
        #chisqared in excel: 0.032
    filter(sublineage_diff, P_U_ratio < 1, class=="LTR") %>% nrow() %>% print() # 8
    filter(sublineage_diff, P_U_ratio > 1, class=="LTR") %>% nrow() %>% print() # 0
        #chisqared in excel: 0.0047
    filter(sublineage_diff, P_U_ratio < 1, class=="DNA") %>% nrow() %>% print() # 5
    filter(sublineage_diff, P_U_ratio > 1, class=="DNA") %>% nrow() %>% print() # 1
        #chisqared in excel: 0.10

# Sizes of repeats
    pdf("repeats_sizes.pdf")
    ggplot(all_samples)+
        geom_histogram(aes(log10(length)), bins=100, color = "black") +  
        theme_bw() +
        theme(axis.text=element_text(size=12,face="bold"),
                axis.title=element_text(size=14,face="bold"),
                text=element_text(size=14,face="bold")) +
        xlab("Repeat lengths (log10)") +
        ggtitle("All repeats at 1 copy or more in all clams")
    ggplot(all_samples)+
        geom_histogram(aes(length), bins=100, color = "black") +  
        theme_bw() +
        theme(axis.text=element_text(size=12,face="bold"),
                axis.title=element_text(size=14,face="bold"),
                text=element_text(size=14,face="bold")) +
        xlab("Repeat lengths") +
        ggtitle("All repeats at 1 copy or more in all clams")
    ggplot(significant)+
        geom_histogram(aes(log10(length)), bins=100, color = "black") +  
        theme_bw() +
        theme(axis.text=element_text(size=12,face="bold"),
                axis.title=element_text(size=14,face="bold"),
                text=element_text(size=14,face="bold")) +
        xlab("Repeat lengths (log10)") +
        ggtitle("Significantly expanded repeats (any amount)")
    ggplot(significant5x)+
        geom_histogram(aes(log10(length)), bins=100, color = "black") +  
        theme_bw() +
        theme(axis.text=element_text(size=12,face="bold"),
                axis.title=element_text(size=14,face="bold"),
                text=element_text(size=14,face="bold")) +
        xlab("Repeat lengths (log10)") +
        ggtitle("5x significantly expanded repeats")
    dev.off()

# plot PEI vs PEI and USA vs USA to check how similar they are
    # pdf("scatterplot_comparisons_repeats.pdf") #, width=4, height=4)
    # for(A in c("MELC.2E11", "MELC.A9", "PEI.DF490", "DF.488", "DN.HL03", "PEI.DN08_S3" ,"MELC.A11_S1", "NYTC.C9_S2", "FFM.19G1", "FFM.20B2", "FFM.22F10")){
    #     for(B in c("MELC.2E11", "MELC.A9", "PEI.DF490", "DF.488", "DN.HL03", "PEI.DN08_S3" ,"MELC.A11_S1", "NYTC.C9_S2", "FFM.19G1", "FFM.20B2", "FFM.22F10")){
    #         if(A != B){
    #             ploidy_regression <- lm(get(A) ~ get(B), data=all_samples)
    #             plot1 <- ggplot(all_samples, aes(get(A),get(B))) +
    #                 geom_point(size = 0.1, alpha = 0.1) +
    #                 theme_classic() +
    #                 theme(axis.text=element_text(size=12,face="bold"),
    #                         axis.title=element_text(size=16,face="bold"),
    #                         text=element_text(size=18,face="bold")) +
    #                 ggtitle(paste(A,"vs",B,"\nr.squared:",round(summary(ploidy_regression)$r.squared,3)))+
    #                 xlim(0,500)+
    #                 ylim(0,500)+
    #                 xlab(A)+
    #                 ylab(B) 
    #             print(plot1)
    #         }
    #     }  
    # }
    # dev.off()

#plots
    pdf("repeat_expansion_plots_final.pdf")
# main plot
    ggplot()+
        geom_point(data = all_samples, aes(log10(C_H_ratio), log10(C_H_ttest), color = "black"), color = "black") +  
        geom_point(data = significant5x, aes(log10(C_H_ratio), log10(C_H_ttest)), color = "black", size =3) +
        geom_point(data = steamer, aes(log10(C_H_ratio), log10(C_H_ttest)), color = "dark green", size =3) +
        geom_point(data = DNA, aes(log10(C_H_ratio), log10(C_H_ttest)), color = "blue") +
        geom_point(data = significant5xDNA,  aes(log10(C_H_ratio), log10(C_H_ttest)), color = "blue", size =3) +
        geom_text(aes(x=1.94, y=-6.35, label = "Steamer"), color = "dark green", size = 4,fontface = "bold") +
        geom_text(aes(x=2.03, y=-8.35, label = "DNA Transposons"), color = "blue", size = 4,fontface = "bold") +
        geom_text(aes(x=-2.1, y=-5.1, label = "corrected p=0.05"), color = "black", size = 4,fontface = "bold") +
        geom_text(aes(x=0.8, y=0, label = "5x"), color = "black", size = 4,fontface = "bold") +
        geom_text(aes(x=-0.8, y=0, label = "5x"), color = "black", size = 4,fontface = "bold") +
        theme_bw() +
        theme(axis.text=element_text(size=12,face="bold"),
                axis.title=element_text(size=14,face="bold"),
                text=element_text(size=14,face="bold")) +
        ylim(0, -11) +
        xlim(-2.5, 2.5) +
        xlab("Ratio of cancer (n=8) to healthy (n-3) repeat element copies") +
        ylab("p value (log10)") +
        ggtitle("Evidence for somatic repeat expansions \n(bonferroni correction)")+
        geom_hline(yintercept = log10(0.05/nrow(all_samples)), linetype="dashed") +
        geom_vline(xintercept = log10(5), linetype="dashed") +
        geom_vline(xintercept = log10(0.2), linetype="dashed") #+
        #geom_vline(xintercept = log10(2), linetype="dashed") #, size=1.5)
# plot with FDR correction rather than bonferroni
    ggplot()+
        geom_point(data = all_samples, aes(log10(C_H_ratio), log10(C_H_ttest_fdr), color = "black"), color = "black") +  
        geom_point(data = significant5x, aes(log10(C_H_ratio), log10(C_H_ttest_fdr)), color = "black", size =3) +
        geom_point(data = steamer, aes(log10(C_H_ratio), log10(C_H_ttest_fdr)), color = "dark green", size =3) +
        geom_point(data = DNA, aes(log10(C_H_ratio), log10(C_H_ttest_fdr)), color = "blue") +
        geom_point(data = significant5xDNA,  aes(log10(C_H_ratio), log10(C_H_ttest_fdr)), color = "blue", size =3) +
        #geom_text(aes(x=1.94, y=-6.35, label = "Steamer"), color = "dark green", size = 4,fontface = "bold") +
        #geom_text(aes(x=2.03, y=-8.35, label = "DNA Transposons"), color = "blue", size = 4,fontface = "bold") +
        geom_text(aes(x=-2.1, y=log10(0.04), label = "p=0.05"), color = "black", size = 4,fontface = "bold") +
        #geom_text(aes(x=0.8, y=0, label = "5x"), color = "black", size = 4,fontface = "bold") +
        #geom_text(aes(x=-0.8, y=0, label = "5x"), color = "black", size = 4,fontface = "bold") +
        theme_bw() +
        theme(axis.text=element_text(size=12,face="bold"),
                axis.title=element_text(size=14,face="bold"),
                text=element_text(size=14,face="bold")) +
        ylim(0, -8) +
        xlim(-2.5, 2.5) +
        xlab("Ratio of cancer (n=8) to healthy (n-3) repeat element copies") +
        ylab("p value (log10,fdr-corrected)") +
        ggtitle("Evidence for somatic repeat expansions \n(fdr correction)")+
        geom_hline(yintercept = log10(0.05), linetype="dashed") +
        geom_vline(xintercept = log10(5), linetype="dashed") +
        geom_vline(xintercept = log10(0.2), linetype="dashed") #+
        #geom_vline(xintercept = log10(2), linetype="dashed") #, size=1.5)

    for(elements in c("DNA","LTR", "LINE","RC","rRNA","Simple_repeat","SINE","snRNA","tRNA")){
        plotting <- filter(all_samples, class == elements) %>%
            ggplot()+
            geom_point(aes(log10(C_H_ratio), log10(C_H_ttest), color = "black"), color = "black") +  
            theme_bw() +
            theme(axis.text=element_text(size=12,face="bold"),
                    axis.title=element_text(size=16,face="bold"),
                    text=element_text(size=18,face="bold")) +
            ylim(0, -11) +
            xlim(-2.5, 2.5) +
            xlab("Log10(BTN/healthy copies)") +
            ylab("Log10(p value)") +
            ggtitle(paste("Element type: ", elements, sep=""))+
            geom_hline(yintercept = log10(0.05/nrow(all_samples)), linetype="dashed") +
            geom_vline(xintercept = log10(5), linetype="dashed") +
            geom_vline(xintercept = log10(0.2), linetype="dashed")
        print(plotting)
    }
    ggplot(all_samples)+
        geom_point(data =all_samples, aes(log10(P_U_ratio), log10(P_U_ttest), color = "black"), color = "black") +  
        geom_point(data = steamer, aes(log10(P_U_ratio), log10(P_U_ttest)), color = "dark green", size =3) +
        geom_point(data = DNA, aes(log10(P_U_ratio), log10(P_U_ttest)), color = "blue") +
        geom_text(aes(x=-2.1, y=-5.1, label = "Significance"), color = "black", size = 4,fontface = "bold") +
        geom_text(aes(x=log10(5), y=0, label = "5x"), color = "black", size = 4,fontface = "bold") +
        geom_text(aes(x=log10(0.2), y=0, label = "5x"), color = "black", size = 4,fontface = "bold") +
        theme_bw() +
        theme(axis.text=element_text(size=12,face="bold"),
                axis.title=element_text(size=14,face="bold"),
                text=element_text(size=14,face="bold")) +
        ylim(0, -10) +
        xlim(-2.5, 2.5) +
        xlab("Log10(PEI/USA copies)") +
        ylab("Log10(p value)") +
        ggtitle("All repeats: PEI vs USA")+
        geom_hline(yintercept = log10(0.05/nrow(all_samples)), linetype="dashed") +
        geom_vline(xintercept = log10(5), linetype="dashed") +
        geom_vline(xintercept = log10(0.2), linetype="dashed")
    for(elements in c("DNA","LTR", "LINE","RC","rRNA","Simple_repeat","SINE","snRNA","tRNA")){
        plotting <- filter(all_samples, class == elements) %>%
            ggplot()+
            geom_point(aes(log10(P_U_ratio), log10(P_U_ttest), color = "black"), color = "black") +  
            theme_bw() +
            theme(axis.text=element_text(size=12,face="bold"),
                    axis.title=element_text(size=16,face="bold"),
                    text=element_text(size=18,face="bold")) +
            ylim(0, -10) +
            xlim(-2.5, 2.5) +
            xlab("Log10(PEI/USA copies)") +
            ylab("Log10(p value)") +
            ggtitle(paste("Element type: ", elements, sep=""))+
            geom_hline(yintercept = log10(0.05/nrow(all_samples)), linetype="dashed") +
            geom_vline(xintercept = log10(5), linetype="dashed") +
            geom_vline(xintercept = log10(0.2), linetype="dashed") #+
        print(plotting)
    }
    dev.off()
