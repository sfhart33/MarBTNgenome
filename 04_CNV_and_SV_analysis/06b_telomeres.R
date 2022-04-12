# original file here: C:\Users\shart\Metzger Lab Dropbox\Sam_data\Mya_genome\telomeres.R

library(tidyverse)

# PREVIOUS KMER COUNTING

    # setwd("/ssd3/Mar_genome_analysis/telomeres")

    # #read.delim("outputcount.txt", header = FALSE)
    # names <- scan("outputnames.txt", character())
    # TotalBasecount <- scan("TotalBaseCount.txt")
    # TTAGGGcount <- scan("TTAGGGcount.txt")
    # CCCTAAcount <- scan("CCCTAAcount.txt")
    # TTAGGGTTAGGGcount <- scan("TTAGGGTTAGGGcount.txt")
    # CCCTAACCCTAAcount <- scan("CCCTAACCCTAAcount.txt")


    # telomeres <- data.frame(name = names, bases = TotalBasecount, TTAGGG = TTAGGGcount, CCCTAA = CCCTAAcount, TTAGGGTTAGGG = TTAGGGTTAGGGcount, CCCTAACCCTAA = CCCTAACCCTAAcount) %>%
    #                 separate(name, c(NA,NA,NA,NA,NA,"name"), "/") %>%
    #                 separate(name, c("name",NA), "_001") %>%
    #                 separate(name, c("name","read"), "_R")

    # telomeres1 <- filter(telomeres, read == 1) %>%
    #     rename(TTAGGG_F = TTAGGG) %>%
    #     rename(CCCTAA_F = CCCTAA) %>%
    #     rename(TTAGGGTTAGGG_F = TTAGGGTTAGGG) %>%
    #     rename(CCCTAACCCTAA_F = CCCTAACCCTAA) %>%
    #     rename(bases_F = bases)
    # telomeres2 <- filter(telomeres, read == 2) #%>%
    #     #rename(TTAGGG_R = TTAGGG)

    # telomeres <- telomeres1
    # telomeres$TTAGGG_R <- telomeres2$TTAGGG
    # telomeres$CCCTAA_R <- telomeres2$CCCTAA
    # telomeres$TTAGGGTTAGGG_R <- telomeres2$TTAGGGTTAGGG
    # telomeres$CCCTAACCCTAA_R <- telomeres2$CCCTAACCCTAA
    # telomeres$bases_R <- telomeres2$bases
    # telomeres <- mutate(telomeres, TTAGGG = TTAGGG_F+TTAGGG_R,
    #                             CCCTAA = CCCTAA_F+CCCTAA_R,
    #                             TTAGGGTTAGGG = TTAGGGTTAGGG_F+TTAGGGTTAGGG_R,
    #                             CCCTAACCCTAA = CCCTAACCCTAA_F+CCCTAACCCTAA_R,
    #                             single = TTAGGG_F+TTAGGG_R+CCCTAA_F+CCCTAA_R,
    #                             double = TTAGGGTTAGGG_F+TTAGGGTTAGGG_R+CCCTAACCCTAA_F+CCCTAACCCTAA_R,
    #                             bases = bases_F+bases_R,  ) %>%
    #                 select(name, bases, single, double)
    # ref_single_norm = filter(telomeres, name=="MELC-2E11")[1,3]/filter(telomeres, name=="MELC-2E11")[1,2]
    # ref_double_norm = filter(telomeres, name=="MELC-2E11")[1,4]/filter(telomeres, name=="MELC-2E11")[1,2]
    # telomeres <- mutate(telomeres,
    #                     single_norm = single/bases/ref_single_norm,
    #                     double_norm = double/bases/ref_double_norm)      
    # telomeres$region <- factor(levels = c("Healthy","PEI","USA"), c("PEI","PEI","USA","USA","Healthy","USA","USA","Healthy","USA","Healthy","PEI","USA"))           

    # pdf("telomeres.pdf")
    # ggplot(telomeres, aes(x=name,y=single_norm,fill=region))+
    #     geom_col()+
    #     scale_fill_manual(values=c("black","red","blue"))+
    #     ylab("TTAGGG/CCCTAA 6mer count relative reference")+
    #     theme_classic() +
    #     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

    # ggplot(telomeres, aes(x=name,y=double_norm,fill=region))+
    #     geom_col()+
    #     scale_fill_manual(values=c("black","red","blue"))+
    #     ylab("TTAGGGTTAGGG/CCCTAACCCTAA 12mer count relative reference")+
    #     theme_classic() +
    #     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    # dev.off()

# Plot TelSeq output

setwd("/ssd3/Mar_genome_analysis/telomeres/outputs")
list.files(pattern = "*telseq")

bams <- list.files(pattern = "*telseq")[2:14]
telos <- read.delim("01.MELC-2E11.bam.telseq", skip = 2, header = TRUE) %>% mutate(sample = "01.MELC-2E11.bam.telseq") %>% select(sample, everything())
for(i in bams){
    telo <- read.delim(i, skip = 2, header = TRUE) %>% mutate(sample = i) %>% select(sample, everything())
    telos <- rbind(telos, telo)
}
telos <- telos %>%
    separate(sample, c(NA,"sample",NA,NA), "\\.") %>%
    select(-ReadGroup,-Library, -Sample, -X) %>%
    mutate(depth = Total/telos[1,"Total"])
telos$region <- factor(levels = c("Healthy","PEI","USA"), c("Healthy","Healthy","Healthy","PEI","PEI","PEI","PEI","USA","USA","USA","USA","USA","USA","USA"))  

telos 

setwd("/ssd3/Mar_genome_analysis/telomeres")
pdf("telseq_outputs.pdf")
ggplot(telos, aes(x=sample,y=LENGTH_ESTIMATE,fill=region))+
    geom_col()+
    scale_fill_manual(values=c("black","red","blue"))+
    ylab("Telseq estimated telomere length (kB)")+
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
for(i in 1:16){
    plot1 <- ggplot(telos, aes(x=sample,y=get(paste0("TEL",i))/depth,fill=region))+
        geom_col()+
        scale_fill_manual(values=c("black","red","blue"))+
        ylab(paste("Telseq estimated TTAGGG kmers (norm depth), count =", i))+
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    print(plot1)
}
dev.off()

pdf("telseq_output_final.pdf")
ggplot(telos[c(1:5,7:9,11,13,14),], aes(x=sample,y=LENGTH_ESTIMATE,fill=region))+
    geom_col()+
    scale_fill_manual(values=c("black","red","blue"))+
    ylab("Telseq estimated telomere length (kB)")+
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


    telos <- telos[c(-6,-10,-12),]
    telos$name <- factor(levels = c("H1","H2","H3","PEI1","PEI2","PEI3","USA1","USA2","USA3","USA4","USA5"), c("H1","H2","H3","PEI1","PEI2","PEI3","USA1","USA2","USA3","USA4","USA5"))    

    telos_summary <- data.frame(matrix(ncol=3,nrow=0, dimnames=list(NULL, c("region", "mean","sd"))))
    count = 0
    for(sub in c("Healthy","USA","PEI")){
        count = count + 1
        subset <- filter(telos, region==sub)
        subset_mean <- mean(subset$LENGTH_ESTIMATE)
        subset_sd <- sd(subset$LENGTH_ESTIMATE)
        subset_sem <- sd(subset$LENGTH_ESTIMATE)
        telos_summary[count,] <- c(sub,subset_mean,subset_sd)
    }
    telos_summary$region <- as.factor(telos_summary$region) 
    telos_summary$mean <- as.numeric(telos_summary$mean) 
    telos_summary$sd <- as.numeric(telos_summary$sd) 
    telos_summary

    pdf("telo_summary_plot.pdf",width=2, height=2)
    ggplot(telos_summary, aes(x=region,y=mean,fill=region)) +
        geom_bar(position="dodge", stat="identity")+
        scale_fill_manual(values=c("black","red", "blue"))+ # 
        geom_point(data=telos,
                    aes(x=region,y=LENGTH_ESTIMATE,fill=region,position=region),
                    position=position_dodge(.9), fill="grey10", color="grey10", size=1)+
        geom_errorbar(aes(ymin=(mean-sd), ymax=(mean+sd)), width=.2,position=position_dodge(.9), color="grey46")+
        #geom_text(data=filter(SV_summary,region=="H"),aes(x=sv,y=-0.05,label=paste0("x",round(mean))))+
        xlab(NA)+
        ylab("Telomere length (kB)")+
        theme_classic() +
        theme(axis.text=element_text(size=8,face="bold"),
                axis.title=element_text(size=8,face="bold"),
                text=element_text(size=8,face="bold"),
                legend.title = element_blank())
    dev.off()



    pdf("telo_summary_plot2.pdf",width=1.5, height=2)  
    ggplot(telos) + 
        geom_pointrange(aes(x=name,y= LENGTH_ESTIMATE,ymin=depth,ymax=depth,color=region), size = 1, fatten = 0.75)+
        ylab("Telomere length (kB)")+
        xlab(NULL)+
        #ylimscale_y_continuous+
        #scale_y_continuous(breaks=seq(1,11,by=1), limits=c(0,12))+
        scale_color_manual(values=c("black", "red", "blue"))+
        theme_classic() +
            theme(axis.text=element_text(size=8,face="bold"),
            axis.title=element_text(size=8,face="bold"),
            text=element_text(size=8,face="bold"),
            legend.position = "none") +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    dev.off()


########## PLOT WITH TISSUES
    # setwd("/ssd3/Mar_genome_analysis/telomeres/outputs")
    # list.files(pattern = "*telseq")

    # bams <- list.files(pattern = "*telseq")[2:22]
    # telos <- read.delim("01.MELC-2E11.bam.telseq", skip = 2, header = TRUE) %>% mutate(sample = "01.MELC-2E11.bam.telseq") %>% select(sample, everything())
    # for(i in bams){
    #     telo <- read.delim(i, skip = 2, header = TRUE) %>% mutate(sample = i) %>% select(sample, everything())
    #     telos <- rbind(telos, telo)
    # }
    # telos <- telos %>%
    #     separate(sample, c(NA,"sample",NA,NA), "\\.") %>%
    #     select(-ReadGroup,-Library, -Sample, -X) %>%
    #     mutate(depth = Total/telos[1,"Total"])
    # telos$region <- factor(levels = c("Healthy","PEI","USA","Tissue"), c("Healthy","Healthy","Healthy","PEI","PEI","PEI","PEI","USA","USA","USA","USA","USA","USA","USA","Tissue","Tissue","Tissue","Tissue","Tissue","Tissue","Tissue","Tissue"))    
    # telos 

    # setwd("/ssd3/Mar_genome_analysis/telomeres")
    # pdf("telseq_outputs2.pdf")
    # ggplot(telos, aes(x=sample,y=LENGTH_ESTIMATE,fill=region))+
    #     geom_col()+
    #     scale_fill_manual(values=c("black","red","blue", "grey"))+
    #     ylab("Telseq estimated telomere length (kB)")+
    #     theme_classic() +
    #     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    # # for(i in 1:16){
    # #     plot1 <- ggplot(telos, aes(x=sample,y=get(paste0("TEL",i))/depth,fill=region))+
    # #         geom_col()+
    # #         scale_fill_manual(values=c("black","red","blue"))+
    # #         ylab(paste("Telseq estimated TTAGGG kmers (norm depth), count =", i))+
    # #         theme_classic() +
    # #         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    # #     print(plot1)
    # # }
    # dev.off()
