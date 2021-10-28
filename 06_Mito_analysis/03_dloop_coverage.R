
# also found here: C:\Users\shart\Metzger Lab Dropbox\Sam_data\Mya_genome\mitochondria\mito_dloop_analysis.r
    library(tidyverse)
    setwd("/ssd3/Mar_genome_analysis/bwa_mapping/mito/all_samples/depth")

# Load data
    depth_files <- read.delim("MELC-A10.mtgenome.depth", head=FALSE) %>% select(V2)%>% rename(pos=V2)
    for(sample in c("DF-488", "DN-HL03","FFM-19G1", "FFM-20B2", "FFM-22F10", "MELC-2E11", "MELC-A11_S1", "MELC-A9", "NYTC-C9_S2", "PEI-DF490", "PEI-DN08_S3")){
        depth_file <- read.delim(paste0(sample,".mtgenome.depth"), head=FALSE) %>% select(V3) %>% rename(!!sample:=V3)
        depth_files <- cbind(depth_files, depth_file)
    }

# Divide dloop and non-dloop regions
    dloop <- depth_files %>% filter(pos > 12164 & pos < 12870)
    notdloop <- depth_files %>% filter(pos < 12164 | pos > 12870)

# Calculate max and mean coverages
    colMax <- function (colData) {
        apply(colData, MARGIN=c(2), max)
    }
    colMeans(dloop)
    colMeans(notdloop)
    colMax(dloop)
    colMax(notdloop)
    colMeans(dloop)/colMeans(notdloop)
    colMax(dloop)/colMax(notdloop)
    colMax(dloop)/colMeans(notdloop)

# Normalize to non-dloop region
    depth_files2 <- depth_files
    depth_files2[,2:12] <- as.data.frame((t(t(as.matrix(depth_files[,2:12]))/colMeans(notdloop[,2:12], na.rm=TRUE))))
    head(depth_files2)

    dloop2 <- depth_files2 %>% filter(pos > 12164 & pos < 12870)

    pdf("depth_norm.pdf", width=4, height=4)
    ggplot(depth_files2) +
        geom_line(aes(pos, get("DF-488")), color = "red")+
        geom_line(aes(pos, get("DN-HL03")), color = "red")+
        # geom_line(aes(pos, get("DN-HL07")), color = "red")+
        geom_line(aes(pos, get("PEI-DN08_S3")), color = "red")+
        geom_line(aes(pos, get("FFM-19G1")), color = "blue")+
        geom_line(aes(pos, get("FFM-20B2")), color = "blue")+
        # geom_line(aes(pos, get("FFM-22A10")), color = "blue")+
        geom_line(aes(pos, get("FFM-22F10")), color = "blue")+
        # geom_line(aes(pos, get("MELC-A10")), color = "blue")+
        geom_line(aes(pos, get("MELC-A11_S1")), color = "blue")+
        geom_line(aes(pos, get("NYTC-C9_S2")), color = "blue")+
        geom_line(aes(pos, get("MELC-2E11")), color = "black")+
        geom_line(aes(pos, get("MELC-A9")), color = "black")+
        geom_line(aes(pos, get("PEI-DF490")), color = "black")+
        geom_text(aes(13500, 3.5, label = "Healthy (n=3)"), size = 4, color = "black", hjust="left") +
        geom_text(aes(13500, 5.5, label = "PEI (n=3)"), size = 4, color = "red", hjust="left") +
        geom_text(aes(13500, 7.5, label = "USA (n=5)"), size = 4, color = "blue", hjust="left") + # , fontface = "bold"
        ylab("Read depth (normalized)")+
        xlab("Mitochonria position")+
            theme_classic() +
                theme(axis.text=element_text(size=8,face="bold"),
                axis.title=element_text(size=10,face="bold"),
                text=element_text(size=12,face="bold"))      
    # for(sample in c("DF-488", "DN-HL03", "DN-HL07", "FFM-19G1", "FFM-20B2", "FFM-22A10", "FFM-22F10", "MELC-2E11", "MELC-A10", "MELC-A11_S1", "MELC-A9", "NYTC-C9_S2", "PEI-DF490", "PEI-DN08_S3")){
    # plot1 <- ggplot(depth_files2) +
    # 	geom_line(aes(pos, get(sample)))+
    # 	ylab("read depth")+
    #         theme_classic() +
    #             theme(axis.text=element_text(size=12,face="bold"),
    #             axis.title=element_text(size=16,face="bold"),
    #             text=element_text(size=18,face="bold")) +
    #         ggtitle(paste("Mt genome coverage", sample))
    # print(plot1)
    # }

    dev.off()
 

# estimate copy number based on peak position, which lines up with integers for healthys
    dlooppeak <- depth_files2 %>% filter(pos > 12300 & pos < 12500) %>% colMeans()
    # dlooppeak <- depth_files2 %>% filter(pos > 12300 & pos < 12500) %>% colMax()
    # dlooppeak <- depth_files2 %>% colMax()
    dlooppeak <- dlooppeak[c("MELC-2E11","MELC-A9","PEI-DF490","DF-488","DN-HL03","PEI-DN08_S3","FFM-19G1","FFM-20B2","FFM-22F10","MELC-A11_S1","NYTC-C9_S2")]

    dlooppeak_norm <- dlooppeak/dlooppeak["MELC-2E11"]*3 #assuming reference clam is CN3 based of pbjelly gap filling
    dloopH <- dlooppeak_norm[c("MELC-2E11","MELC-A9","PEI-DF490")]
    dloopC <- dlooppeak_norm[c("DF-488","DN-HL03","PEI-DN08_S3","FFM-19G1","FFM-20B2","FFM-22F10","MELC-A11_S1","NYTC-C9_S2")]
    dloopPEI <- dlooppeak_norm[c("DF-488","DN-HL03","PEI-DN08_S3")]
    dloopUSA <- dlooppeak_norm[c("FFM-19G1","FFM-20B2","FFM-22F10","MELC-A11_S1","NYTC-C9_S2")]

# PEI vs USA
    t.test(dloopPEI,dloopUSA)$p.value #  0.000169876 p<0.001
# HEALTHY vs CANCER
    t.test(dloopH,dloopC)$p.value # 0.0002933025 p<0.001
# HEALTHY vs USA
    t.test(dloopH,dloopUSA)$p.value # 0.0002164924 p<0.001
# HEALTHY vs PEI
    t.test(dloopH,dloopPEI)$p.value # 0.01285444 p<0.05

dlooppeakD <- as.data.frame(dlooppeak_norm)
dlooppeakD$region <- as.factor(c("H","H","H","PEI","PEI","PEI","USA","USA","USA","USA","USA"))
dlooppeakD
dlooppeak_summary <- dlooppeakD[c(1,4,7),]
dlooppeak_summary[1,1] <- mean(dlooppeakD[c("MELC-2E11","MELC-A9","PEI-DF490"),1])
dlooppeak_summary[2,1] <- mean(dlooppeakD[c("DF-488","DN-HL03","PEI-DN08_S3"),1])
dlooppeak_summary[3,1] <- mean(dlooppeakD[c("FFM-19G1","FFM-20B2","FFM-22F10","MELC-A11_S1","NYTC-C9_S2"),1])
dlooppeak_summary$sd <- c(sd(dlooppeakD[c("MELC-2E11","MELC-A9","PEI-DF490"),1]),
                        sd(dlooppeakD[c("DF-488","DN-HL03","PEI-DN08_S3"),1]),
                        sd(dlooppeakD[c("FFM-19G1","FFM-20B2","FFM-22F10","MELC-A11_S1","NYTC-C9_S2"),1]))
dlooppeak_summary$sem <- c(sd(dlooppeakD[c("MELC-2E11","MELC-A9","PEI-DF490"),1])/sqrt(3),
                        sd(dlooppeakD[c("DF-488","DN-HL03","PEI-DN08_S3"),1])/sqrt(3),
                        sd(dlooppeakD[c("FFM-19G1","FFM-20B2","FFM-22F10","MELC-A11_S1","NYTC-C9_S2"),1])/sqrt(5))



pdf("dloop_CN_plot.pdf", width=2, height=2) #,width=2, height=4)
    ggplot(dlooppeak_summary, aes(x=region,y=dlooppeak_norm,fill=region)) +
        geom_bar(position="dodge", stat="identity")+
        scale_fill_manual(values=c("black","red", "blue"))+ # 
        geom_point(data=dlooppeakD, aes(x=region,y=dlooppeak_norm),, color="grey", size=1)+ # position=position_dodge(.9)
        geom_errorbar(aes(ymin=(dlooppeak_norm-2*sem), ymax=(dlooppeak_norm+2*sem)), width=.2, color="grey")+ # grey46 ,position=position_dodge(.9)
        #xlab(NULL)+
        ylab("Estimated dloop copies")+
        scale_y_continuous(breaks = seq(0, 10, by = 1))+
        theme_classic() +
        theme(legend.position = "none",
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text=element_text(size=8,face="bold"),
            axis.title=element_text(size=10,face="bold"),
            text=element_text(size=10,face="bold"))#+
        #ggtitle("Estimated dloop copies")
dev.off()
