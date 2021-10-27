
library(tidyverse)
setwd("/ssd3/Mar_genome_analysis/bwa_mapping/mito/all_samples/depth")

depth_files <- read.delim("MELC-A10.mtgenome.depth", head=FALSE) %>% select(V2)%>% rename(pos=V2)

for(sample in c("DF-488", "DN-HL03", "DN-HL07", "FFM-19G1", "FFM-20B2", "FFM-22A10", "FFM-22F10", "MELC-2E11", "MELC-A10", "MELC-A11_S1", "MELC-A9", "NYTC-C9_S2", "PEI-DF490", "PEI-DN08_S3")){
	depth_file <- read.delim(paste0(sample,".mtgenome.depth"), head=FALSE) %>% select(V3) %>% rename(!!sample:=V3)
	depth_files <- cbind(depth_files, depth_file)
}


head(depth_files)

colMeans(depth_files)

pdf("depth.pdf")
ggplot(depth_files) +
	geom_line(aes(pos, get("DF-488")), color = "red")+
	geom_line(aes(pos, get("DN-HL03")), color = "red")+
	geom_line(aes(pos, get("DN-HL07")), color = "red")+
	geom_line(aes(pos, get("PEI-DN08_S3")), color = "red")+
	geom_line(aes(pos, get("FFM-19G1")), color = "blue")+
	geom_line(aes(pos, get("FFM-20B2")), color = "blue")+
	geom_line(aes(pos, get("FFM-22A10")), color = "blue")+
	geom_line(aes(pos, get("FFM-22F10")), color = "blue")+
	geom_line(aes(pos, get("MELC-A10")), color = "blue")+
	geom_line(aes(pos, get("MELC-A11_S1")), color = "blue")+
	geom_line(aes(pos, get("NYTC-C9_S2")), color = "blue")+
	geom_line(aes(pos, get("MELC-2E11")), color = "green")+
	geom_line(aes(pos, get("MELC-A9")), color = "green")+
	geom_line(aes(pos, get("PEI-DF490")), color = "green")+
	ylab("read depth")+
        theme_classic() +
            theme(axis.text=element_text(size=12,face="bold"),
            axis.title=element_text(size=16,face="bold"),
            text=element_text(size=18,face="bold")) +
        ggtitle("Mt genome coverage (PEI=red,USA=blue,Healthy=green)")
for(sample in c("DF-488", "DN-HL03", "DN-HL07", "FFM-19G1", "FFM-20B2", "FFM-22A10", "FFM-22F10", "MELC-2E11", "MELC-A10", "MELC-A11_S1", "MELC-A9", "NYTC-C9_S2", "PEI-DF490", "PEI-DN08_S3")){
plot1 <- ggplot(depth_files) +
	geom_line(aes(pos, get(sample)))+
	ylab("read depth")+
        theme_classic() +
            theme(axis.text=element_text(size=12,face="bold"),
            axis.title=element_text(size=16,face="bold"),
            text=element_text(size=18,face="bold")) +
        ggtitle(paste("Mt genome coverage", sample))
print(plot1)
}
dev.off()


dloop <- depth_files %>% filter(pos > 12164 & pos < 12870)
notdloop <- depth_files %>% filter(pos < 12164 | pos > 12870)

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



snvs <- read.delim("/ssd3/Mar_genome_analysis/bwa_mapping/mito/all_samples/somatypus/Somatypus_SNVs_final.counts", header = TRUE)
head(snvs)

#not in dloop region
snvs_notdloop <- filter(snvs, pos < 12060 | pos > 12971)
snvs_dloop <- filter(snvs, pos > 12164 & pos < 12870)


pdf("dloop_snp_freq.pdf")
for(i in c(	"Href_f", "Husa_f", "Hpei_f", 
	"Cpei0_f", "Cpei1_f", "Cpei2_f", "Cpei3_f",
	"Cusa0a_f", "Cusa0b_f", "Cusa1_f", "Cusa2_f", "Cusa3_f", "Cusa4_f", "Cusa5_f")){ #	"Tpei1_f", "Tpei2_f", "Tpei3_f", "Tusa1_f", "Tusa2_f", "Tusa3_f", "Tusa4_f", "Tusa5_f"
	snvs_dloop1 <- filter(snvs_dloop, get(i) > 0.05)
	plot1 <- ggplot(snvs_dloop1, aes(get(i)))+
        geom_histogram(binwidth=0.01)+
        ggtitle(i)+
        theme_classic() +
        theme(axis.text=element_text(size=12,face="bold"),
                axis.title=element_text(size=16,face="bold"),
                text=element_text(size=18,face="bold"))
    print(plot1)
}
dev.off()


snvs_dloop %>% mutate(usa_avg = (Cusa0a_f + Cusa0b_f + Cusa1_f + Cusa2_f + Cusa3_f + Cusa4_f + Cusa5_f)/7, pei_avg = (Cpei0_f + Cpei1_f + Cpei2_f + Cpei3_f)/4)

write.table(snvs_dloop, "dloop_snps.txt", append = FALSE, row.names=FALSE, col.names=TRUE, quote=FALSE, sep = '\t')

depth_files2 <- depth_files
depth_files2[,2:15] <- depth_files[,2:15]/colMeans(depth_files[,2:15])
head(depth_files2)
depth_files2[,2:15] <- as.data.frame((t(t(as.matrix(depth_files[,2:15]))/colMeans(depth_files[,2:15], na.rm=TRUE))))

pdf("depth_norm.pdf")
ggplot(depth_files2) +
	geom_line(aes(pos, get("DF-488")), color = "red")+
	geom_line(aes(pos, get("DN-HL03")), color = "red")+
	geom_line(aes(pos, get("DN-HL07")), color = "red")+
	geom_line(aes(pos, get("PEI-DN08_S3")), color = "red")+
	geom_line(aes(pos, get("FFM-19G1")), color = "blue")+
	geom_line(aes(pos, get("FFM-20B2")), color = "blue")+
	geom_line(aes(pos, get("FFM-22A10")), color = "blue")+
	geom_line(aes(pos, get("FFM-22F10")), color = "blue")+
	geom_line(aes(pos, get("MELC-A10")), color = "blue")+
	geom_line(aes(pos, get("MELC-A11_S1")), color = "blue")+
	geom_line(aes(pos, get("NYTC-C9_S2")), color = "blue")+
	geom_line(aes(pos, get("MELC-2E11")), color = "green")+
	geom_line(aes(pos, get("MELC-A9")), color = "green")+
	geom_line(aes(pos, get("PEI-DF490")), color = "green")+
	ylab("Mt read depth normalized to avg depth")+
        theme_classic() +
            theme(axis.text=element_text(size=12,face="bold"),
            axis.title=element_text(size=16,face="bold"),
            text=element_text(size=18,face="bold")) +
        ggtitle("PEI=red,USA=blue,Healthy=green")
for(sample in c("DF-488", "DN-HL03", "DN-HL07", "FFM-19G1", "FFM-20B2", "FFM-22A10", "FFM-22F10", "MELC-2E11", "MELC-A10", "MELC-A11_S1", "MELC-A9", "NYTC-C9_S2", "PEI-DF490", "PEI-DN08_S3")){
plot1 <- ggplot(depth_files2) +
	geom_line(aes(pos, get(sample)))+
	ylab("read depth")+
        theme_classic() +
            theme(axis.text=element_text(size=12,face="bold"),
            axis.title=element_text(size=16,face="bold"),
            text=element_text(size=18,face="bold")) +
        ggtitle(paste("Mt genome coverage", sample))
print(plot1)
}
dev.off()

