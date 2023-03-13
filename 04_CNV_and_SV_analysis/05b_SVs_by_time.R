# Full notes here: C:\Users\shart\Metzger Lab Dropbox\Sam_data\Mya_genome\revision\genome_stability_by_time.r

setwd("/ssd3/Mar_genome_analysis/RNAseq/paper_revisions")
    # combined manually into single file:
    all_data <- read.table("stability_by_date.txt", header=T)
        #        sample region  BND  DUP     dloop steamer   date
        # 1   PEI.DF488    PEI  747  725  5.858938     290 -10.83
        # 2    PEI.DN03    PEI  796  755  6.264484     291 -11.28
        # 3 PEI.DN08_S3    PEI  633  676  5.249071     275 -11.28
        # 4   FFM.22F10    USA 1146 1328 11.051997     438  -0.40
        # 5    FFM.19G1    USA 1273 1416  9.351337     460  -0.98
        # 6    FFM.20B2    USA 1225 1423  8.587608     455  -0.98
        # 7 MELC.A11_S1    USA 1132 1357 10.030630     437  -7.94
        # 8  NYTC.C9_S2    USA 1091 1392 10.134541     424  -7.37

    usa_data <- filter(all_data, region == "USA")

    pdf("stability_over_time.pdf",width=3,height=3)
        for(i in c("BND", "DUP","dloop","steamer")){
        plot1 <- ggplot(usa_data, aes(x = jitter(date+2021.73), y = get(i))) + 
                scale_color_manual(values=c("red","blue"))+
                stat_smooth(method = "lm",color="blue",fullrange = T) + 
                geom_point(data = all_data, aes(x = jitter(date+2021.73), y =  get(i), color=region))+
                scale_x_continuous(breaks=seq(2010,2022,2))+
                xlab("Year sampled")+
                ylab(i)+   
                scale_fill_manual(values=c("red","blue"))+
                theme_classic() +
                theme(axis.text=element_text(size=8,face="bold"),
                    axis.title=element_text(size=8,face="bold"),
                    text=element_text(size=8,face="bold"),
                    title=element_text(size=8, face='bold'),
                    plot.title = element_text(size=8, face='bold'),
                    legend.position = "none")
            print(plot1) 
        }
    dev.off()

