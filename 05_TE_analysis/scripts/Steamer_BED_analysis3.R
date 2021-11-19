# install.packages("tidyverse")
library(tidyverse)

#upload file, name columns, get rid of unneccessery columns, filter qual
intreads <- read.delim("intreads.bed", header = FALSE)
intreads <- intreads[,1:9]
colnames(intreads) <- c("contig", "start", "end", "read", "qual", "strand", "flank", "pos", "name")
# intreads$dunno <- NULL
intreads$read <- NULL
intreads$pos <- NULL
intreads2 <- arrange(intreads, contig, start)
intreads2$qual <- NULL

# Combine duplicates and count how many reads for each of up and down
grouped_reads <- group_by(intreads2, name) %>%
  mutate(count = n())
grouped_reads2 <- group_by(intreads2, name, flank) %>%
  mutate(count = n()) %>%
  distinct()
up_reads <- filter(grouped_reads2, flank == "up") %>%
  rename(up = count)
dn_reads <- filter(grouped_reads2, flank == "down") %>%
  rename(down = count)
up_reads$flank <- NULL
dn_reads$flank <- NULL
final_reads <- full_join(up_reads, dn_reads, by = c("contig", "start", "end", "strand", "name")) %>%
  arrange(contig, start)
final_reads[is.na(final_reads)] <- 0
# final_reads <- as.data.frame(final_reads) 
final_reads$total <- (final_reads$up + final_reads$down)
#print(final_reads)
both_reads <- filter(final_reads, up > 0 & down > 0) # %>% print()
many_reads <- filter(final_reads, up + down > 1)     # %>% print()
#ten_reads <- filter(final_reads, total > 9)          # %>% print()
write.table(many_reads, file = "insertion_sites_multiple.bed", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(both_reads, file = "insertion_sites_updown.bed", sep = "\t", quote = FALSE, row.names = FALSE)
