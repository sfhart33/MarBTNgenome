# Calculate trinucleotide freq in genome, exome, CDS, etc to correct for mutational oppertunities

library(Biostrings)
library(sigfit)
data("cosmic_signatures_v2")

setwd("/ssd3/Mar_genome_analysis/mut_sig/trinuc_freq")
genome <- readDNAStringSet("Mar.3.4.6.p1.genome.fasta")
cds <- readDNAStringSet("Mar.3.4.6.p1.CDS.fasta")
exon <- readDNAStringSet("Mar.3.4.6.p1.exon.fasta")
gene <- readDNAStringSet("Mar.3.4.6.p1.gene.fasta")
five_prime <- readDNAStringSet("Mar.3.4.6.p1.five_prime_UTR.fasta")
three_prime <- readDNAStringSet("Mar.3.4.6.p1.three_prime_UTR.fasta")

genome2 <- trinucleotideFrequency(genome, as.prob=TRUE)
cds2 <- trinucleotideFrequency(cds, as.prob=TRUE)
exon2 <- trinucleotideFrequency(exon, as.prob=TRUE)
gene2 <- trinucleotideFrequency(gene, as.prob=TRUE)
five_prime2 <- trinucleotideFrequency(five_prime, as.prob=TRUE)
three_prime2 <- trinucleotideFrequency(three_prime, as.prob=TRUE)
mut_prob <- rbind(genome2, gene2, exon2, cds2, five_prime2, three_prime2)
rownames(mut_prob) <- c("genome", "gene", "exon", "cds","fiveprimeUTR","threeprimeUTR")
mut <- as.data.frame(mut_prob)
mut_C <- mutate(mut,
	"ACA>" = (ACA + TGT)/3, 
	"ACC>" = (ACC + GGT)/3, 
	"ACG>" = (ACG + CGT)/3, 
	"ACT>" = (ACT + AGT)/3,
	"CCA>" = (CCA + TGG)/3, 
	"CCC>" = (CCC + GGG)/3, 
	"CCG>" = (CCG + CGG)/3, 
	"CCT>" = (CCT + AGG)/3,
	"GCA>" = (GCA + TGC)/3, 
	"GCC>" = (GCC + GGC)/3, 
	"GCG>" = (GCG + CGC)/3, 
	"GCT>" = (GCT + AGC)/3,
	"TCA>" = (TCA + TGA)/3, 
	"TCC>" = (TCC + GGA)/3, 
	"TCG>" = (TCG + CGA)/3, 
	"TCT>" = (TCT + AGA)/3)
mut_T <- mutate(mut,
	"ATA>" = (ATA + TAT)/3, 
	"ATC>" = (ATC + GAT)/3, 
	"ATG>" = (ATG + CAT)/3, 
	"ATT>" = (ATT + AAT)/3,
	"CTA>" = (CTA + TAG)/3, 
	"CTC>" = (CTC + GAG)/3, 
	"CTG>" = (CTG + CAG)/3, 
	"CTT>" = (CTT + AAG)/3,
	"GTA>" = (GTA + TAC)/3, 
	"GTC>" = (GTC + GAC)/3, 
	"GTG>" = (GTG + CAC)/3, 
	"GTT>" = (GTT + AAC)/3,
	"TTA>" = (TTA + TAA)/3, 
	"TTC>" = (TTC + GAA)/3, 
	"TTG>" = (TTG + CAA)/3, 
	"TTT>" = (TTT + AAA)/3) 
mut_opp <- cbind(mut_C[,65:80],mut_C[,65:80],mut_C[,65:80],mut_T[,65:80],mut_T[,65:80],mut_T[,65:80])

colnames(mut_opp) <- colnames(cosmic_signatures_v2)
rownames(mut_opp) <-  c("genome", "gene", "exon", "cds","fiveprimeUTR","threeprimeUTR")
mut_opp
rowSums(mut_opp) # all equal 1
plot_spectrum(mut_opp, pdf_path = "mutational_oppertunities.pdf")

mut_opp <- mut_opp[order(row.names(mut_opp)), ]
setwd("/ssd3/Mar_genome_analysis/mut_sig/trinuc_freq")
saveRDS(mut_opp, file = "mutational_oppertunities_all.rds")
saveRDS(mut_opp[5,], file = "mutational_oppertunities_genome.rds")
write.table(mut_opp, file = "mutational_oppertunities_all.txt")