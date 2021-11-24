# Calculate trinucleotide freq in genome, exome, CDS, etc to correct for mutational oppertunities

# Make fasta files for gene regions
module load bedtools
GENOME=/ssd3/Mar_genome_analysis/genomes/Mar.3.4.6.p1_Q30Q30A.fasta
REGIONS=/ssd3/Mar_genome_analysis/genomes/maker
OUTPUT=/ssd3/Mar_genome_analysis/mut_sig/trinuc_freq
TYPES="gene CDS exon five_prime_UTR three_prime_UTR"

for type in $TYPES
do
bedtools getfasta -fi $GENOME -bed $REGIONS/2020-09-11_Mar_genome_snap02.all.$type.bed |\
awk 'BEGIN{print ">sequence"} $0 !~ /^>.*/ {print $1}' > $OUTPUT/Mar.3.4.6.p1.$type.fasta
done

awk 'BEGIN{print ">sequence"} $0 !~ /^>.*/ {print $1}' $GENOME > $OUTPUT/Mar.3.4.6.p1.genome.fasta
