######## Mapping trim20 reads to genome ################

GENOME=/ssd3/Mar_genome_analysis/genomes/Mar.3.4.6.p1_Q30Q30A.fasta
LIST="MELC-A11_S1 PEI-DN08_S3 NYTC-C9_S2 MELC-2E11 MELC-A10 MELC-A9 PEI-DF490 DF-488 DN-HL03 FFM-19G1 FFM-20B2"
WORK=/ssd3/Mar_genome_analysis/bwa_mapping/Mar.3.4.6.p1/all_samples
ILLUMINA=/ssd2/Illumina_data/dedupe/trim20

cd $WORK
for i in $LIST
do
	bwa mem -t 50 $GENOME $ILLUMINA/$i"_R1_001.fastq.gz" $ILLUMINA/$i"_R2_001.fastq.gz" | samtools view -b -h -@ 10 | samtools sort -O bam -@ 10 > Mar.3.4.6.p1.$i.bam
done

for i in $LIST
do
	samtools index -b -@ 20 Mar.3.4.6.p1.$i.bam
done

############## Map new samples seperately

LIST="DN-HL07 FFM-22A10 FFM-22F10 FFM-22A10-adductor FFM-22F10-adductor MELC-A10-siphon MELC-A11-siphon NYTC-C9-mantle PEI-DN03-siphon PEI-DN07-siphon PEI-DN08-siphon"
# LIST1="FFM-22A10-adductor FFM-22F10-adductor MELC-A10-siphon MELC-A11-siphon NYTC-C9-mantle PEI-DN03-siphon PEI-DN07-siphon PEI-DN08-siphon"
# LIST2="DN-HL07 FFM-22A10 FFM-22F10"
GENOME=/ssd3/Mar_genome_analysis/genomes/Mar.3.4.6.p1_Q30Q30A.fasta
WORK=/ssd3/Mar_genome_analysis/bwa_mapping/Mar.3.4.6.p1/all_samples
ILLUMINA=/ssd2/Illumina_data/dedupe/2021_newsamples_trim20

cd $WORK
for i in $LIST
do
	bwa mem -t 80 -v 1 $GENOME $ILLUMINA/$i"_R1_001.fastq.gz" $ILLUMINA/$i"_R2_001.fastq.gz" | samtools view -b -h -@ 10 | samtools sort -O bam -@ 10 > Mar.3.4.6.p1.$i.bam
done

for i in $LIST
do
	samtools index -b -@ 10 Mar.3.4.6.p1.$i.bam &
done
wait 

