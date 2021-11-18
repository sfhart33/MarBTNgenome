# Find location of Steamer in new haploid genomes

GENOMES=/ssd3/Mar_genome_analysis/genomes
# two haploid genomes to compare
LIST="Mar.3.3.4.p0 Mar.3.3.4.p1"
WORK=/ssd3/Mar_genome_analysis/steamer

################ make bwa indexes for genomes
cd $WORK
for i in $LIST
do
	bwa index $GENOMES/$i.fasta
done

################ map Steamer to genomes
for i in $LIST
do
	bwa mem -t 5 $GENOMES/$i.fasta enS6a1-Steamer.fasta > $i.enS6a1-Steamer.sam
	bwa mem -t 5 -a $GENOMES/$i.fasta enS6a1-Steamer.fasta > $i.enS6a1-SteamerALL.sam
done


### Find location in final genome: Mar.3.4.6.p1_Q30Q30A
GENOME=/ssd3/Mar_genome_analysis/genomes/Mar.3.4.6.p1_Q30Q30A.fasta
bwa index $GENOME
STEAMER=/ssd2/Steamer_pipeline/references/Steamer_full.fasta
LTR=/ssd2/Steamer_pipeline/references/SteamerLTRonly.fasta
WORK=/ssd3/Mar_genome_analysis/steamer
cd $WORK
bwa mem -t 5 -a $GENOME $STEAMER > Mar.3.4.6.p1_Q30Q30A.Steamer_full.fasta.sam
bwa mem -t 5 -a $GENOME $LTR > Mar.3.4.6.p1_Q30Q30A.SteamerLTRonly.fasta.sam

