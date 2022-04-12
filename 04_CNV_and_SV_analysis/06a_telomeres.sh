# From here: /ssd3/Mar_genome_analysis/telomeres/README

module load telseq
BAMFILES=/ssd3/Mar_genome_analysis/bwa_mapping/Mar.3.4.6.p1/all_samples
TELOMERES=/ssd3/Mar_genome_analysis/telomeres
OUTPUT=/ssd3/Mar_genome_analysis/telomeres/outputs
cd $BAMFILES

# TEST RUN
	# telseq 15.PEI-DN03-siphon.bam 21.MELC-A11-siphon.bam > $TELOMERES/telseq_test_output

# RUN ALL IN ONE RUN
	# ls *.bam | telseq > $TELOMERES/telseq_full_output

# RUN ALL SEPERATELY IN PARELLEL - DIFFERENT OUTPUT FILES
FILES=$(ls *.bam)
for i in $FILES
do
telseq $i > $OUTPUT/$i.telseq &
done
