# original file here: /ssd3/Mar_genome_analysis/delly/final_run/README

module load delly/


GENOME=/ssd3/Mar_genome_analysis/genomes/Mar.3.4.6.p1_Q30Q30A.fasta
INPUT=/ssd3/Mar_genome_analysis/bwa_mapping/Mar.3.4.6.p1/all_samples
OUTPUT=/ssd3/Mar_genome_analysis/delly/final_run

cd $INPUT
SAMPLES=$(ls *.bam | head -14)
echo $SAMPLES
cd $OUTPUT


for sample in $SAMPLES
do
	delly call -o $OUTPUT/$sample.delly.bcf -q 30 -g $GENOME $INPUT/$sample &
done
wait

#################################################

INPUT=/ssd3/Mar_genome_analysis/delly/final_run
OUTPUT=/ssd3/Mar_genome_analysis/delly/final_run
cd $INPUT
SAMPLES=$(ls *.bcf)
for sample in $SAMPLES
do
	echo "Starting " $sample
	bcftools view $INPUT/$sample | \
	awk ' $0 ~ /^#.*/ && $0 !~ /^##.*/ { split($10, name, "." ); print "CHR", $2, "TYPE", $4, $5, "END", "INS_LEN", name[2] "_NR", name[2] "_NV"};
	$0 !~ /^#.*/ {
		split($8, info, ";" );
		split(info[2], type, "=" );
		split(info[4], end, "=" );
		split(info[5], inslen, "=" );
		split($10, counts, ":" );
		total = counts[11] + counts[12];
		if($7 ~ "PASS" && info[1] ~ /^PRECISE.*/ && type[2] ~ "INS"){ print $1, $2, type[2], $4, $5, end[2], inslen[2], total, counts[12]}
		if($7 ~ "PASS" && info[1] ~ /^PRECISE.*/ && type[2] !~ "INS"){ print $1, $2, type[2], $4, $5, end[2], 0, total, counts[12]}

} ' OFS='\t' > $OUTPUT/$sample.delly.table &
done
cd $OUTPUT

