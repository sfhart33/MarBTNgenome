# original file here: C:\Users\shart\Metzger Lab Dropbox\Sam_data\Mya_genome\downsampling_illumina.sh

# downsample all to ~600 million reads before running delly
    
    # TEST_BAM=/ssd3/Mar_genome_analysis/bwa_mapping/Mar.3.4.6.p1/all_samples/01.MELC-2E11.bam
    # COUNT=$(samtools view -@ 50 -c $TEST_BAM)
    # FRACTION=$(awk -v COUNT=$COUNT "BEGIN {print 600000000 / COUNT}") # use awk to calculate fraction
    # echo $FRACTION
    # echo $COUNT

    # samtools view -@ 50 --subsample $FRACTION --subsample-seed 1234 $TEST_BAM


module load delly/
module load samtools/1.9


GENOME=/ssd3/Mar_genome_analysis/genomes/Mar.3.4.6.p1_Q30Q30A.fasta
INPUT=/ssd3/Mar_genome_analysis/bwa_mapping/Mar.3.4.6.p1/all_samples
OUTPUT=/ssd3/Mar_genome_analysis/delly/downsampled

cd $INPUT
SAMPLES=$(ls *.bam | head -14)
echo $SAMPLES
cd $OUTPUT

# sub-sample and generate new bam file
for sample in $SAMPLES
do
	echo $sample 
    COUNT=$(samtools view -@ 10 -c $INPUT/$sample)
    echo $COUNT
    FRACTION=$(awk -v COUNT=$COUNT "BEGIN {print 600000000 / COUNT}") # use awk to calculate fraction
    echo $FRACTION
    samtools view -b -h -@ 5 -s $FRACTION $INPUT/$sample | samtools sort -O bam -@ 5 > $OUTPUT/$sample &
done
wait

# index bam files
echo "Starting Index"
for sample in $SAMPLES
do
	samtools index -b -@ 10 $OUTPUT/$sample &
done
wait 

# run delly
echo "Starting delly"
for sample in $SAMPLES
do
	delly call -o $OUTPUT/$sample.delly.bcf -q 30 -g $GENOME $OUTPUT/$sample &
done
wait 

################### Ran above

# Check on the read counts
for sample in $SAMPLES
do
	echo $sample >> $OUTPUT/final.counts
    samtools view -@ 5 -c $OUTPUT/$sample >> $OUTPUT/final.counts
done
wait
# 600003992
# 599997749
# 600014207
# 599992801
# 599990107
# 599962640
# 599977262
# 600000378
# 600014316
# 599998507
# 600026747
# 600014103
# 600000573

# delete bam and index files after 
for sample in $SAMPLES
do
	rm $OUTPUT/$sample
	rm $OUTPUT/$sample.bai
done
wait 

# re-run downstream delly processes

INPUT=/ssd3/Mar_genome_analysis/delly/downsampled
OUTPUT=/ssd3/Mar_genome_analysis/delly/downsampled
cd $INPUT
SAMPLES=$(ls *.bcf)
echo $SAMPLES

for sample in $SAMPLES
do
	echo "Starting " $sample
	bcftools view $INPUT/$sample | \
	awk ' $0 ~ /^#.*/ && $0 !~ /^##.*/ { split($10, name, "." ); print "CHR", $2, "TYPE", $4, $5, "END", "INS_LEN", name[2] "_NR", name[2] "_NV"};
    $0 !~ /^#.*/ { split($8, info, ";" );
        split(info[2], type, "=" );
        split(info[4], end, "=" );
        split(info[5], inslen, "=" );
        split($10, counts, ":" );
        total = counts[11] + counts[12];
        if($7 ~ "PASS" && info[1] ~ /^PRECISE.*/ && type[2] ~ "INS"){print $1, $2, type[2], $4, $5, end[2], inslen[2], total, counts[12]};
        if($7 ~ "PASS" && info[1] ~ /^PRECISE.*/ && type[2] !~ "INS"){ print $1, $2, type[2], $4, $5, end[2], 0, total, counts[12]}
    } ' OFS='\t' > $OUTPUT/$sample.delly.table
done

cd $OUTPUT


