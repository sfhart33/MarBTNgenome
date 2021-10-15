# Running dedupe and clumpify programs to get rid of PCR duplicates from the data

# Actual run ####################################################
SAMPLES="MELC-2E11 MELC-A11_S1 PEI-DN08_S3 PEI-DF490 MELC-A9 NYTC-C9_S2 MELC-A10"
FULLDATA=/ssd2/Illumina_data
cd /ssd2/Illumina_data/dedupe
for sample in $SAMPLES
do
/ssd2/Illumina_data/bbmap/clumpify.sh in1=$FULLDATA/$sample"_R1_001.fastq.gz" in2=$FULLDATA/$sample"_R2_001.fastq.gz" out1=$sample"_R1_001.fastq.gz" out2=$sample"_R2_001.fastq.gz" dedupe optical dupedist=2500
done
# origionally ran in /ssd2/Steamer_pipeline/dedupe
#################################################################
# run for new samples
SAMPLES="DF-488 DN-HL03 FFM-19G1 FFM-20B2"
FULLDATA=/ssd2/Illumina_data
cd /ssd2/Illumina_data/dedupe
for sample in $SAMPLES
do
/ssd2/Illumina_data/bbmap/clumpify.sh in1=$FULLDATA/$sample"_R1_001.fastq.gz" in2=$FULLDATA/$sample"_R2_001.fastq.gz" out1=$sample"_R1_001.fastq.gz" out2=$sample"_R2_001.fastq.gz" dedupe optical dupedist=2500
done

#at end, deleted original illumina files to save space

# trim 20
#SAMPLES="MELC-2E11 MELC-A11_S1 PEI-DN08_S3 PEI-DF490 MELC-A9 NYTC-C9_S2 MELC-A10"
SAMPLES="DF-488 DN-HL03 FFM-19G1 FFM-20B2"
FULLDATA=/ssd2/Illumina_data/dedupe
cd /ssd2/Illumina_data/dedupe/trim20
for sample in $SAMPLES
do
java -jar /ssd2/Illumina_data/dedupe/trim5/trimmomatic-0.36.jar PE -threads 20 -phred33 $FULLDATA/$sample"_R1_001.fastq.gz" $FULLDATA/$sample"_R2_001.fastq.gz" $sample"_R1_001.fastq.gz" $sample"_R1_001_unpaired.fastq.gz" $sample"_R2_001.fastq.gz" $sample"_R2_001_unpaired.fastq.gz" LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:25
done
