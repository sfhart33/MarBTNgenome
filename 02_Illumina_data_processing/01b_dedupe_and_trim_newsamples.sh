# New samples - dedupe and trim in one go

SAMPLES="DN-HL07 FFM-22A10 FFM-22F10"
DATA=/ssd2/Illumina_data/30-503541615/00_fastq
OUTPUT=/ssd2/Illumina_data/dedupe/2021_newsamples_trim20

cd $OUTPUT
for sample in $SAMPLES
do
/ssd2/Illumina_data/bbmap/clumpify.sh in1=$DATA/$sample"-hemocytes_R1_001.fastq.gz" in2=$DATA/$sample"-hemocytes_R2_001.fastq.gz" out1=$sample"_R1_001.dedupe.fastq.gz" out2=$sample"_R2_001.dedupe.fastq.gz" dedupe optical dupedist=2500 &
done 
wait
for sample in $SAMPLES
do
java -jar /ssd2/Illumina_data/dedupe/trim5/trimmomatic-0.36.jar PE -threads 10 -phred33 $sample"_R1_001.dedupe.fastq.gz" $sample"_R2_001.dedupe.fastq.gz" $sample"_R1_001.fastq.gz" $sample"_R1_001_unpaired.fastq.gz" $sample"_R2_001.fastq.gz" $sample"_R2_001_unpaired.fastq.gz" LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:25 &
done 



################### tissues

SAMPLES="FFM-22A10-adductor FFM-22F10-adductor MELC-A10-siphon MELC-A11-siphon NYTC-C9-mantle PEI-DN03-siphon PEI-DN07-siphon PEI-DN08-siphon"
DATA=/ssd2/Illumina_data/30-503543126/00_fastq
OUTPUT=/ssd2/Illumina_data/dedupe/2021_newsamples_trim20

cd $OUTPUT
for sample in $SAMPLES
do
/ssd2/Illumina_data/bbmap/clumpify.sh in1=$DATA/$sample"_R1_001.fastq.gz" in2=$DATA/$sample"_R2_001.fastq.gz" out1=$sample"_R1_001.dedupe.fastq.gz" out2=$sample"_R2_001.dedupe.fastq.gz" dedupe optical dupedist=2500 # &
done 
wait
for sample in $SAMPLES
do
java -jar /ssd2/Illumina_data/dedupe/trim5/trimmomatic-0.36.jar PE -threads 5 -phred33 $sample"_R1_001.dedupe.fastq.gz" $sample"_R2_001.dedupe.fastq.gz" $sample"_R1_001.fastq.gz" $sample"_R1_001_unpaired.fastq.gz" $sample"_R2_001.fastq.gz" $sample"_R2_001_unpaired.fastq.gz" LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:25 # &
done 

############### some failed for unknown reason: and needed to be rerun. Perhaps running in parallel caused an issue: in future don't use &

