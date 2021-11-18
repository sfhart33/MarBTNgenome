# Merge libaries from repeatmodeler (reference genome), repdenovo (cancer)
# Merge repdenovo (healthy too - to control for bias?)

RM_LIB=/ssd3/Mar_genome_analysis/repeatmodeler/Mar.3.4.6.p1_Q30Q30A-families.fa
RD_LIB_USA=/ssd2/REPdenovo/MELC_A11_output/contigs.fa.classified
RD_LIB_PEI=/ssd2/REPdenovo/PEI-DN08_output/contigs.fa.classified
RD_LIB_REF=/ssd2/REPdenovo/MELC-2E11_output/contigs.fa.classified
COMMAND=/ssd3/Mar_genome_analysis/repeat_lib/cd-hit-v4.8.1-2019-0228
WORKD=/ssd3/Mar_genome_analysis/repeat_lib

# # Install CD-HIT - program used by repeatmodeler2 to merge libraries
# cd /ssd3/Mar_genome_analysis/repeat_lib
# wget https://github.com/weizhongli/cdhit/releases/download/V4.8.1/cd-hit-v4.8.1-2019-0228.tar.gz
# tar xvf cd-hit-v4.8.1-2019-0228.tar.gz --gunzip
# cd cd-hit-v4.8.1-2019-0228
# make
# cd cd-hit-auxtools
# make

##  Run algorithm
# $COMMAND/cd-hit-est-2d -i $RD_LIB_USA -i2 $RD_LIB_PEI -o $WORKD/repdenovo_usa-pei_merge
# cat $RD_LIB_USA $RD_LIB_PEI $RD_LIB_REF > $WORKD/repdenovo_UPR_cat
# $COMMAND/cd-hit-est -aS 0.8 -c 0.8 -g 1 -G 0 -A 80 -M 10000 -i $WORKD/repdenovo_UPR_cat -o $WORKD/repdenovo_UPR_combined -T 20


############################ 4/17/21
#copy over 
cd $WORKD
cp $RM_LIB ref.genome.repeatmoduler
cp $RD_LIB_USA usa.repdenovo
cp $RD_LIB_PEI pei.repdenovo
cp $RD_LIB_REF ref.repdenovo
awk '$0 !~ /^>.*/ {print} $0 ~ /^>.*/ {split($1, name, "#" ); print name[1] "#usaRD#" name[2]}' usa.repdenovo > usa.repdenovo.named
awk '$0 !~ /^>.*/ {print} $0 ~ /^>.*/ {split($1, name, "#" ); print name[1] "#peiRD#" name[2]}' pei.repdenovo > pei.repdenovo.named
awk '$0 !~ /^>.*/ {print} $0 ~ /^>.*/ {split($1, name, "#" ); print name[1] "#refRD#" name[2]}' ref.repdenovo > ref.repdenovo.named
awk '$0 !~ /^>.*/ {print} $0 ~ /^>.*/ {split($1, name, "#" ); print name[1] "#refRM#" name[2]}' ref.genome.repeatmoduler > ref.genome.repeatmoduler.named
cat usa.repdenovo.named pei.repdenovo.named ref.repdenovo.named ref.genome.repeatmoduler.named > all_repeat_lib_cat
$COMMAND/cd-hit-est -aS 0.8 -c 0.8 -g 1 -G 0 -A 80 -M 10000 -i $WORKD/all_repeat_lib_cat -o $WORKD/all_repeat_lib.cdhit -T 50
FILES="usa.repdenovo.named pei.repdenovo.named ref.repdenovo.named ref.genome.repeatmoduler.named all_repeat_lib.cdhit"
for i in $FILES
do
awk '$0 ~ /^>.*/ {count += 1} END{print count}' $i
done

# bwa map reads to library, keep only mapping reads
LIST="MELC-A11_S1 PEI-DN08_S3 NYTC-C9_S2 MELC-2E11 MELC-A9 PEI-DF490 DF-488 DN-HL03 FFM-19G1 FFM-20B2"
ILLUMINA=/ssd2/Illumina_data/dedupe/trim20
bwa index all_repeat_lib.cdhit

for i in $LIST
do
	bwa mem -t 10 -v 2 all_repeat_lib.cdhit $ILLUMINA/$i"_R1_001.fastq.gz" $ILLUMINA/$i"_R2_001.fastq.gz" | samtools view -b -h -F 4 | samtools sort -O bam > $WORKD/repeatLIB.$i.bam &
done
wait
for i in $LIST
do
	samtools index -b -@ 10 $WORKD/repeatLIB.$i.bam &
done
wait

# get total nucleotides mapped to each 
#samtools depth (would then need to add up for every repeat)
#samtools coverage (computes summary stats on its own)
cd $WORKD
for i in $LIST
do
	# samtools coverage repeatLIB.$i.bam -o repeatLIB.$i.coverage &
	samtools depth -a -a -d 0 repeatLIB.$i.bam > repeatLIB.$i.coverage &
done
wait

# Use awk of R to count up total coverage of each 
for i in $LIST
do
awk -v sample=$i -v FS='\t' -v OFS='\t' 'BEGIN{print "repeat_name", "length", sample;
	name = "NEW_CONTIG_MERGE_2#usaRD#Unknown";
	if(sample ~ "FFM-19G1" ){ avgcov = 96.8552};
	if(sample ~ "FFM-20B2" ){ avgcov = 96.8232};
	if(sample ~ "MELC-A11_S1" ){ avgcov = 65.6252};
	if(sample ~ "NYTC-C9_S2" ){ avgcov = 61.4114};
	if(sample ~ "DF-488" ){ avgcov = 98.3213};
	if(sample ~ "DN-HL03" ){ avgcov = 92.4672};
	if(sample ~ "PEI-DN08_S3" ){ avgcov = 66.9775};
	if(sample ~ "MELC-2E11" ){ avgcov = 58.7932};
	if(sample ~ "MELC-A9" ){ avgcov = 60.0382};
	if(sample ~ "PEI-DF490" ){ avgcov = 60.0636}} 
	name ~ $1 {rep_length += 1; coverage += $3};
	name !~ $1 {print name, rep_length, coverage/rep_length/avgcov;
		name = $1; rep_length = 1; coverage = $3}' OFS='\t' repeatLIB.$i.coverage > repeatLIB.$i.depth &
done
