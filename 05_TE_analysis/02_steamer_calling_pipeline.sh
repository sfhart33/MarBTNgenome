#!/bin/bash

# Original file here: C:\Users\shart\Metzger Lab Dropbox\Sam_data\Mya_genome\steamer\final_pipeline_run.sh

REF_DIR=/ssd2/Steamer_pipeline/references
WORK_DIR=/ssd3/Mar_genome_analysis/steamer/final_pipeline
FASTQ_DIR=/ssd2/Illumina_data/dedupe/trim20
SCRIPTS_DIR=/ssd2/Steamer_pipeline
GENOME=/ssd3/Mar_genome_analysis/genomes/Mar.3.4.6.p1_Q30Q30A.fasta
REPEAT_LIB=/ssd3/Mar_genome_analysis/repeatmodeler/Mar.3.4.6.p1_Q30Q30A-families.fa
#REPEAT_LIB=/ssd3/Mar_genome_analysis/repeat_lib/all_repeat_lib.cdhit
SAMPLES="MELC-2E11 MELC-A9 PEI-DF490 FFM-19G1 FFM-20B2 MELC-A11_S1 NYTC-C9_S2 DF-488 DN-HL03 PEI-DN08_S3 FFM-22F10" # MELC-A10 
#SAMPLES="NYTC-C9_S2 DF-488 DN-HL03 PEI-DN08_S3 FFM-22F10" # SUBSET


#Set quality threshold and other variables here so they can be easily changed
#MIN_MAPQ=30
MIN_QUAL=30

######################################################
echo "Mapping raw reads to steamer for discordant reads and LTR to extract split reads"
######################################################

for read in $SAMPLES
do
    cd $WORK_DIR
    mkdir $read
    cd $read
    if [ $read = "FFM-22F10" ]; then
        FASTQ_DIR=/ssd2/Illumina_data/dedupe/2021_newsamples_trim20
    fi
# only reads that dont map, but with mapped mate: for discordant reads
    echo "Mapping" $read "discordant reads"
    bwa mem -t 70 -v 1 $REF_DIR/Steamer_full.fasta \
        $FASTQ_DIR/$read"_R1_001.fastq.gz" $FASTQ_DIR/$read"_R2_001.fastq.gz" | \
        samtools fastq -f 4 -F 8 -@ 10 > discordant.fastq # without -n option to keep /1 or /2 at the end
    bwa mem -v 1 $GENOME discordant.fastq | samtools view -F 4 -q $MIN_QUAL > discordant_genome_mapped.sam
# only reads that map: to extract split reads
    echo "Mapping" $read "split reads"
    bwa mem -t 70 -v 1 -SP $REF_DIR/SteamerLTRonly.fasta \
        $FASTQ_DIR/$read"_R1_001.fastq.gz" $FASTQ_DIR/$read"_R2_001.fastq.gz" | \
        samtools view -F 4 -h -@ 10 > Steamer_LTRonly.sam 
    echo $read "complete"
done
echo "Mapping complete"

# Load R
module load R/3.6.0
# Package tidyverse must be installed in R

for sample in $SAMPLES
do
    cd $WORK_DIR/$sample
    echo $sample": extract flanking reads and map to genome"
# extract soft clipped
    perl $SCRIPTS_DIR/extractsoftclipseq2_SH.pl Steamer_LTRonly.sam
# Map flanking reads to Steamer and only keep reads that DONT map
    bwa mem $REF_DIR/Steamer_full.fasta softclips_down.fastq > dn_steamer_mapped.sam
    samtools fastq -f 4 dn_steamer_mapped.sam > softclips_down_flanking.fastq
    rm dn_steamer_mapped.sam 
    bwa mem $REF_DIR/Steamer_full.fasta softclips_up.fastq > up_steamer_mapped.sam
    samtools fastq -f 4 up_steamer_mapped.sam > softclips_up_flanking.fastq
    rm up_steamer_mapped.sam 
# Map flanking reads to full genome
    bwa mem $GENOME softclips_down_flanking.fastq > down_genome_mapped.sam
    bwa mem $GENOME softclips_up_flanking.fastq > up_genome_mapped.sam
# Filter out only reads with MAPQ < 30 or unmapped, and longer than 10bp
    awk '{if ($1 ~ /^@.*/) {print $0} 
        if ($1 !~ /^@.*/ && ($2 == 4 || $5 < 30) && length($10) > 10) {print $0}}' \
        down_genome_mapped.sam | samtools fastq | \
        awk '{if ($1 !~ /^@.*/) {print $0}
            if ($1 ~ /^@.*/) {split($0, read, "*"); print read[1] "/" read[2]}}' \
            > down_genome_mapped_bad.fastq
    awk '{if ($1 ~ /^@.*/) {print $0} 
        if ($1 !~ /^@.*/ && ($2 == 4 || $5 < 30) && length($10) > 10) {print $0}}' \
        up_genome_mapped.sam | samtools fastq | \
        awk '{if ($1 !~ /^@.*/) {print $0}
            if ($1 ~ /^@.*/) {split($0, read, "*"); print read[1] "/" read[2]}}' \
            > up_genome_mapped_bad.fastq
# Filter reads only with MAPQ >= 30
    awk '{if ($1 ~ /^@.*/) {print $0} 
        if ($1 !~ /^@.*/ && $5 >= 30) {print $0}}' \
        down_genome_mapped.sam > down_genome_mapped_good.sam
    awk '{if ($1 ~ /^@.*/) {print $0} 
        if ($1 !~ /^@.*/ && $5 >= 30) {print $0}}' \
        up_genome_mapped.sam > up_genome_mapped_good.sam
# find mates for low quality mapping reads
    perl $SCRIPTS_DIR/fastq_matefinder.pl down_genome_mapped_bad.fastq discordant.fastq
    perl $SCRIPTS_DIR/fastq_matefinder.pl up_genome_mapped_bad.fastq discordant.fastq
    seqtk seq -r up_genome_mapped_bad_paired.fastq > up_genome_mapped_bad_pairedRC.fastq
    rm up_genome_mapped_bad_paired.fastq    
done


#################################### stop halfway after rate-limiting step to re-run
for sample in $SAMPLES
do
    cd $WORK_DIR/$sample
    echo $sample
# Use sampe to map reads with mates and "rescue" some reads
    bwa aln -t 20 $GENOME down_genome_mapped_bad_mates.fastq > down_genome_mapped_bad_mates.sai
    bwa aln -t 20 $GENOME down_genome_mapped_bad_paired.fastq > down_genome_mapped_bad_paired.sai
    bwa aln -t 20 $GENOME up_genome_mapped_bad_mates.fastq > up_genome_mapped_bad_mates.sai
    bwa aln -t 20 $GENOME up_genome_mapped_bad_pairedRC.fastq > up_genome_mapped_bad_pairedRC.sai
    bwa sampe $GENOME \
        down_genome_mapped_bad_mates.sai down_genome_mapped_bad_paired.sai \
        down_genome_mapped_bad_mates.fastq down_genome_mapped_bad_paired.fastq \
        > down_genome_bad_remappedSAMPE.sam
    samtools view -f 2 -h down_genome_bad_remappedSAMPE.sam > down_genome_bad_remappedSAMPE2.sam
    samtools fastq -f 128 -F 2 down_genome_bad_remappedSAMPE.sam > down_genome_bad_paired_unrescued.fastq # just read 1 that don't map in proper pair
    bwa sampe $GENOME \
        up_genome_mapped_bad_mates.sai up_genome_mapped_bad_pairedRC.sai \
        up_genome_mapped_bad_mates.fastq up_genome_mapped_bad_pairedRC.fastq \
        > up_genome_bad_remappedSAMPE.sam
    samtools view -f 2 -h up_genome_bad_remappedSAMPE.sam > up_genome_bad_remappedSAMPE2.sam
    samtools fastq -f 128 -F 2 up_genome_bad_remappedSAMPE.sam > up_genome_bad_paired_unrescued.fastq # just read 1 that don't map in proper pair
    awk '{if ($1 ~ /^@/) print $0;                              # print @ headers
        else {if(and($2,0x40) && $5 >= 30 ) {                    # filter when read 1 has MAPQ > minimun
            getline; if (and($2,0x80)) {print $0}               # check next line is read 2 and print
            }
        }
    }' up_genome_bad_remappedSAMPE2.sam > up_genome_rescuedSAMPE.sam
    awk '{if ($1 ~ /^@/) print $0;                              # print @ headers
        else {if(and($2,0x40) && $5 < 30) {                     # filter when read 1 has MAPQ > minimun
            getline; if (and($2,0x80)) {print $0}               # check next line is read 2 and print
            }
        }
    }' up_genome_bad_remappedSAMPE2.sam | samtools fastq > up_genome_unrescued_quality.fastq
    awk '{if ($1 ~ /^@/) print $0;                              # print @ headers
        else {if(and($2,0x40) && $5 >= 30 ) {                    # filter when read 1 has MAPQ > minimum
            getline; if (and($2,0x80)) {print $0}               # check next line is read 2 and print
            }
        }
    }' down_genome_bad_remappedSAMPE2.sam > down_genome_rescuedSAMPE.sam
    awk '{if ($1 ~ /^@/) print $0;                              # print @ headers
        else {if(and($2,0x40) && $5 < 30) {                     # filter when read 1 has MAPQ > minimun
            getline; if (and($2,0x80)) {print $0}               # check next line is read 2 and print
            }
        }
    }' down_genome_bad_remappedSAMPE2.sam  | samtools fastq > down_genome_unrescued_quality.fastq
    rm up_genome_bad_remappedSAMPE.sam
    rm up_genome_bad_remappedSAMPE2.sam
    rm down_genome_bad_remappedSAMPE.sam
    rm down_genome_bad_remappedSAMPE2.sam
# Generate BED-like files from SAM files
    perl $SCRIPTS_DIR/intflanktorefSAMtoBEDplus3.pl  up_genome_mapped_good.sam down_genome_mapped_good.sam
    mv intreads.bed intreads_good.bed
    # reverse up direction mappings for correct output (were reverse complemented previously)
        awk '{if ($1 ~ /^@/) print $0;                              # print @ headers
            else {
                $2=xor($2, 0x10);                                   # flip whether read or pair is in reverse direction
                $2=xor($2, 0x20);
                print $0
            }
        }' OFS='\t' up_genome_rescuedSAMPE.sam > up_genome_rescuedSAMPE_corrected.sam
    perl $SCRIPTS_DIR/intflanktorefSAMtoBEDplus3.pl  up_genome_rescuedSAMPE_corrected.sam down_genome_rescuedSAMPE.sam
    #perl $SCRIPTS_DIR/intflanktorefSAMtoBEDplus3.pl  up_genome_rescuedSAMPE.sam down_genome_rescuedSAMPE.sam
    mv intreads.bed intreads_rescued.bed
    rm *.sai
# Compile non-rescued reads and map to repeat library
    cat down_genome_bad_paired_unrescued.fastq down_genome_mapped_bad_unpaired.fastq \
        down_genome_unrescued_quality.fastq > down_unmapped.fastq
    cat up_genome_bad_paired_unrescued.fastq up_genome_mapped_bad_unpaired.fastq \
        up_genome_unrescued_quality.fastq > up_unmapped.fastq
    bwa mem $REPEAT_LIB down_unmapped.fastq > down_repeatlibrary.sam
    bwa mem $REPEAT_LIB up_unmapped.fastq > up_repeatlibrary.sam
    samtools view -F 4 -q 30 -h down_repeatlibrary.sam > down_repeatlibrary_mapped.sam
    samtools view -F 4 -q 30 -h up_repeatlibrary.sam > up_repeatlibrary_mapped.sam
    perl $SCRIPTS_DIR/intflanktorefSAMtoBEDplus3.pl up_repeatlibrary_mapped.sam down_repeatlibrary_mapped.sam
    mv intreads.bed intreads_repeats.bed
    cat intreads_good.bed intreads_rescued.bed intreads_repeats.bed > intreads2.bed
# Change information for insertions found in refrence genome
    awk '{if($9=="Mar.3.4.6.p1_scaffold0_9821207_F") print $1,"9816237","9816242",$4,$5,$6,$7,"9816238","Mar.3.4.6.p1_scaffold0_9816238_F";
        else if($9=="Mar.3.4.6.p1_scaffold11_53805234_R") print $1,"53810202","53810207",$4,$5,$6,$7,"53810207","Mar.3.4.6.p1_scaffold11_53810207_R";
        else if($9=="Mar.3.4.6.p1_scaffold3_74958302_F") print $1,"74953332","74953337",$4,$5,$6,$7,"74953333","Mar.3.4.6.p1_scaffold3_74953333_F";
        else if($9=="Mar.3.4.6.p1_scaffold8_47484110_F") print $1,"47479138","47479143",$4,$5,$6,$7,"47479139","Mar.3.4.6.p1_scaffold8_47479139_F";
        else if($9=="Mar.3.4.6.p1_scaffold8_61542916_F") print $1,"61537946","61537951",$4,$5,$6,$7,"61537947","Mar.3.4.6.p1_scaffold8_61537947_F";
        else if($9=="Mar.3.4.6.p1_scaffold9_23229502_R") print $1,"23229677","23229682",$4,$5,$6,$7,"23229682","Mar.3.4.6.p1_scaffold9_23229682_R";
        else print $0}'  OFS='\t' intreads2.bed > intreads.bed
    rm intreads2.bed
done


for sample in $SAMPLES
do
# Run R script to aggregate read counts for each insertion site
    cd $WORK_DIR/$sample
    Rscript $SCRIPTS_DIR/Steamer_BED_analysis3.R 
    cp insertion_sites_updown.bed $WORK_DIR/$sample"_updown.bed"
    cp insertion_sites_multiple.bed $WORK_DIR/$sample"_multiple.bed"
done

# Summrize output
# for sample in MELC-2E11 MELC-A11_S1 PEI-DN08_S3 PEI-DF490 MELC-A9 NYTC-C9_S2 MELC-A10
echo "SUMMARY OF COUNTS" > $WORK_DIR/summary_counts
for sample in $SAMPLES
do
    cd $WORK_DIR/$sample
    echo $sample >> $WORK_DIR/summary_counts
    awk '$6 == "+" {plus++}
        $6 == "-" {minus++}
        $7 == "up" {up++}
        $7 == "down" {down++}
        END {print "genome: " plus "=plus, " minus "=minus, " up "=up, " down "=down" }' intreads_good.bed >> $WORK_DIR/summary_counts
    awk '$6 == "+" {plus++}
        $6 == "-" {minus++}
        $7 == "up" {up++}
        $7 == "down" {down++}
        END {print "rescued: " plus "=plus, " minus "=minus, " up "=up, " down "=down" }' intreads_rescued.bed >> $WORK_DIR/summary_counts
    awk '$6 == "+" {plus++}
        $6 == "-" {minus++}
        $7 == "up" {up++}
        $7 == "down" {down++}
        END {print "repeats: " plus "=plus, " minus "=minus, " up "=up, " down "=down" }' intreads_repeats.bed >> $WORK_DIR/summary_counts
    awk '$6 == "+" {plus++}
        $6 == "-" {minus++}
        $7 == "up" {up++}
        $7 == "down" {down++}
        END {print "ALL: " plus "=plus, " minus "=minus, " up "=up, " down "=down" }' intreads.bed >> $WORK_DIR/summary_counts
done


################
# Calculate coverage at steamer insertion sites
FILES=/ssd3/Mar_genome_analysis/steamer/final_pipeline
SAMPLES="MELC-2E11 MELC-A9 PEI-DF490 FFM-19G1 FFM-20B2 FFM-22F10 MELC-A11_S1 NYTC-C9_S2 DF-488 DN-HL03 PEI-DN08_S3"
BAMFILES=/ssd3/Mar_genome_analysis/bwa_mapping/Mar.3.4.6.p1/all_samples
INPUT=/ssd3/Mar_genome_analysis/steamer/final_pipeline
OUTPUT=/ssd3/Mar_genome_analysis/steamer/final_pipeline/coverage
TYPES="multiple updown"
cd $INPUT
for type in $TYPES
do
    echo $type
    for sample in $SAMPLES
    do
        if [ $sample = "MELC-2E11" ]; then number="01.MELC-2E11.bam"; fi
        if [ $sample = "MELC-A9" ]; then number="02.MELC-A9.bam"; fi
        if [ $sample = "PEI-DF490" ]; then number="03.PEI-DF490.bam"; fi
        if [ $sample = "FFM-19G1" ]; then number="08.FFM-19G1.bam"; fi
        if [ $sample = "FFM-20B2" ]; then number="09.FFM-20B2.bam"; fi
        if [ $sample = "FFM-22F10" ]; then number="11.FFM-22F10.bam"; fi
        if [ $sample = "MELC-A11_S1" ]; then number="13.MELC-A11_S1.bam"; fi
        if [ $sample = "NYTC-C9_S2" ]; then number="14.NYTC-C9_S2.bam"; fi
        if [ $sample = "DF-488" ]; then number="04.PEI-DF488.bam"; fi
        if [ $sample = "DN-HL03" ]; then number="05.PEI-DN03.bam"; fi
        if [ $sample = "PEI-DN08_S3" ]; then number="07.PEI-DN08_S3.bam"; fi
        echo $number started
        awk '$1 ~ /^M.*/ {chr=$1; start=$2; end=$3; name=$5; print chr, start-10, start, name}' OFS='\t' $INPUT/$sample"_"$type".bed" > $OUTPUT/sites_up.bed   # window 10bp up and downstream insertion
        awk '$1 ~ /^M.*/ {chr=$1; start=$2; end=$3; name=$5; print chr, end, end+10, name}' OFS='\t' $INPUT/$sample"_"$type".bed" > $OUTPUT/sites_dn.bed   # window 10bp up and downstream insertion
        samtools bedcov -Q 30 $OUTPUT/sites_up.bed $BAMFILES/$number | awk '{count = $5/10; print $4, count}' OFS='\t' > $OUTPUT/$sample"_"$type"_depthUP.bed"
        samtools bedcov -Q 30 $OUTPUT/sites_dn.bed $BAMFILES/$number | awk '{count = $5/10; print $4, count}' OFS='\t' > $OUTPUT/$sample"_"$type"_depthDN.bed"
        join $OUTPUT/$sample"_"$type"_depthUP.bed" $OUTPUT/$sample"_"$type"_depthDN.bed" | awk '{print $1, ($2 + $3)/2}' OFS='\t' > $OUTPUT/$sample"_"$type"_depth.bed"
        rm $OUTPUT/$sample"_"$type"_depthUP.bed"
        rm $OUTPUT/$sample"_"$type"_depthDN.bed"
    done
done

rm $OUTPUT/sites_up.bed
rm $OUTPUT/sites_dn.bed

