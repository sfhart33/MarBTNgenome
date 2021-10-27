# MAP TO PUBLISHED MT GENOME AND INDEX

GENOME=/ssd3/Mar_genome_analysis/genomes/mito/mt_genome.fasta
LIST="MELC-A10 MELC-A11_S1 NYTC-C9_S2 PEI-DN08_S3 MELC-2E11 MELC-A9 PEI-DF490 DF-488 DN-HL03 FFM-19G1 FFM-20B2"
WORK=/ssd3/Mar_genome_analysis/bwa_mapping/mito/all_samples
ILLUMINA=/ssd2/Illumina_data/dedupe/trim20

cd $WORK
for i in $LIST
do
	# map to mt genome, convert to bam and only kep mapped reads, sort
	bwa mem -t 10 $GENOME $ILLUMINA/$i"_R1_001.fastq.gz" $ILLUMINA/$i"_R2_001.fastq.gz" | samtools view -b -h -F 4 -@ 1 | samtools sort -O bam -@ 1 > $i.mtgenome.bam &
done
wait


LIST="DN-HL07 FFM-22A10 FFM-22F10 FFM-22A10-adductor FFM-22F10-adductor MELC-A10-siphon MELC-A11-siphon NYTC-C9-mantle PEI-DN03-siphon PEI-DN07-siphon PEI-DN08-siphon"
ILLUMINA=/ssd2/Illumina_data/dedupe/2021_newsamples_trim20

cd $WORK
for i in $LIST
do
	# map to mt genome, convert to bam and only kep mapped reads, sort
	bwa mem -t 10 $GENOME $ILLUMINA/$i"_R1_001.fastq.gz" $ILLUMINA/$i"_R2_001.fastq.gz" | samtools view -b -h -F 4 -@ 1 | samtools sort -O bam -@ 1 > $i.mtgenome.bam &
done
wait

LIST="DF-488 DN-HL03 DN-HL07 FFM-19G1 FFM-20B2 FFM-22A10-adductor FFM-22A10 FFM-22F10-adductor FFM-22F10 MELC-2E11 MELC-A10 MELC-A10-siphon MELC-A11_S1 MELC-A11-siphon MELC-A9 NYTC-C9-mantle NYTC-C9_S2 PEI-DF490 PEI-DN03-siphon PEI-DN07-siphon PEI-DN08_S3 PEI-DN08-siphon"

for i in $LIST
do
	# index
	samtools index -b -@ 3 $i.mtgenome.bam &
done
wait

# RUN SOMATYPUS

module load somatypus/
OUTPUT=/ssd3/Mar_genome_analysis/bwa_mapping/mito/all_samples/somatypus
cd $WORK
somatypus -i $WORK -o $OUTPUT -g $GENOME -c 8

# NOTES ON VCF STRUCTURE:
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	DF-488.mtgenome	DN-HL03.mtgenome	DN-HL07.mtgenome	FFM-19G1.mtgenome	FFM-20B2.mtgenome	FFM-22A10-adductor.mtgenome	FFM-22A10.mtgenome	FFM-22F10-adductor.mtgenome	FFM-22F10.mtgenome	MELC-2E11.mtgenome	MELC-A10-siphon.mtgenome	MELC-A10.mtgenome	MELC-A11-siphon.mtgenome	MELC-A11_S1.mtgenome	MELC-A9.mtgenome	NYTC-C9-mantle.mtgenome	NYTC-C9_S2.mtgenome	PEI-DF490.mtgenome	PEI-DN03-siphon.mtgenome	PEI-DN07-siphon.mtgenome	PEI-DN08-siphon.mtgenome	PEI-DN08_S3.mtgenome
    # order of samples
    # $10 = DF-488.mtgenome     Cpei0
    # $11 = DN-HL03.mtgenome    Cpei1
    # $12 = DN-HL07.mtgenome    Cpei2
    # $13 = FFM-19G1.mtgenome   Cusa0a
    # $14 = FFM-20B2.mtgenome   Cusa0b
    # $15 = FFM-22A10-adductor.mtgenome Tusa1
    # $16 = FFM-22A10.mtgenome  Cusa1
    # $17 = FFM-22F10-adductor.mtgenome Tusa2
    # $18 =	FFM-22F10.mtgenome   Cusa2
    # $19 =	MELC-2E11.mtgenome  Href
    # $20 =	MELC-A10-siphon.mtgenome Tusa3	
    # $21 =	MELC-A10.mtgenome      Cusa3	
    # $22 =	MELC-A11-siphon.mtgenome Tusa4	
    # $23 =	MELC-A11_S1.mtgenome   Cusa4	
    # $24 =	MELC-A9.mtgenome	   Husa
    # $25 =	NYTC-C9-mantle.mtgenome	 Tusa5
    # $26 =	NYTC-C9_S2.mtgenome   Cusa5
    # $27 =	PEI-DF490.mtgenome	Hpei
    # $28 =	PEI-DN03-siphon.mtgenome	Tpei1
    # $29 =	PEI-DN07-siphon.mtgenome	Tpei2	
    # $30 =	PEI-DN08-siphon.mtgenome	Tpei3	
    # $31 =	PEI-DN08_S3.mtgenome    Cpei3


# MAKE SUMMARY FILE TO PLOT IN R

TYPES="SNVs Indels"
for type in $TYPES
do
awk ' $0 ~ /^#CHROM.*/ {print "chr","pos", 
	"Href_a", "Husa_a", "Hpei_a", 
	"Cpei0_a", "Cpei1_a", "Cpei2_a", "Cpei3_a",
	"Cusa0a_a", "Cusa0b_a", "Cusa1_a", "Cusa2_a", "Cusa3_a", "Cusa4_a", "Cusa5_a",
	"Tpei1_a", "Tpei2_a", "Tpei3_a",
	"Tusa1_a", "Tusa2_a", "Tusa3_a", "Tusa4_a", "Tusa5_a",
	"Href_t", "Husa_t", "Hpei_t", 
	"Cpei0_t", "Cpei1_t", "Cpei2_t", "Cpei3_t",
	"Cusa0a_t", "Cusa0b_t", "Cusa1_t", "Cusa2_t", "Cusa3_t", "Cusa4_t", "Cusa5_t",
	"Tpei1_t", "Tpei2_t", "Tpei3_t",
	"Tusa1_t", "Tusa2_t", "Tusa3_t", "Tusa4_t", "Tusa5_t",
	"Href_f", "Husa_f", "Hpei_f", 
	"Cpei0_f", "Cpei1_f", "Cpei2_f", "Cpei3_f",
	"Cusa0a_f", "Cusa0b_f", "Cusa1_f", "Cusa2_f", "Cusa3_f", "Cusa4_f", "Cusa5_f",
	"Tpei1_f", "Tpei2_f", "Tpei3_f",
	"Tusa1_f", "Tusa2_f", "Tusa3_f", "Tusa4_f", "Tusa5_f"}
    $0 !~ /^#.*/ {
        split($10, Cpei0, ":" );
        split($11, Cpei1, ":" );
        split($12, Cpei2, ":" );
        split($13, Cusa0a, ":" );
        split($14, Cusa0b, ":" );
        split($15, Tusa1, ":" );
        split($16, Cusa1, ":" );
        split($17, Tusa2, ":" );
        split($18, Cusa2, ":" );
        split($19, Href, ":" );
        split($20, Tusa3, ":" );
        split($21, Cusa3, ":" );
        split($22, Tusa4, ":" );
        split($23, Cusa4, ":" );
        split($24, Husa, ":" );
        split($25, Tusa5, ":" );
        split($26, Cusa5, ":" );
        split($27, Hpei, ":" );
        split($28, Tpei1, ":" );
        split($29, Tpei2, ":" );
        split($30, Tpei3, ":" );
        split($31, Cpei3, ":" ); 
        print $1, $2, 
	Href[6], Husa[6], Hpei[6], 
	Cpei0[6], Cpei1[6], Cpei2[6], Cpei3[6],
	Cusa0a[6], Cusa0b[6], Cusa1[6], Cusa2[6], Cusa3[6], Cusa4[6], Cusa5[6],
	Tpei1[6], Tpei2[6], Tpei3[6],
	Tusa1[6], Tusa2[6], Tusa3[6], Tusa4[6], Tusa5[6],
	Href[5], Husa[5], Hpei[5], 
	Cpei0[5], Cpei1[5], Cpei2[5], Cpei3[5],
	Cusa0a[5], Cusa0b[5], Cusa1[5], Cusa2[5], Cusa3[5], Cusa4[5], Cusa5[5],
	Tpei1[5], Tpei2[5], Tpei3[5],
	Tusa1[5], Tusa2[5], Tusa3[5], Tusa4[5], Tusa5[5],
	Href[6]/Href[5], Husa[6]/Husa[5], Hpei[6]/Hpei[5], 
	Cpei0[6]/Cpei0[5], Cpei1[6]/Cpei1[5], Cpei2[6]/Cpei2[5], Cpei3[6]/Cpei3[5],
            Cusa0a[6]/Cusa0a[5], Cusa0b[6]/Cusa0b[5], Cusa1[6]/Cusa1[5], Cusa2[6]/Cusa2[5], Cusa3[6]/Cusa3[5], Cusa4[6]/Cusa4[5], Cusa5[6]/Cusa5[5],
            Tpei1[6]/Tpei1[5], Tpei2[6]/Tpei2[5], Tpei3[6]/Tpei3[5],
            Tusa1[6]/Tusa1[5], Tusa2[6]/Tusa2[5], Tusa3[6]/Tusa3[5], Tusa4[6]/Tusa4[5], Tusa5[6]/Tusa5[5]
}'  OFS='\t'  $OUTPUT/"Somatypus_"$type"_final.vcf" > $OUTPUT/"Somatypus_"$type"_final.counts"
done


# CACULATE READ DEPTH ACROSS MT GENOME
INPUT=/ssd3/Mar_genome_analysis/bwa_mapping/mito/all_samples
OUTPUT=/ssd3/Mar_genome_analysis/bwa_mapping/mito/all_samples/depth
LIST="DF-488 DN-HL03 DN-HL07 FFM-19G1 FFM-20B2 FFM-22A10 FFM-22F10 MELC-2E11 MELC-A10 MELC-A11_S1 MELC-A9 NYTC-C9_S2 PEI-DF490 PEI-DN08_S3"
cd $OUTPUT
for i in $LIST
do
samtools depth -a -d 0 $INPUT/$i.mtgenome.bam > $OUTPUT/$i.mtgenome.depth
done

