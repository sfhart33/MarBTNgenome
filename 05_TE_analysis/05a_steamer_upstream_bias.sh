GENES=/ssd3/Mar_genome_analysis/genomes/maker/2020-09-11_Mar_genome_snap02.all.genes.gff
UPSTREAM=/ssd3/Mar_genome_analysis/genomes/maker/2020-09-11_Mar_genome_snap02.all.upstream10000.gff
UPSTREAM2=/ssd3/Mar_genome_analysis/genomes/maker/2020-09-11_Mar_genome_snap02.all.upstream1000.gff
UPSTREAM3=/ssd3/Mar_genome_analysis/genomes/maker/2020-09-11_Mar_genome_snap02.all.upstream2000.gff
CDS=/ssd3/Mar_genome_analysis/genomes/maker/2020-09-11_Mar_genome_snap02.all.CDS.gff
GENOMEF=/ssd3/Mar_genome_analysis/genomes/maker/Mar.3.4.6.p1_Q30Q30A.genome
FIVEP=/ssd3/Mar_genome_analysis/genomes/maker/2020-09-11_Mar_genome_snap02.all.five_prime_UTR.gff
INPUT=/ssd3/Mar_genome_analysis/steamer/final_pipeline/merge
OUTPUT=/ssd3/Mar_genome_analysis/steamer/final_pipeline/genes
# STEAMER=all_sites.bed
module load bedtools/
cd $OUTPUT

# Test which upstream region has highest enrichment: 1000 and 2000bp
SIZES="100 1000 2000 2500 3000 4000 5000 10000"
FILES="all_sites anyCnoH allCnoH_allPEI_allUSA allCnoH allUSA allPEI"
echo -e "SAMPLES\tSIZEupstream\tSTEAMERins\BP" > size_tests.txt
for bed in $FILES
do
    for size in $SIZES
    do
        COUNTING=$(bedtools intersect -u -a $INPUT/$bed.bed -b /ssd3/Mar_genome_analysis/genomes/maker/"2020-09-11_Mar_genome_snap02.all.upstream"$size".gff" | wc -l)
        REGION_SIZE=$(bedtools sort -i /ssd3/Mar_genome_analysis/genomes/maker/"2020-09-11_Mar_genome_snap02.all.upstream"$size".gff" | bedtools merge | awk '{width = $3-$2; total += width} END{print total}'	)
        echo -e "$bed\t$size\t$COUNTING\t$REGION_SIZE" >> size_tests.txt
    done
done

FILES="all_sites anyCnoH allCnoH_allPEI_allUSA allCnoH allUSA allPEI"
echo -e "Steamer\tgenome\tgenes\tupstream10000\tupstream1000\tupstream2000\tcds\tfiveprime" > regions_intersect.txt
for bed in $FILES
do
    S_TOTAL=$(cat $INPUT/$bed.bed | wc -l)
    S_GENES=$(bedtools intersect -u -a $INPUT/$bed.bed -b $GENES | wc -l)
    S_UP=$(bedtools intersect -u -a $INPUT/$bed.bed -b $UPSTREAM | wc -l)
    S_UP2=$(bedtools intersect -u -a $INPUT/$bed.bed -b $UPSTREAM2 | wc -l)
    S_UP3=$(bedtools intersect -u -a $INPUT/$bed.bed -b $UPSTREAM3 | wc -l)
    S_CDS=$(bedtools intersect -u -a $INPUT/$bed.bed -b $CDS | wc -l)
    S_FIVEP=$(bedtools intersect -u -a $INPUT/$bed.bed -b $FIVEP | wc -l)
    echo -e "$bed\t$S_TOTAL\t$S_GENES\t$S_UP\t$S_UP2\t$S_UP3\t$S_CDS\t$S_FIVEP" >> regions_intersect.txt
done

G_TOTAL=$(awk '{total += $2} END{print total}' $GENOMEF	)
G_GENES=$(bedtools sort -i $GENES | bedtools merge | awk '{width = $3-$2; total += width} END{print total}')
G_UP=$(bedtools sort -i $UPSTREAM | bedtools merge | awk '{width = $3-$2; total += width} END{print total}' )
G_UP2=$(bedtools sort -i $UPSTREAM2 | bedtools merge | awk '{width = $3-$2; total += width} END{print total}'	)
G_UP3=$(bedtools sort -i $UPSTREAM3 | bedtools merge | awk '{width = $3-$2; total += width} END{print total}'	)
G_CDS=$(bedtools sort -i $CDS | bedtools merge | awk '{width = $3-$2; total += width} END{print total}' )
G_FIVEP=$(bedtools sort -i $FIVEP | bedtools merge | awk '{width = $3-$2; total += width} END{print total}'	)
echo -e "SIZE\t$G_TOTAL\t$G_GENES\t$G_UP\t$G_UP2\t$G_UP3\t$G_CDS\t$G_FIVEP" >> regions_intersect.txt

head regions_intersect.txt

# generate proximity files to analyze in R
    for bed in $FILES
    do
        bedtools sort -i $GENES | \
            bedtools closest -D b -io -a $INPUT/$bed.bed -b - | \
            awk '{split($15, id, ";" ); split(id[1], name, "=" ); print $5, name[2], $16}'  OFS='\t' > $bed.proximity
        bedtools sort -i $GENES | \
            bedtools closest -D b -io -s -a $INPUT/$bed.bed -b - | \
            awk '{split($15, id, ";" ); split(id[1], name, "=" ); print $5, name[2], $16}'  OFS='\t' > $bed.proximity.samestrand
        bedtools sort -i $GENES | \
            bedtools closest -D b -io -S -a $INPUT/$bed.bed -b - | \
            awk '{split($15, id, ";" ); split(id[1], name, "=" ); print $5, name[2], $16}'  OFS='\t' > $bed.proximity.difstrand
    done

# generate gene list files to compare to COSMIC genes
    FILE_SUBSET="allCnoH_allPEI_allUSA allCnoH allUSA allPEI"
    for bed in $FILE_SUBSET
    do
        awk '$3 > -2000 && $3 < 0 {print $2}'  OFS='\t' $bed.proximity > $bed.2000bp_upstream
        awk '$3 > -2000 && $3 < 0 {print $2}'  OFS='\t' $bed.proximity.samestrand > $bed.samestrand.2000bp_upstream
        bedtools intersect -wb -a $INPUT/$bed.bed -b $CDS | \
            awk '{split($15, id, ";" ); split(id[1], name, "=" ); print name[2]}'  OFS='\t' > $bed.cds
        bedtools intersect -wb -a $INPUT/$bed.bed -b $GENES | \
            awk '{split($15, id, ";" ); split(id[1], name, "=" ); print name[2]}'  OFS='\t' > $bed.genes
    done

# Test to correct for possibilty of read mapping bias toward upstream genes 
    TEST_BAM=/ssd3/Mar_genome_analysis/bwa_mapping/Mar.3.4.6.p1/all_samples/01.MELC-2E11.bam
    # Take only one out of every 1,000 reads, only first in pair and generate fake steamer insertion file
        samtools view -q 30 -@ 50 -f 64 -s 0.001 $TEST_BAM > MELC-2E11.subset.sam 
            awk '{end=$4+5; print $3, $4, end}' OFS='\t' MELC-2E11.subset.sam  | \
            bedtools sort -i - > MELC-2E11.subset.bed
        wc -l MELC-2E11.subset.bed 
        bedtools sort -i $GENES | \
            bedtools closest -D b -io -a MELC-2E11.subset.bed -b - | \
            awk '{split($12, id, ";" ); split(id[1], name, "=" ); print $1":"$2"-"$3, name[2], $13}'  OFS='\t' > MELC-2E11.subset.proximity
    # count intersections between fake file and regions
        echo -e "genome\tgenes\tupstream10000\tupstream1000\tupstream2000\tcds\tfiveprime" > regions_intersect_tests.txt
        S_TOTAL=$(cat MELC-2E11.subset.bed | wc -l)
        S_GENES=$(bedtools intersect -u -a MELC-2E11.subset.bed -b $GENES | wc -l)
        S_UP=$(bedtools intersect -u -a MELC-2E11.subset.bed -b $UPSTREAM | wc -l)
        S_UP2=$(bedtools intersect -u -a MELC-2E11.subset.bed -b $UPSTREAM2 | wc -l)
        S_UP3=$(bedtools intersect -u -a MELC-2E11.subset.bed -b $UPSTREAM3 | wc -l)
        S_CDS=$(bedtools intersect -u -a MELC-2E11.subset.bed -b $CDS | wc -l)
        S_FIVEP=$(bedtools intersect -u -a MELC-2E11.subset.bed -b $FIVEP | wc -l)
        echo -e "INTERSECTIONS\t$S_TOTAL\t$S_GENES\t$S_UP\t$S_UP2\t$S_UP3\t$S_CDS\t$S_FIVEP" >> regions_intersect_tests.txt
        G_TOTAL=$(awk '{total += $2} END{print total}' $GENOMEF	)
        G_GENES=$(bedtools sort -i $GENES | bedtools merge | awk '{width = $3-$2; total += width} END{print total}')
        G_UP=$(bedtools sort -i $UPSTREAM | bedtools merge | awk '{width = $3-$2; total += width} END{print total}' )
        G_UP2=$(bedtools sort -i $UPSTREAM2 | bedtools merge | awk '{width = $3-$2; total += width} END{print total}'	)
        G_UP3=$(bedtools sort -i $UPSTREAM3 | bedtools merge | awk '{width = $3-$2; total += width} END{print total}'	)
        G_CDS=$(bedtools sort -i $CDS | bedtools merge | awk '{width = $3-$2; total += width} END{print total}' )
        G_FIVEP=$(bedtools sort -i $FIVEP| bedtools merge | awk '{width = $3-$2; total += width} END{print total}'	)
        echo -e "SIZE\t$G_TOTAL\t$G_GENES\t$G_UP\t$G_UP2\t$G_UP3\t$G_CDS\t$G_FIVEP" >> regions_intersect_tests.txt
        head regions_intersect_tests.txt

# compare to ATG rather than gene
# generated in R from dndscv buildref: start_codons.bed
    ATG_SITES=/ssd3/Mar_genome_analysis/steamer/final_pipeline/genes/start_codons_atg.bed
    START_SITES=/ssd3/Mar_genome_analysis/steamer/final_pipeline/genes/start_codons.bed
    OUTPUT=/ssd3/Mar_genome_analysis/steamer/final_pipeline/genes
    # check to make sure I got ATG sites

        bedtools getfasta -s -fi /ssd3/Mar_genome_analysis/genomes/Mar.3.4.6.p1_Q30Q30A.fasta -bed $START_SITES | head -100
        bedtools getfasta -s -fi /ssd3/Mar_genome_analysis/genomes/Mar.3.4.6.p1_Q30Q30A.fasta -bed $ATG_SITES | head -20
        bedtools getfasta -s -fi /ssd3/Mar_genome_analysis/genomes/Mar.3.4.6.p1_Q30Q30A.fasta -bed $ATG_SITES | tail -20
    # generate proximity files to analyze in R
    #FILE_SUBSET="allCnoH allUSA allPEI"
    for bed in allCnoH_allPEI_allUSA
    do
        for codon in start_codons_atg start_codons
        do
            bedtools sort -i $OUTPUT/$codon.bed | \
                bedtools closest -D b -io -a $INPUT/$bed.bed -b - | \
                awk '{print $4, $10, $13}'  OFS='\t' > $bed.$codon.proximity
            bedtools sort -i $OUTPUT/$codon.bed | \
                bedtools closest -D b -io -s -a $INPUT/$bed.bed -b - | \
                awk '{print $4, $10, $13}'  OFS='\t' > $bed.$codon.proximity.samestrand
            bedtools sort -i $OUTPUT/$codon.bed | \
                bedtools closest -D b -io -S -a $INPUT/$bed.bed -b - | \
                awk '{print $4, $10, $13}'  OFS='\t' > $bed.$codon.proximity.difstrand
        done
    done

