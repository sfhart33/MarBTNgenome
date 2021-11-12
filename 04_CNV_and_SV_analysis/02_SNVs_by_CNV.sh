    module load bedtools
    SNPS=/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/bins
    # allPEInoUSAnoH_freq.bed
    # allUSAnoPEInoH_freq.bed
    CNV=/ssd3/Mar_genome_analysis/CNV_calling/FINAL/output_bed
    LOH=/ssd3/Mar_genome_analysis/LOH/july_2021/output/new
    # USA_LOH_10_hetero_counts.bed.merge
    # PEI_LOH_10_hetero_counts.bed.merge
    OUTPUT=/ssd3/Mar_genome_analysis/CNV_calling/FINAL/SNVs
    RANGE="0 1 2 3 4 5 6 7plus"
    FILES_USA="allUSAnoPEInoH allCanyH allCnoH"
    FILES_PEI="allPEInoUSAnoH allCanyH allCnoH"
    for method in "" "_noH0-1"
    do
        for i in $RANGE
        do
            for file in $FILES_PEI
            do
                bedtools intersect -v -a $SNPS/$file"_freq.bed" -b $LOH/PEI_LOH_10_hetero_counts.bed.merge | \
                    bedtools intersect -a - -b $CNV/"PEI"$method"_CN"$i".bed" > $OUTPUT/"PEI_"$file"_noLOH"$method"_CN"$i".bed"
                bedtools intersect -a $SNPS/$file"_freq.bed" -b $LOH/PEI_LOH_10_hetero_counts.bed.merge | \
                    bedtools intersect -a - -b $CNV/"PEI"$method"_CN"$i".bed" > $OUTPUT/"PEI_"$file"_LOH"$method"_CN"$i".bed"
                bedtools intersect -a $SNPS/$file"_freq.bed" -b $CNV/"PEI"$method"_CN"$i".bed" > $OUTPUT/"PEI_"$file$method"_CN"$i".bed"
            done
            for file in $FILES_USA
            do
                bedtools intersect -v -a $SNPS/$file"_freq.bed" -b $LOH/USA_LOH_10_hetero_counts.bed.merge | \
                    bedtools intersect -a - -b $CNV/"USA"$method"_CN"$i".bed" > $OUTPUT/"USA_"$file"_noLOH"$method"_CN"$i".bed"
                bedtools intersect -a $SNPS/$file"_freq.bed" -b $LOH/USA_LOH_10_hetero_counts.bed.merge | \
                    bedtools intersect -a - -b $CNV/"USA"$method"_CN"$i".bed" > $OUTPUT/"USA_"$file"_LOH"$method"_CN"$i".bed"
                bedtools intersect -a $SNPS/$file"_freq.bed" -b $CNV/"USA"$method"_CN"$i".bed" > $OUTPUT/"USA_"$file$method"_CN"$i".bed"
            done
        done
    done