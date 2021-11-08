# Original file here: C:\Users\shart\Metzger Lab Dropbox\Sam_data\Mya_genome\LOH\somatypus_output_LOH_merge-helmsman.sh
module load bedtools/

    cd /ssd3/Mar_genome_analysis/LOH/july_2021
    BEDS=$(ls *.bed)
    OUTPUT=/ssd3/Mar_genome_analysis/LOH/july_2021/output/new
    SNVS=/ssd3/Mar_genome_analysis/LOH/july_2021/output

# Merge overlapping LOH regions, count size
    echo -e "Sample\tRegions\tbp" > $OUTPUT/LOH.count.test.txt
    for bed in $BEDS
    do
        echo $bed
            bedtools merge -d 0 -i $bed > $OUTPUT/$bed.merge
    done
    for bed in $BEDS
    do
        echo $bed
        awk -v sample=$bed '{width = $3 - $2 - 1; widthsum += width} END{print sample, NR, widthsum}' OFS='\t' $OUTPUT/$bed.merge >> $OUTPUT/LOH.count.test.txt
    done
    rm $BEDS

# Make helmsman input file
# Note we want to look at PEI SNVs in USA LOH regions, as those will likely include founder germline SNV contamination
    cd $OUTPUT
    BEDS_USA=$(ls USA*.merge)
    BEDS_PEI=$(ls PEI*.merge)
    #MERGES="0 10000"
    rm helmsman.txt
    for bed in $BEDS_USA
    do
        # for merge in $MERGES
        # do
            echo $bed
            bedtools intersect -wa -a $SNVS/allPEInoUSAnoH.bed -b $OUTPUT/$bed | \
                awk -v name=$bed.PEIloh '{print $1, $3, $4, $5, name}' OFS='\t' >> $OUTPUT/helmsman.txt
            bedtools intersect -v -wa -a $SNVS/allPEInoUSAnoH.bed -b $OUTPUT/$bed | \
                awk -v name=$bed.PEIsom '{print $1, $3, $4, $5, name}' OFS='\t' >> $OUTPUT/helmsman.txt
        # done
    done
    for bed in $BEDS_PEI
    do
        # for merge in $MERGES
        # do
            echo $bed
            bedtools intersect -wa -a $SNVS/allUSAnoPEInoH.bed -b $OUTPUT/$bed | \
                awk -v name=$bed.USAloh '{print $1, $3, $4, $5, name}' OFS='\t' >> $OUTPUT/helmsman.txt
            bedtools intersect -v -wa -a $SNVS/allUSAnoPEInoH.bed -b $OUTPUT/$bed | \
                awk -v name=$bed.USAsom '{print $1, $3, $4, $5, name}' OFS='\t' >> $OUTPUT/helmsman.txt
        # done
    done

# Run helmsman
    SNPS=/ssd3/Mar_genome_analysis/LOH/july_2021/output/new/helmsman.txt
    OUTPUT=/ssd3/Mar_genome_analysis/LOH/july_2021/output/new/helmsman
    GENOME=/ssd3/Mar_genome_analysis/genomes/Mar.3.4.6.p1_Q30Q30A.fasta
    module load conda/
    source activate helmsman
    python /ssd3/Mar_genome_analysis/mut_sig/helmsman/helmsman.py --verbose --decomp nmf --mode txt --input $SNPS --fastafile $GENOME --projectdir $OUTPUT --length 3


# SPLIT VCF FILES INTO LOH AND NON-LOH REGIONS
    WORK_DIR=/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/bins
    OUTPUT=/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/bins/LOH
    LOH=/ssd3/Mar_genome_analysis/LOH/july_2021/output/new
    SNVs="SNVs Indels"
    REGIONS="USA PEI"
    COUNTS=$(echo {1..50}) #"5 25" # $(echo {1..50}) #
    STRINGENCY="regular stringent"
    cd $WORK_DIR

    for count in $COUNTS
    do
    for region in $REGIONS
    do
    cp $LOH/$region"_LOH_"$count"_hetero_counts.bed.merge" $LOH/$region"_hetero_LOH_"$count"_counts.bed"
    done
    done

# split regions into new bed files
    for region in $REGIONS
    do
        echo $region
        if [ $region = "USA" ]; then
            otherregion="PEI"
        else
            otherregion="USA"
        fi
        echo $otherregion
        ###
        for count in $COUNTS
        do
            echo "   "$count
            for stringency in $STRINGENCY
            do
            echo "      "$stringency
            bedtools intersect -wa -a $WORK_DIR/"all"$otherregion"no"$region"noH_"$stringency".bed" -b $LOH/$region"_hetero_LOH_"$count"_counts.bed" | \
                awk -v name="all"$otherregion"no"$region"noH_"$stringency"_"$count"_LOH" '{print name, $1, $3, $4, $5}' OFS='\t' > $OUTPUT/"all"$otherregion"no"$region"noH_"$stringency"_"$count"_LOH.dnds" &
            bedtools intersect -v -wa -a $WORK_DIR/"all"$otherregion"no"$region"noH_"$stringency".bed" -b $LOH/$region"_hetero_LOH_"$count"_counts.bed" | \
                awk -v name="all"$otherregion"no"$region"noH_"$stringency"_"$count"_nonLOH" '{print name, $1, $3, $4, $5}' OFS='\t' > $OUTPUT/"all"$otherregion"no"$region"noH_"$stringency"_"$count"_nonLOH.dnds" &
            done
        done
    done

# count all regions
    echo -e "loh_region\tsnv_region\tstringency\tcount\tLOH_SNVs\tnonLOH_SNVs\tLOH_SIZE" > $WORK_DIR/LOH_abberent_SNVs_counts.txt
    for region in $REGIONS
    do
        echo $region
        if [ $region = "USA" ]; then
            otherregion="PEI"
        else
            otherregion="USA"
        fi
        echo $otherregion
        ###
        for count in $COUNTS
        do
            echo "   count = "$count
            for stringency in $STRINGENCY
            do
            echo "      "$stringency
            LOH_SNVs=$(bedtools intersect -wa -a $WORK_DIR/"all"$otherregion"no"$region"anyH_"$stringency".bed" -b $LOH/$region"_hetero_LOH_"$count"_counts.bed" | \
                awk -v name="all"$otherregion"no"$region"anyH_"$stringency"_"$count"_LOH" '{print name, $1, $3, $4, $5}' OFS='\t' | \
                wc -l)
            nonLOH_SNVs=$(bedtools intersect -v -wa -a $WORK_DIR/"all"$otherregion"no"$region"anyH_"$stringency".bed" -b $LOH/$region"_hetero_LOH_"$count"_counts.bed" | \
                awk -v name="all"$otherregion"no"$region"anyH_"$stringency"_"$count"_nonLOH" '{print name, $1, $3, $4, $5}' OFS='\t'  | \
                wc -l)
            LOH_SIZE=$(awk '{width = $3 - $2 - 1; widthsum += width} END{print widthsum}' $LOH/$region"_hetero_LOH_"$count"_counts.bed")
            echo -e $region"\t"$otherregion"\t"$stringency"\t"$count"\t"$LOH_SNVs"\t"$nonLOH_SNVs"\t"$LOH_SIZE # >> $WORK_DIR/LOH_abberent_SNVs_counts.txt
            done
        done
    done




# Get rid of extra columns of the vcf files to make a more compact version and compress output
    WORK_DIR=/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run
    for snvs in $SNVs
    do
        #testing
        # snvs="SNVs"
        # head -1000 $WORK_DIR/"Somatypus_"$snvs"_final.vcf" > $WORK_DIR/"Somatypus_"$snvs"_final_1000.vcf"
        awk ' $0 ~ /^#.*/ {print} 
            $0 !~ /^#.*/ {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23}
            ' OFS='\t' $WORK_DIR/"Somatypus_"$snvs"_final.vcf" | \
            bgzip -c > $WORK_DIR/"Somatypus_"$snvs"_final_compact.vcf.gz"
        tabix -fp vcf $WORK_DIR/"Somatypus_"$snvs"_final_compact.vcf.gz"
    done

# split regions into new vcf files
    for snvs in $SNVs
    do
    for region in $REGIONS
    do
        bcftools filter -O v -T $LOH/$region"_LOH_5_counts.bed" $WORK_DIR/"Somatypus_"$snvs"_final_compact.vcf.gz" > $WORK_DIR/"Somatypus_"$snvs"_compact_"$region"_LOH_5.vcf" &
        bcftools filter -O v -T $LOH/$region"_LOH_10_counts.bed" $WORK_DIR/"Somatypus_"$snvs"_final_compact.vcf.gz" > $WORK_DIR/"Somatypus_"$snvs"_compact_"$region"_LOH_10.vcf" &
        bcftools filter -O v -T ^$LOH/$region"_LOH_5_counts.bed" $WORK_DIR/"Somatypus_"$snvs"_final_compact.vcf.gz" > $WORK_DIR/"Somatypus_"$snvs"_compact_"$region"_nonLOH_5.vcf" &
        bcftools filter -O v -T ^$LOH/$region"_LOH_10_counts.bed" $WORK_DIR/"Somatypus_"$snvs"_final_compact.vcf.gz" > $WORK_DIR/"Somatypus_"$snvs"_compact_"$region"_nonLOH_10.vcf" &
    done
    done
