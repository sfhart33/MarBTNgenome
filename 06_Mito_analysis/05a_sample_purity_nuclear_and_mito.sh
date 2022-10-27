WORKD=/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run
CNV=/ssd3/Mar_genome_analysis/CNV_calling/FINAL/output_bed
COMMAND=/ssd3/Mar_genome_analysis/somatypus
OUTPUT=/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/validation

cd $OUTPUT

bcftools filter -O b -T $CNV/PEI_CN2.bed $WORKD/Somatypus_SNVs_final.vcf.gz | \
    bcftools filter -O v -T $CNV/USA_CN2.bed - > $OUTPUT/SNVs_bothCN2b.vcf &
bcftools filter -O b -T $CNV/PEI_CN4.bed $WORKD/Somatypus_SNVs_final.vcf.gz | \
    bcftools filter -O v -T $CNV/USA_CN4.bed - > $OUTPUT/SNVs_bothCN4b.vcf &
wait

wc - l *vcf

awk '$0 !~ /^##contig=.*/ && $0 !~ /^##FILTER=<ID=PASS.*/ && $0 !~ /^##bcftools_filter.*/  {print}' SNVs_bothCN2b.vcf > SNVs_bothCN2.vcf &
awk '$0 !~ /^##contig=.*/ && $0 !~ /^##FILTER=<ID=PASS.*/ && $0 !~ /^##bcftools_filter.*/  {print}' SNVs_bothCN4b.vcf > SNVs_bothCN4.vcf &
wait

#python2 $COMMAND/ExtractVcfData.py SNVs_USA_noLOH.vcf &
python2 $COMMAND/ExtractVcfData.py SNVs_bothCN2.vcf &
#python2 $COMMAND/ExtractVcfData.py SNVs_PEI_noLOH.vcf &
python2 $COMMAND/ExtractVcfData.py SNVs_bothCN4.vcf &
wait

rm *.vcf
