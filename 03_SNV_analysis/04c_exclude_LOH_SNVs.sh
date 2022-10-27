# This was at the bottom of 04_run_somatypus
########################

module load R/3.6.0
R
# See here: C:\Users\shart\Metzger Lab Dropbox\Sam_data\Mya_genome\SNPs\somatypus_output_BINNING.r


#####################
# outside_LOH regions
WORKD=/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run
LOH=/ssd3/Mar_genome_analysis/LOH/july_2021/output/new
COMMAND=/ssd3/Mar_genome_analysis/somatypus
cd $WORKD
bcftools filter -O b -T ^$LOH/PEI_hetero_LOH_10_counts.bed Somatypus_SNVs_final.vcf.gz | \
	bcftools filter -O v -T ^$LOH/USA_hetero_LOH_10_counts.bed > Somatypus_SNVs_final_noLOH.vcf
awk '$0 !~ /^##contig=.*/ && $0 !~ /^##FILTER=<ID=PASS.*/ && $0 !~ /^##bcftools_filter.*/  {print}' Somatypus_SNVs_final_noLOH.vcf > Somatypus_SNVs_final_noLOH2.vcf
$COMMAND/ExtractVcfData.py Somatypus_SNVs_final_noLOH2.vcf
rm Somatypus_SNVs_final_noLOH.vcf
rm Somatypus_SNVs_final_noLOH2.vcf

cat $LOH/PEI_hetero_LOH_10_counts.bed $LOH/USA_hetero_LOH_10_counts.bed | sort -k1,1 -k2,2n -k3,3n | bedtools merge -d 0 -i - > $LOH/BOTH_SUBLINEAGES_LOH_10_counts.bed

########################
# Original file here: C:\Users\shart\Metzger Lab Dropbox\Sam_data\Mya_genome\LOH\LOH_revision_validation.sh

OUTPUT=/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/validation

cd $OUTPUT

#bcftools filter -O v -T ^$LOH/USA_hetero_LOH_10_counts.bed $WORKD/Somatypus_SNVs_final.vcf.gz > $OUTPUT/Somatypus_SNVs_USA_noLOH.vcf &
bcftools filter -O v -T $LOH/USA_hetero_LOH_10_counts.bed $WORKD/Somatypus_SNVs_final.vcf.gz > $OUTPUT/Somatypus_SNVs_USA_LOH.vcf &
#bcftools filter -O v -T ^$LOH/PEI_hetero_LOH_10_counts.bed $WORKD/Somatypus_SNVs_final.vcf.gz > $OUTPUT/Somatypus_SNVs_PEI_noLOH.vcf &
bcftools filter -O v -T $LOH/PEI_hetero_LOH_10_counts.bed $WORKD/Somatypus_SNVs_final.vcf.gz > $OUTPUT/Somatypus_SNVs_PEI_LOH.vcf &
wait


#awk '$0 !~ /^##contig=.*/ && $0 !~ /^##FILTER=<ID=PASS.*/ && $0 !~ /^##bcftools_filter.*/  {print}' Somatypus_SNVs_USA_noLOH.vcf > SNVs_USA_noLOH.vcf &
awk '$0 !~ /^##contig=.*/ && $0 !~ /^##FILTER=<ID=PASS.*/ && $0 !~ /^##bcftools_filter.*/  {print}' Somatypus_SNVs_USA_LOH.vcf > SNVs_USA_LOH.vcf &
#awk '$0 !~ /^##contig=.*/ && $0 !~ /^##FILTER=<ID=PASS.*/ && $0 !~ /^##bcftools_filter.*/  {print}' Somatypus_SNVs_PEI_noLOH.vcf > SNVs_PEI_noLOH.vcf &
awk '$0 !~ /^##contig=.*/ && $0 !~ /^##FILTER=<ID=PASS.*/ && $0 !~ /^##bcftools_filter.*/  {print}' Somatypus_SNVs_PEI_LOH.vcf > SNVs_PEI_LOH.vcf &
wait

#python2 $COMMAND/ExtractVcfData.py SNVs_USA_noLOH.vcf &
python2 $COMMAND/ExtractVcfData.py SNVs_USA_LOH.vcf &
#python2 $COMMAND/ExtractVcfData.py SNVs_PEI_noLOH.vcf &
python2 $COMMAND/ExtractVcfData.py SNVs_PEI_LOH.vcf &
wait

rm *.vcf

