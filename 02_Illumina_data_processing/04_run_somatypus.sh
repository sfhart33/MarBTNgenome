# MOVE ALL FILES INTO THIS FOLDER FIRST:
BAM_FILES=/ssd3/Mar_genome_analysis/bwa_mapping/Mar.3.4.6.p1/all_samples
WORK_DIR=/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run
GENOME=/ssd3/Mar_genome_analysis/genomes/Mar.3.4.6.p1_Q30Q30A.fasta

cd $WORK_DIR
module load somatypus/
somatypus -i $BAM_FILES -o $WORK_DIR/MAPQ30_new -g $GENOME -c 8

# started 6/29/21 at 4:41pm
# Lost connection 7/4/21 at 7:40pm: resume run using checkpoint file by re-running same command, finished shortly after
# deleted intermidiate files since they took up a to of space

# try the tools that comes with it for extracting info

# test here: 

COMMAND=/ssd3/Mar_genome_analysis/somatypus
cd $COMMAND
./ExtractVcfData.py $WORK_DIR/Somatypus_SNVs_final.vcf
./ExtractVcfData.py $WORK_DIR/Somatypus_Indels_final.vcf
cd $WORK_DIR

# divide VCF  regions
GENES=/ssd3/Mar_genome_analysis/genomes/maker
TYPES="gene CDS exon five_prime_UTR three_prime_UTR"

# compress vcf to vcf.gz for bcftools
bgzip -c $WORK_DIR/Somatypus_SNVs_final.vcf > $WORK_DIR/Somatypus_SNVs_final.vcf.gz
tabix -fp vcf $WORK_DIR/Somatypus_SNVs_final.vcf.gz

########################
for type in $TYPES
do
bcftools filter -R $GENES/2020-09-11_Mar_genome_snap02.all.$type.bed $WORK_DIR/Somatypus_SNVs_final.vcf.gz > $WORK_DIR/Somatypus_SNVs_final.$type.vcf &
done
rm $WORK_DIR/Somatypus_SNVs_final.vcf.gz
rm $WORK_DIR/Somatypus_SNVs_final.vcf.gz.tbi

# reformat so commented lines match that expected by ExtractVcfData.py
for type in $TYPES
do
awk '$0 !~ /^##contig=.*/ && $0 !~ /^##FILTER=<ID=PASS.*/  {print}' $WORK_DIR/Somatypus_SNVs_final.$type.vcf > $WORK_DIR/Somatypus_SNVs_final.$type.vcf2 &
done
wait
for type in $TYPES
do
mv $WORK_DIR/Somatypus_SNVs_final.$type.vcf2 $WORK_DIR/Somatypus_SNVs_final.$type.vcf
done


COMMAND=/ssd3/Mar_genome_analysis/somatypus
cd $WORK_DIR
for type in $TYPES
do
$COMMAND/ExtractVcfData.py Somatypus_SNVs_final.$type.vcf &
done

for type in $TYPES
do
rm $WORK_DIR/Somatypus_SNVs_final.$type.vcf
done

