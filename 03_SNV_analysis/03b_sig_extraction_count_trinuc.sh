# original file here: C:\Users\shart\Metzger Lab Dropbox\Sam_data\Mya_genome\SNPs\somatypus_output_sigS_helmsman.sh

## helmsman run: done outside R
SNPS=/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run
OUTPUT=/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/helmsman
GENOME=/ssd3/Mar_genome_analysis/genomes/Mar.3.4.6.p1_Q30Q30A.fasta

cd $SNPS
module load conda/
source activate helmsman
python /ssd3/Mar_genome_analysis/mut_sig/helmsman/helmsman.py --verbose --decomp nmf --mode txt --input $SNPS/helmsman.txt --fastafile $GENOME --projectdir $OUTPUT --length 3 
#python /ssd3/Mar_genome_analysis/mut_sig/helmsman/helmsman.py --verbose --decomp nmf --mode txt --input $SNPS/helmsman.txt --fastafile $GENOME --projectdir $OUTPUT/length5 --length 5
#python /ssd3/Mar_genome_analysis/mut_sig/helmsman/helmsman.py --verbose --decomp nmf --mode txt --input $SNPS/helmsman.txt --fastafile $GENOME --projectdir $OUTPUT/length7 --length 7
