    SNPS=/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/pairwise/noLOH
    OUTPUT=/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/pairwise/noLOH/helmsman_div_output
    OUTPUT2=/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/pairwise/noLOH/helmsman_all_output
    GENOME=/ssd3/Mar_genome_analysis/genomes/Mar.3.4.6.p1_Q30Q30A.fasta

    cd $SNPS
    module load conda/
    source activate helmsman
    python /ssd3/Mar_genome_analysis/mut_sig/helmsman/helmsman.py --verbose --decomp nmf --mode txt --input $SNPS/helmsman.divergence.txt --fastafile $GENOME --projectdir $OUTPUT --length 3 
    python /ssd3/Mar_genome_analysis/mut_sig/helmsman/helmsman.py --verbose --decomp nmf --mode txt --input $SNPS/helmsman.output.txt --fastafile $GENOME --projectdir $OUTPUT2 --length 3 

    SNPS=/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/pairwise/withLOH
    OUTPUT=/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/pairwise/withLOH/helmsman_div_output
    # OUTPUT2=/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/pairwise/withLOH/helmsman_all_output
    GENOME=/ssd3/Mar_genome_analysis/genomes/Mar.3.4.6.p1_Q30Q30A.fasta

    cd $SNPS
    python /ssd3/Mar_genome_analysis/mut_sig/helmsman/helmsman.py --verbose --decomp nmf --mode txt --input $SNPS/helmsman.divergence.txt --fastafile $GENOME --projectdir $OUTPUT --length 3 

