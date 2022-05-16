# From here: /ssd3/Mar_genome_analysis/genomes/stats/README

# Two options for calculating N50 from Michael
	# module load conda2
	# conda activate pbbioconda-0.0.5
	# falconc stats-assembly -f ~/MELC-2E11/4-polish/cns-output/cns_p_ctg.fasta
	# 
	# conda activate denovo_asm
	# python scripts/get_asm_stats.py 4-polish/cns-output/cns_p_ctg.fasta

INPUT=/ssd3/Mar_genome_analysis/genomes
OUTPUT=/ssd3/Mar_genome_analysis/genomes/stats
GENOMES="Mya.genome.v1.01 Mar.1.1.1 Mar.3.1.1 Mar.3.2.3_curated.FALC Mar.3.3.2_p1_PGA_assembly Mar.3.3.3.p1 Mar.3.4.6.p1_Q30Q30A"
GENOME=/ssd3/Mar_genome_analysis/genomes/Mar.3.4.6.p1_Q30Q30A.fasta

# LIST GENOMES
for i in $GENOMES
do
	echo $i
	head -1 $INPUT/$i.fasta
	wc -l $INPUT/$i.fasta
done

# CALCULATE SIZE, GC, N50, CONTIGS
for i in $GENOMES
do
	echo $i
	fold -w 50 $INPUT/$i.fasta > $OUTPUT/$i.wrap.fasta
	/ssd2/Illumina_data/bbmap/stats.sh $OUTPUT/$i.wrap.fasta > $OUTPUT/$i.bbmap.stats
done

# RUN REPEATMASKER
module load repeatmodeler
FULL_LB=/ssd3/Mar_genome_analysis/repeat_lib/all_repeat_lib.cdhit
RM_LB=/ssd3/Mar_genome_analysis/repeatmodeler/Mar.3.4.6.p1_Q30Q30A-families.fa
# GENOMES="Mar.3.1.1 Mar.3.2.3_curated.FALC Mar.3.3.2_p1_PGA_assembly Mar.3.3.3.p1 Mar.3.4.6.p1_Q30Q30A"
for i in $GENOMES
do 
	echo $i
	# RepeatMasker -pa 20 -lib $FULL_LB -dir $OUTPUT/FULL_LB $INPUT/$i.fasta
	RepeatMasker -pa 10 -lib $RM_LB -dir $OUTPUT/RM_LB $INPUT/$i.fasta &
done
# DELETE OUTPUTS BESIDES .tbl

# to do:



# RUN BUSCO

INPUT=/ssd3/Mar_genome_analysis/genomes
OUTPUT=/ssd3/Mar_genome_analysis/genomes/stats/busco
# GENOMES="Mya.genome.v1.01 Mar.1.1.1 Mar.3.1.1"
BUSCO=/biotools/busco/scripts
DATA=/home/shart/busco
export AUGUSTUS_CONFIG_PATH=/ssd2/augustus_busco/Mar_BUSCO/myAugconfig

#for i in $GENOMES
for i in $GENOMES
do
	python3.5 $BUSCO/run_BUSCO.py -i $INPUT/$i.fasta -o $i.metazoa_odb10.busco -l $DATA/metazoa_odb10 -m geno -f -z -c 20 &
done

# MOVE FROM /var/tmp/busco
# cp -R /var/tmp/busco/*.metazoa_odb10.busco $OUTPUT


	# python3 $BUSCO/run_BUSCO.py -i $GENOMES/$i -o $i.mollusca_odb10.busco -l $DATA/mollusca_odb10 -m geno -f -z &> $i.mollusca.logfile &
	# MOLLUSCA RUNS ARE TAKING AWHILE, PERHAPS GOT STUCK SOMEWHERE?