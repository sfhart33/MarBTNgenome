
module load repeatmodeler
GENOMES=/ssd3/Mar_genome_analysis/genomes
WORKD=/ssd3/Mar_genome_analysis/repeatmodeler
i=Mar.3.4.6.p1_Q30Q30A
	cd $GENOMES
	BuildDatabase -name $i $i.fasta
	cd $WORKD
	# RepeatModeler -database $GENOMES/$i -pa 50 -LTRStruct &> $i.RM.logfile
	# RepeatModeler -database $GENOMES/$i -pa 30 -LTRStruct &> $i.RM.rerun.logfile
	############ Failed^ maybe used too many threads ############
	RepeatModeler -database $GENOMES/$i -pa 20 -LTRStruct &> $i.RM.rerun2.logfile
	awk ' $1 ~ /^>/ {count += 1} END {print count}' $i.consensi.fa.classified > $i.RM.elementcount
	# Mar.3.4.6.p1_Q30Q30A-families.fa
	awk ' $1 ~ /^>/ {count += 1} END {print count}' Mar.3.4.6.p1_Q30Q30A-families.fa > Mar.3.4.6.p1_Q30Q30A-families.elementcount
	RepeatMasker -pa 20 -lib $i-families.fa $GENOMES/$i.fasta &> $i.RMasker.logfile


