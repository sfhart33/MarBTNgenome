# From here: /ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/bins/dnds/README

module load R/3.6.0

# DETERMINING LOH THRESHOLD

DNDS=/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/bins/dnds/dndscv_script.R
INPUT=/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/bins/LOH
OUPUT=/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/bins/dnds/LOH

cd $INPUT
FILES=$(ls *.dnds)
FILECOUNT=$(ls *.dnds | wc -l)
FILESUBSET1="allUSAnoPEInoH_stringent"
FILESUBSET2="allUSAnoPEInoH_regular"
FILESUBSET3="allPEInoUSAnoH_stringent"
FILESUBSET4="allPEInoUSAnoH_regular"
COUNT=0
for subset in $FILESUBSET1 $FILESUBSET2 $FILESUBSET3 $FILESUBSET4
do
	cd $INPUT
	ls $subset* | wc -l
	FILES=$(ls $subset*)
	echo $FILES
	cd $OUPUT
	for file in $FILES
	do
	{ Rscript $DNDS $file $INPUT; let "COUNT+=1"; echo -e $COUNT/$FILECOUNT "\t" $file finished; } & 
	done
	wait
	echo "#################################### SET DONE ##################################"
done

# RUN FOR ALL BINS

DNDS=/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/bins/dnds/dndscv_script.R
INPUT=/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/bins/dnds/outputs
OUPUT=/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/bins/dnds/outputs

cd $INPUT
FILES=$(ls *.dnds)
FILECOUNT=$(ls *.dnds | wc -l)
echo $FILES
echo $FILECOUNT
cd $OUPUT
for file in $FILES
do
	{ Rscript $DNDS $file $INPUT; echo -e "#######################" $file finished "######################"; } & 
done
wait

# RUN FOR TWO POSSIBLE LOH THRESHOLDS FOR USA AND PEI SUBLINEAGES

# make somatic bins
SOMTIC=/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/bins/dnds/somatic
LOH_INPUT=/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/bins/LOH
INPUT=/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/bins/dnds/outputs
cd $INPUT
UNIQUE=$(ls uniq*dnds)
COUNTS="10 25"
STRINGENCY="stringent regular"
for count in $COUNTS
do
for stringency in $STRINGENCY
do
cat $LOH_INPUT/"allPEInoUSAnoH_"$stringency"_"$count"_nonLOH.dnds" $LOH_INPUT/"allUSAnoPEInoH_"$stringency"_"$count"_nonLOH.dnds" > $SOMTIC/"sublineages_"$stringency"_"$count"_nonLOH_alone.dnds"
cat $UNIQUE $SOMTIC/"sublineages_"$stringency"_"$count"_nonLOH.dnds" > $SOMTIC/"sublineages_"$stringency"_"$count"_nonLOH_uniq.dnds"
cat $SOMTIC/"sublineages_"$stringency"_"$count"_nonLOH_uniq.dnds" multipleU1234noPEInoH.dnds multipleU12noPEInoH.dnds > $SOMTIC/"sublineages_"$stringency"_"$count"_nonLOH_uniqmult.dnds"
done
done


DNDS=/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/bins/dnds/dndscv_script.R
INPUT=/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/bins/dnds/somatic
OUPUT=/ssd3/Mar_genome_analysis/somatypus/Mar.3.4.6.p1/final_run/bins/dnds/somatic

cd $INPUT
FILES=$(ls *.dnds)
FILECOUNT=$(ls *.dnds | wc -l)
echo $FILES
echo $FILECOUNT
cd $OUPUT
for file in $FILES
do
	{ Rscript $DNDS $file $INPUT; echo -e "#######################" $file finished "######################"; } & 
done
wait
