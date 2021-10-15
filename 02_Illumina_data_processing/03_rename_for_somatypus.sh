cd /ssd3/Mar_genome_analysis/bwa_mapping/Mar.3.4.6.p1/all_samples

# rename PEI samples to match naming convension of other samples
mv Mar.3.4.6.p1.DN-HL07.bam Mar.3.4.6.p1.PEI-DN07.bam
mv Mar.3.4.6.p1.DN-HL03.bam Mar.3.4.6.p1.PEI-DN03.bam
mv Mar.3.4.6.p1.DF-488.bam Mar.3.4.6.p1.PEI-DF488.bam

# rename so that somatypus will run okay

LIST="MELC-2E11 MELC-A9 PEI-DF490 PEI-DF488 PEI-DN03 PEI-DN07 PEI-DN08_S3 FFM-19G1 FFM-20B2 FFM-22A10 FFM-22F10 MELC-A10 MELC-A11_S1 NYTC-C9_S2 PEI-DN03-siphon PEI-DN07-siphon PEI-DN08-siphon FFM-22A10-adductor FFM-22F10-adductor MELC-A10-siphon MELC-A11-siphon NYTC-C9-mantle"

# test to make sure I have all correct names
for i in $LIST
do
	stat --printf="%s" Mar.3.4.6.p1.$i.bams
done
# test renaming loop
count=1
for i in $LIST
do
	echo $count	
	printf -v counttwodigits  "%02d" $count
	echo $counttwodigits.$i.bam
	count=$(expr $count + 1)
	
done


count=1
rm *.bam.bai
for i in $LIST
do
	echo $count	
	printf -v counttwodigits  "%02d" $count
	mv Mar.3.4.6.p1.$i.bam $counttwodigits.$i.bam
	count=$(expr $count + 1)
	
done

count=1
for i in $LIST
do
	echo $count	
	printf -v counttwodigits  "%02d" $count
	samtools index -b -@ 5 $counttwodigits.$i.bam &
	count=$(expr $count + 1)
done
wait 


