SAMPLES="MELC-2E11 MELC-A9 PEI-DF490 FFM-19G1 FFM-20B2 FFM-22F10 MELC-A11_S1 NYTC-C9_S2 DF-488 DN-HL03 PEI-DN08_S3"
WORKD=/ssd3/Mar_genome_analysis/steamer/final_pipeline/logos
cd $WORKD

echo -e "sample\tbp\tnuc\tcount\ttotal" > steamer_bias_down.txt
echo -e "sample\tbp\tnuc\tcount\ttotal" > steamer_bias_up.txt

for sample in $SAMPLES
do
    for i in {1..20}
    do
        for n in A T G C
        do
            awk -v sample=$sample -v count=$i -v nuc=$n ' BEGIN{Ncount=0}
            $0 ~ /^@/ {
                getline;
                if (length($1)>20){
                    NRcount +=1;
                    if ( substr($1, count, 1) ~ nuc){
                        Ncount+=1
                    }
                }
            }
            END{print sample, count, nuc, Ncount, NRcount}' OFS='\t' /ssd3/Mar_genome_analysis/steamer/final_pipeline/$sample/softclips_down_flanking.fastq >> steamer_bias_down.txt
            # $DOWNSTREAM
        done
    done
done

for sample in $SAMPLES
do
    for i in {1..20}
    do
        for n in A T G C
        do
            awk -v sample=$sample -v count=$i -v nuc=$n ' BEGIN{Ncount=0}
            $0 ~ /^@/ {
                getline;
                if (length($1)>20){
                    NRcount +=1;
                    if ( substr($1, (length($1)+1-count), 1) ~ nuc){
                        Ncount+=1
                    }
                }
            }
            END{print sample, count, nuc, Ncount, NRcount}' OFS='\t' /ssd3/Mar_genome_analysis/steamer/final_pipeline/$sample/softclips_up_flanking.fastq >> steamer_bias_up.txt
            # $DOWNSTREAM
        done
    done
done

# Get site biases directly from genome


FILES="all_sites anyCnoH allCnoH_allPEI_allUSA allCnoH allUSA allPEI"
INPUT=/ssd3/Mar_genome_analysis/steamer/final_pipeline/merge
GENOME=/ssd3/Mar_genome_analysis/genomes/Mar.3.4.6.p1_Q30Q30A.fasta
cd $WORKD
for bed in $FILES
do
    awk '{start = $2-15; end = $3+15; print $1, start, end, $4, $5, $6}'  OFS='\t' $INPUT/$bed.bed | \
        bedtools getfasta -bedOut -s -fi $GENOME -bed - > $WORKD/$bed.seq.bed
done

# create R input 

echo -e "sample\tbp\tnuc\tcount\ttotal" > steamer_bias_genome.txt

for sample in $FILES
do
    for i in {1..35}
    do
        for n in A T G C
        do
            awk -v sample=$sample -v count=$i -v nuc=$n ' BEGIN{Ncount=0}
            {if ( substr($7, count, 1) ~ nuc){
                    Ncount+=1
                }
            }
            END{print sample, count, nuc, Ncount, NR}' OFS='\t' $WORKD/$sample.seq.bed >> steamer_bias_genome.txt
        done
    done
done

# count number with ATG in positions 7/8/9
    wc -l allCnoH.seq.bed
    awk '$7 ~ "..................................." {count +=1} END{print count}' anyCnoH.seq.bed # total: 550
    awk '$7 ~ "...........CAT....................." {count +=1} END{print count}' anyCnoH.seq.bed # CAT upstream: 103
    awk '$7 ~ ".....................ATG..........." {count +=1} END{print count}' anyCnoH.seq.bed # ATG downstream: 90
    awk '$7 ~ "...........CAT.......ATG..........." {count +=1} END{print count}' anyCnoH.seq.bed # CAT up and ATG down: 12
# Expected calculated in R
