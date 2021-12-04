# originally here: C:\Users\shart\Metzger Lab Dropbox\Sam_data\Mya_genome\repeats\generate_repeat_library.sh

# Merge libaries from repeatmodeler (reference genome), repdenovo (1xHEALTHY, 1xPEI, 1xUSA)
    RM_LIB=/ssd3/Mar_genome_analysis/repeatmodeler/Mar.3.4.6.p1_Q30Q30A-families.fa
    RD_LIB_USA=/ssd2/REPdenovo/MELC_A11_output/contigs.fa.classified
    RD_LIB_PEI=/ssd2/REPdenovo/PEI-DN08_output/contigs.fa.classified
    RD_LIB_REF=/ssd2/REPdenovo/MELC-2E11_output/contigs.fa.classified
    COMMAND=/ssd3/Mar_genome_analysis/repeat_lib/cd-hit-v4.8.1-2019-0228
    WORKD=/ssd3/Mar_genome_analysis/repeat_lib

# Install CD-HIT - program used by repeatmodeler2 to merge libraries
    # cd /ssd3/Mar_genome_analysis/repeat_lib
    # wget https://github.com/weizhongli/cdhit/releases/download/V4.8.1/cd-hit-v4.8.1-2019-0228.tar.gz
    # tar xvf cd-hit-v4.8.1-2019-0228.tar.gz --gunzip
    # cd cd-hit-v4.8.1-2019-0228
    # make
    # cd cd-hit-auxtools
    # make
    # # TEST run
    # cat $RD_LIB_USA $RD_LIB_PEI $RD_LIB_REF > $WORKD/repdenovo_UPR_cat
    # $COMMAND/cd-hit-est -aS 0.8 -c 0.8 -g 1 -G 0 -A 80 -M 10000 -i $WORKD/repdenovo_UPR_cat -o $WORKD/repdenovo_UPR_combined -T 20

# Rename input files
    cd $WORKD
    cp $RM_LIB ref.genome.repeatmoduler
    cp $RD_LIB_USA usa.repdenovo
    cp $RD_LIB_PEI pei.repdenovo
    cp $RD_LIB_REF ref.repdenovo
    awk '$0 !~ /^>.*/ {print} $0 ~ /^>.*/ {split($1, name, "#" ); print name[1] "#usaRD#" name[2]}' usa.repdenovo > usa.repdenovo.named
    awk '$0 !~ /^>.*/ {print} $0 ~ /^>.*/ {split($1, name, "#" ); print name[1] "#peiRD#" name[2]}' pei.repdenovo > pei.repdenovo.named
    awk '$0 !~ /^>.*/ {print} $0 ~ /^>.*/ {split($1, name, "#" ); print name[1] "#refRD#" name[2]}' ref.repdenovo > ref.repdenovo.named
    awk '$0 !~ /^>.*/ {print} $0 ~ /^>.*/ {split($1, name, "#" ); print name[1] "#refRM#" name[2]}' ref.genome.repeatmoduler > ref.genome.repeatmoduler.named
    cat usa.repdenovo.named pei.repdenovo.named ref.repdenovo.named ref.genome.repeatmoduler.named > all_repeat_lib_cat

# Run CD-HIT and count number of elements in input and output
    $COMMAND/cd-hit-est -aS 0.8 -c 0.8 -g 1 -G 0 -A 80 -M 10000 -i $WORKD/all_repeat_lib_cat -o $WORKD/all_repeat_lib.cdhit -T 50
    bwa index all_repeat_lib.cdhit
    FILES="usa.repdenovo.named pei.repdenovo.named ref.repdenovo.named ref.genome.repeatmoduler.named all_repeat_lib.cdhit"
    for i in $FILES
    do
        awk '$0 ~ /^>.*/ {count += 1} END{print count}' $i
    done

# bwa map reads to library, keep only mapping reads
    LIST="MELC-A11_S1 PEI-DN08_S3 NYTC-C9_S2 MELC-2E11 MELC-A9 PEI-DF490 DF-488 DN-HL03 FFM-19G1 FFM-20B2"
    ILLUMINA=/ssd2/Illumina_data/dedupe/trim20
    for i in $LIST
    do
        bwa mem -t 10 -v 2 all_repeat_lib.cdhit $ILLUMINA/$i"_R1_001.fastq.gz" $ILLUMINA/$i"_R2_001.fastq.gz" | samtools view -b -h -F 4 | samtools sort -O bam > $WORKD/repeatLIB.$i.bam &
    done
    wait
    for i in $LIST
    do
        samtools index -b -@ 10 $WORKD/repeatLIB.$i.bam &
    done
    wait
# FFM-22F10 mapped seperately since it's from newer sequencing run
    bwa mem -t 70 -v 2 all_repeat_lib.cdhit /ssd2/Illumina_data/dedupe/2021_newsamples_trim20/FFM-22F10_R1_001.fastq.gz /ssd2/Illumina_data/dedupe/2021_newsamples_trim20/FFM-22F10_R2_001.fastq.gz | \
        samtools view -b -h -F 4 -@ 10 | samtools sort -@ 10 -O bam > $WORKD/repeatLIB.FFM-22F10.bam
    samtools index -b -@ 10 $WORKD/repeatLIB.FFM-22F10.bam &


# get total nucleotides mapped to each 
    #samtools depth (would then need to add up for every repeat)
    #samtools coverage (computes summary stats on its own)
    cd $WORKD
    LIST="MELC-A11_S1 PEI-DN08_S3 NYTC-C9_S2 MELC-2E11 MELC-A9 PEI-DF490 DF-488 DN-HL03 FFM-19G1 FFM-20B2 FFM-22F10"
    for i in $LIST
    do
        # samtools coverage repeatLIB.$i.bam -o repeatLIB.$i.coverage &
        samtools depth -a -a -d 0 repeatLIB.$i.bam > repeatLIB.$i.coverage &
    done
    wait

# Normalize by repeat length and average coverage
    # average coverage calculated based on average depth at all locations where an SNV is called in ANY sample during somatypus run:
        # C:\Users\shart\Metzger Lab Dropbox\Sam_data\Mya_genome\SNPs\somatypus_output_pairwise.r
    for i in $LIST
    do
    awk -v sample=$i -v FS='\t' -v OFS='\t' '
        BEGIN{print "repeat_name", "length", sample;
            name = "NEW_CONTIG_MERGE_2#usaRD#Unknown";
            if(sample ~ "FFM-19G1" ){ avgcov = 92.81066};
            if(sample ~ "FFM-20B2" ){ avgcov = 92.78821};
            if(sample ~ "FFM-22F10" ){ avgcov = 63.84493};
            if(sample ~ "MELC-A11_S1" ){ avgcov = 63.49766};
            if(sample ~ "NYTC-C9_S2" ){ avgcov = 59.25352};
            if(sample ~ "DF-488" ){ avgcov = 94.64185};
            if(sample ~ "DN-HL03" ){ avgcov = 89.04225};
            if(sample ~ "PEI-DN08_S3" ){ avgcov = 65.3303};
            if(sample ~ "MELC-2E11" ){ avgcov = 59.91782};
            if(sample ~ "MELC-A9" ){ avgcov = 60.0784};
            if(sample ~ "PEI-DF490" ){ avgcov = 60.3944}
        } 
        name ~ $1 {rep_length += 1; coverage += $3};
        name !~ $1 {print name, rep_length, coverage/rep_length/avgcov; name = $1; rep_length = 1; coverage = $3}
        END{print name, rep_length, coverage/rep_length/avgcov}' OFS='\t' repeatLIB.$i.coverage > repeatLIB.$i.depth &
    done

# delete intermidiate bam and .coverage files to save space
    # for i in $LIST
    # do
    #     rm repeatLIB.$i.bam
    #     rm repeatLIB.$i.bam.bai
    #     rm repeatLIB.$i.coverage
    # done
