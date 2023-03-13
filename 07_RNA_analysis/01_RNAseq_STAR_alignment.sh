# notes here: C:\Users\shart\Metzger Lab Dropbox\Sam_data\Mya_genome\RNAseq\RNAseq_STAR_alignment

RNADIR=/ssd3/RNAseq

cd $RNADIR
# make own version of GFF>GTF with correct names
        cd /ssd3/Mar_genome_analysis/genomes/maker/
        awk ' $3 ~ "exon" || $3 ~ "CDS" {split($9, names, ";" );
                                        split(names[2], id, "=" );
                                        genecol = "gene_id " "\x22" id[2] "\x22" "; transcript_id " "\x22" id[2] "\x22" ";";
                                        print $1, $2, $3, $4, $5, $6, $7, $8, genecol
                                        } ' OFS="\t" 2020-09-11_Mar_genome_snap02.all.noseq.gff > 2020-09-11_Mar_genome_snap02.all.manual.gtf

        head 2020-09-11_Mar_genome_snap02.all.noseq.gff 
        head 2020-09-11_Mar_genome_snap02.all.manual.gtf

module load star
# https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf
# GENOME INDEX
STAR \
        --runThreadN 50 \
        --runMode genomeGenerate \
        --genomeDir /ssd3/RNAseq/STAR_genome \
        --readFilesCommand zcat \
        --genomeFastaFiles /ssd3/Mar_genome_analysis/genomes/Mar.3.4.6.p1_Q30Q30A.fasta \
        --sjdbGTFfile /ssd3/Mar_genome_analysis/genomes/maker/2020-09-11_Mar_genome_snap02.all.manual.gtf

# RUN STAR
        BTN_RNA=/ssd3/Mar_genome_analysis/RNAseq
        BTN_RNA_SAMPLES="SH-1 SH-2 SH-3 SH-4"
        TISSUE_RNA=/ssd2/Trinity/All_samples_seperate
        TISSUE_RNA_SAMPLES="1-mantle 2-foot 3-siphon 5-muscle 6-gills 7-hemocytes"

        cd /ssd3/RNAseq
        for sample in $BTN_RNA_SAMPLES
        do  
            STAR \
                --runThreadN 10 \
                --genomeDir /ssd3/RNAseq/STAR_genome \
                --readFilesIn $BTN_RNA/$sample"_R1_001.fastq.gz" $BTN_RNA/$sample"_R2_001.fastq.gz" \
                --outFileNamePrefix  /ssd3/RNAseq/STAR_genome_output/$sample"_" \
                --readFilesCommand zcat \
                --quantMode GeneCounts \
                --outSAMtype None &
        done

        for sample in $TISSUE_RNA_SAMPLES
        do  
            STAR \
                --runThreadN 10 \
                --genomeDir /ssd3/RNAseq/STAR_genome \
                --readFilesIn $TISSUE_RNA/"MELC-2E11-"$sample"_R1_001.fastq.gz" $TISSUE_RNA/"MELC-2E11-"$sample"_R2_001.fastq.gz" \
                --outFileNamePrefix  /ssd3/RNAseq/STAR_genome_output/$sample"_" \
                --readFilesCommand zcat \
                --quantMode GeneCounts \
                --outSAMtype None &
        done
        wait

        NEW_RNASEQ=/ssd3/RNAseq/seq
        cd $NEW_RNASEQ
        NEW_SAMPLES=$(ls *R1_001.fastq.gz | awk '{split($0,name,"_R1"); print name[1]}')
        cd /ssd3/RNAseq
        for sample in $NEW_SAMPLES
        do  
            STAR \
                --runThreadN 5 \
                --genomeDir /ssd3/RNAseq/STAR_genome \
                --readFilesIn $NEW_RNASEQ/$sample"_R1_001.fastq.gz" $NEW_RNASEQ/$sample"_R2_001.fastq.gz" \
                --outFileNamePrefix  /ssd3/RNAseq/STAR_genome_output/$sample"_" \
                --readFilesCommand zcat \
                --quantMode GeneCounts \
                --outSAMtype None &
        done


ls *ReadsPerGene.out.tab | awk '{split($0,name,"_ReadsPerGene"); print name[1]}' > samples_list.txt

