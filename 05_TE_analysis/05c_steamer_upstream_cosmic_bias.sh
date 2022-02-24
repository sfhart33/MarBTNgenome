# cosmic genes associated with cancer
# Manually downloaded cancer gene census and all cosmic genes from here: https://cancer.sanger.ac.uk/cosmic/download
	cd /ssd3/Mar_genome_analysis/cosmic
	gunzip All_COSMIC_Genes.fasta.gz
	wc -l All_COSMIC_Genes.fasta
	wc -l Census_all_09-12-2021_2021.tsv
	head Census_all_09-12-2021_2021.tsv
	grep ACVR1 All_COSMIC_Genes.fasta
	grep "ACVR2A " All_COSMIC_Genes.fasta

	awk 'NR > 1 {print $1 }' Census_all_09-12-2021_2021.tsv > cosmic_cancer_gene_census.txt
	awk 'NR > 1 {print $1 }' Census_all_09-12-2021_2021.tsv | head -10 > cosmic_cancer_gene_census_test.txt
	seqtk subseq All_COSMIC_Genes.fasta cosmic_cancer_gene_census_test.txt > cosmic_cancer_gene_census_test.fasta
	seqtk subseq All_COSMIC_Genes.fasta cosmic_cancer_gene_census.txt > cosmic_cancer_gene_census.fasta

module load blast+/2.10.0

# MICHAEL:
    # # Use Mya arenaria proteins to search uniprot to find likely gene identification
    #     makeblastdb -in uniprot_sprot.fasta -out uniprot_sprot -dbtype prot
    #     blastp -query /home/metzgerm/MAKER_Mya/2020-09-11-Mar.3.4.6.p1-MAKER/snap02/2020-09-11_Mar_genome_snap02.all.maker.proteins.fasta -db uniprot_sprot -evalue 1e-6 -max_hsps 1 -max_target_seqs 1 -outfmt 6 -out 2020-09-11_Mar_genome_snap02.all.maker.proteins.fasta.blastp -num_threads 20

    # # Use Mya arenaria proteins to search a concatenated file of well-annotated bivalves to find likely gene identification
    #     makeblastdb -in CgiCviMcoPmaMye_protein.fasta -out CgiCviMcoPmaMye_protein -dbtype prot
    #     blastp -query /home/metzgerm/MAKER_Mya/2020-09-11-Mar.3.4.6.p1-MAKER/snap02/2020-09-11_Mar_genome_snap02.all.maker.proteins.fasta -db CgiCviMcoPmaMye_protein -evalue 1e-6 -max_hsps 1 -max_target_seqs 1 -outfmt 6 -out 2020-09-11_Mar_genome_snap02.all.maker.proteins.fasta.CgiCviMcoPmaMye_blastp -num_threads 20

    # # Use Mya arenaria proteins to search COSMIC proteins to find which Mya arenaria proteins have COSMIC orthologs
    #     makeblastdb -in All_COSMIC_Proteins.fasta -out All_COSMIC_Proteins -dbtype prot
    #     blastp -query /home/metzgerm/MAKER_Mya/2020-09-11-Mar.3.4.6.p1-MAKER/snap02/2020-09-11_Mar_genome_snap02.all.maker.proteins.fasta -db All_COSMIC_Proteins -evalue 1e-6 -max_hsps 1 -max_target_seqs 1 -outfmt 6 -out 2020-09-11_Mar_genome_snap02.all.maker.proteins.fasta.All_COSMIC_Proteins -num_threads 20

    # # Use COSMIC genes and proteins to search Mya arenaria proteins to find how many COSMIC genes/proteins have Mya arenaria orthologs
    #     makeblastdb -in /home/metzgerm/MAKER_Mya/2020-09-11-Mar.3.4.6.p1-MAKER/snap02/2020-09-11_Mar_genome_snap02.all.maker.proteins.fasta -out 2020-09-11_Mar_genome_snap02.all.maker.proteins -dbtype prot
    #     blastx -query All_COSMIC_Genes_dMT.fasta -db 2020-09-11_Mar_genome_snap02.all.maker.proteins -evalue 1e-6 -max_hsps 1 -max_target_seqs 1 -outfmt 6 -out All_COSMIC_Genes_dMT.fasta.2020-09-11_Mar_genome_snap02.all.maker.proteins -num_threads 20

    #     blastp -query All_COSMIC_Proteins_dMT.fasta -db 2020-09-11_Mar_genome_snap02.all.maker.proteins -evalue 1e-6 -max_hsps 1 -max_target_seqs 1 -outfmt 6 -out All_COSMIC_Proteins_dMT.fasta.2020-09-11_Mar_genome_snap02.all.maker.proteins -num_threads 20

    # # 42360/56461 = 75.0%

# ME:
    MYAGENES=/ssd3/Mar_genome_analysis/genomes/maker
    makeblastdb -in $MYAGENES/2020-09-11_Mar_genome_snap02.all.maker.proteins.fasta -out $MYAGENES/2020-09-11_Mar_genome_snap02.all.maker.proteins -dbtype prot
    makeblastdb -in cosmic_cancer_gene_census.fasta -out cosmic_cancer_gene_census -dbtype nucl
    makeblastdb -in cosmic_cancer_gene_census_test.fasta -out cosmic_cancer_gene_census_test -dbtype nucl


    blastx -query cosmic_cancer_gene_census_test.fasta -db $MYAGENES/2020-09-11_Mar_genome_snap02.all.maker.proteins -evalue 1e-6 -max_hsps 1 -max_target_seqs 1 -outfmt 6 -out test1 -num_threads 20
    tblastn -query $MYAGENES/2020-09-11_Mar_genome_snap02.all.maker.proteins.fasta -db cosmic_cancer_gene_census_test -evalue 1e-6 -max_hsps 1 -max_target_seqs 1 -outfmt 6 -out test2 -num_threads 20

    blastx -query cosmic_cancer_gene_census.fasta -db $MYAGENES/2020-09-11_Mar_genome_snap02.all.maker.proteins -evalue 1e-6 -max_hsps 1 -max_target_seqs 1 -outfmt 6 -out test3 -num_threads 20
        wc -l test3 # 633
        wc -l cosmic_cancer_gene_census.txt # 729
        # 633/729 = 86%
    tblastn -query $MYAGENES/2020-09-11_Mar_genome_snap02.all.maker.proteins.fasta -db cosmic_cancer_gene_census -evalue 1e-6 -max_hsps 1 -max_target_seqs 1 -outfmt 6 -out test4c -num_threads 10
        # only runs with threads of 10 - failed at 20, 50 and 100
        wc -l test4c # 5450
        grep ">" $MYAGENES/2020-09-11_Mar_genome_snap02.all.maker.proteins.fasta | wc -l # 38609
        # 5450/38609 = 14%

    tblastn -query $MYAGENES/2020-09-11_Mar_genome_snap02.all.maker.proteins.fasta -db cosmic_cancer_gene_census -evalue 1e-10 -max_hsps 1 -max_target_seqs 1 -outfmt 6 -out test5 -num_threads 5
    wc -l test5
    mv test5 mya_genes_cosmic_e10_4164.blast
    mv test4c mya_genes_cosmic_e6_5450.blast
    mv test3 top_mya_genes_cosmic_e6_633.blast
    awk '{split($1, gene, "-mRNA-" ); print gene[1]}' mya_genes_cosmic_e10_4164.blast > mya_genes_cosmic_e10_4164.list
    awk '{split($1, gene, "-mRNA-" ); print gene[1]}' mya_genes_cosmic_e6_5450.blast > mya_genes_cosmic_e6_5450.list
    awk '{split($2, gene, "-mRNA-" ); print gene[1]}' top_mya_genes_cosmic_e6_633.blast > top_mya_genes_cosmic_e6_633.list
    awk '$0 ~ /^>.*/ {split($1, gene, "-mRNA-" ); split(gene[1], gene2, ">" ); print gene2[2]}' $MYAGENES/2020-09-11_Mar_genome_snap02.all.maker.proteins.fasta > all_mya_genes.list
    wc -l all_mya_genes.list # 38609 
# To compare to: steamer insertions and CN gains
    CN_GAIN=/ssd3/Mar_genome_analysis/CNV_calling/GENES
    STEAMER=/ssd3/Mar_genome_analysis/steamer/final_pipeline/genes
    # comm to compare two gene lists for overlap: 
    comm -12 <(sort $STEAMER/allCnoH_allPEI_allUSA.2000bp_upstream | uniq) <(sort mya_genes_cosmic_e6_5450.list | uniq) | wc -l # 19
    comm -12 <(sort $STEAMER/allCnoH_allPEI_allUSA.2000bp_upstream | uniq) <(sort top_mya_genes_cosmic_e6_633.list | uniq) | wc -l # 1
    comm -12 <(sort $STEAMER/allCnoH_allPEI_allUSA.2000bp_upstream | uniq) <(sort all_mya_genes.list | uniq) | wc -l # 90
    # 19/90 = 21% (vs 14%)
    # 1/90 = 1.1% (vs 1.6%)

# count number of comsmic genes in other files:
    head -100 $STEAMER/allCnoH_allPEI_allUSA.cds # NOTE THE SAMES ARE IN A DIFFERENT FORMAT SO ALL SHOW UP AS ZERO
    head -100 $STEAMER/allCnoH_allPEI_allUSA.genes
    ls $STEAMER
    LIST="allCnoH_allPEI_allUSA.samestrand.2000bp_upstream allCnoH_allPEI_allUSA.2000bp_upstream allCnoH_allPEI_allUSA.cds allCnoH_allPEI_allUSA.genes allCnoH.samestrand.2000bp_upstream allCnoH.2000bp_upstream allCnoH.cds allCnoH.genes allPEI.samestrand.2000bp_upstream allPEI.2000bp_upstream allPEI.cds allPEI.genes allUSA.samestrand.2000bp_upstream allUSA.2000bp_upstream allUSA.cds allUSA.genes"
    for i in $LIST
    do
        echo $i    
        comm -12 <(sort $STEAMER/$i | uniq) <(sort mya_genes_cosmic_e6_5450.list | uniq) | wc -l
        comm -12 <(sort $STEAMER/$i | uniq) <(sort top_mya_genes_cosmic_e6_633.list | uniq) | wc -l
        comm -12 <(sort $STEAMER/$i | uniq) <(sort all_mya_genes.list | uniq) | wc -l
    done
        # allCnoH_allPEI_allUSA.samestrand.2000bp_upstream
        # 7 = 15%
        # 0
        # 46
        # allCnoH_allPEI_allUSA.2000bp_upstream
        # 19 = 21%
        # 1
        # 90
        # SUBSTRACT TO GET DIFSTRAND: 12/44 = 27%***
        # allCnoH_allPEI_allUSA.cds
        # 0
        # 0
        # 0
        # allCnoH_allPEI_allUSA.genes
        # 29 = 15%
        # 6 = 3.2% ***
        # 190
        # allCnoH.samestrand.2000bp_upstream
        # 2 = 13%
        # 0
        # 15
        # allCnoH.2000bp_upstream
        # 11 = 29% ***
        # 1 = 2.6%
        # 38
        # SUBSTRACT TO GET DIFSTRAND: 9/23 = 39% ***
        # allCnoH.cds
        # 0
        # 0
        # 0
        # allCnoH.genes
        # 13 = 14.4%
        # 3 = 3.3%
        # 90
        # allPEI.samestrand.2000bp_upstream
        # 1 = 25%
        # 0
        # 4
        # allPEI.2000bp_upstream
        # 2 = 20%
        # 0
        # 10
        # allPEI.cds
        # 0
        # 0
        # 0
        # allPEI.genes
        # 4 = 19%
        # 1 = 4.7%
        # 21
        # allUSA.samestrand.2000bp_upstream
        # 4 = 15%
        # 0
        # 27
        # allUSA.2000bp_upstream
        # 6 = 14%
        # 0
        # 43
        # allUSA.cds
        # 0
        # 0
        # 0
        # allUSA.genes
        # 12 = 15%
        # 2
        # 80

