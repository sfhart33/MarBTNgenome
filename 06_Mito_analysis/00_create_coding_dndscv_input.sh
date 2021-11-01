# downloaded from published mya mt genome
head /ssd3/Mar_genome_analysis/genomes/mito/mt_genome.gff3

awk ' $3 ~ "CDS" {
    split($9, names, ";" );
    split(names[6], cds, "=" );
    split($1, genome, "." );
    length_cds = $5 - $4 + 1;
    if($7 ~ "+"){strand = 1};
    if($7 ~ "-"){strand = -1};
    print cds[2], cds[2], cds[2], genome[1], $4, $5, 1, length_cds, length_cds, strand} '  OFS='\t'  /ssd3/Mar_genome_analysis/genomes/mito/mt_genome.gff3 > /ssd3/Mar_genome_analysis/genomes/mito/mt_genome.dnds.input
