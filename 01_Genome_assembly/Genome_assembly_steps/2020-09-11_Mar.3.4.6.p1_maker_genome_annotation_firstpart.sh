#!/bin/bash
# Start log

exec 1>2020-09-11_Mar.3.4.6.p1_maker_genome_annotation.log 2>&1

# Load MAKER module 

module load exonerate/2.2.0
module load repeatmasker/4.1.0
module load maker/2.31.10

# Document programs in PATH (primarily for program version ID)

date >> system_path.log
echo "" >> system_path.log
printf "%0.s-" {1..10} >> system_path.log
echo ${PATH} | tr : \\n >> system_path.log

## Establish variables for more readable code

### Paths to Maker binaries
maker=/opt/pnri/modules/sw/maker/2.31.10/x86_64-Linux-ubuntu-16.04/bin/maker
gff3_merge=/opt/pnri/modules/sw/maker/2.31.10/x86_64-Linux-ubuntu-16.04/bin/gff3_merge
fasta_merge=/opt/pnri/modules/sw/maker/2.31.10/x86_64-Linux-ubuntu-16.04/bin/fasta_merge
maker2zff=/opt/pnri/modules/sw/maker/2.31.10/x86_64-Linux-ubuntu-16.04/bin/maker2zff
fathom=~/programs/snap/fathom
forge=~/programs/snap/forge
hmmassembler=~/programs/snap/hmm-assembler.pl

### Path to Mya arenaria genome fasta file
Mar_genome=~/MAKER_Mya/2020-09-11-Mar.3.4.6.p1-MAKER/Mar.3.4.6.p1_Q30Q30A.fasta

### Path to Mya arenaria transcriptome fasta file
Mar_transcriptome=~/MAKER_Mya/2020-09-11-Mar.3.4.6.p1-MAKER/MELC-2E11_allfiles.Trinity.fasta

### Path to Crassotrea gigas NCBI protein FastA
gigas_proteome=~/MAKER_Mya/BivalveProteomes/GCF_902806645.1_cgigas_uk_roslin_v1_protein.faa

### Path to Crassostrea virginica NCBI protein fasta
virginica_proteome=~/MAKER_Mya/BivalveProteomes/GCF_002022765.2_C_virginica-3.0_protein.faa

### Path to Mytilus coruscus NCBI protein fasta
coruscus_proteome=~/MAKER_Mya/BivalveProteomes/GCA_011752425.2_MCOR1.1_protein.faa

### Path to Pecten maximus NCBI protein fasta
maximus_proteome=~/MAKER_Mya/BivalveProteomes/GCF_902652985.1_xPecMax1.1_protein.faa

### Path to Mizuhopecten yessoensis NCBI protein fasta
yessoensis_proteome=~/MAKER_Mya/BivalveProteomes/GCF_002113885.1_ASM211388v2_protein.faa

### Path to concatenated proteins fasta
bivalve_proteomes=~/MAKER_Mya/BivalveProteomes/CgiCviMcoPmaMye_protein.fasta

### Path to Mya arenaria-specific repeat library (from RepeatModeler2)
Mar_repeat_library=~/MAKER_Mya/2020-09-11-Mar.3.4.6.p1-MAKER/Mar.3.4.6.p1_Q30Q30A-families.fa

## Create Maker control files needed for running Maker
$maker -CTL

## Store path to options and exe control file
maker_opts_file=./maker_opts.ctl

## Create combined proteome FastA file, only if it doesn't already exist.
if [ ! -e "$bivalve_proteomes" ]; then
    touch "$bivalve_proteomes"
    cat "$gigas_proteome" >> "$bivalve_proteomes"
    cat "$virginica_proteome" >> "$bivalve_proteomes"
	cat "$coruscus_proteome" >> "$bivalve_proteomes"
	cat "$maximus_proteome" >> "$bivalve_proteomes"
	cat "$yessoensis_proteome" >> "$bivalve_proteomes"
fi

## Edit options file

### Set paths to Mya arenaria genome and transcriptome.
### Set path to combined proteomes.
## The use of the % symbol sets the delimiter sed uses for arguments.
## Normally, the delimiter that most examples use is a slash "/".
## But, we need to expand the variables into a full path with slashes, which screws up sed.
## Thus, the use of % symbol instead (it could be any character that is NOT present in the expanded variable; doesn't have to be "%").
sed -i "/^genome=/ s% %$Mar_genome %" "$maker_opts_file"
sed -i "/^est=/ s% %$Mar_transcriptome %" "$maker_opts_file"
sed -i "/^protein=/ s% %$bivalve_proteomes %" "$maker_opts_file"
sed -i "/^rmlib=/ s% %$Mar_repeat_library %" "$maker_opts_file"
sed -i "/^est2genome=0/ s/est2genome=0/est2genome=1/" "$maker_opts_file"
sed -i "/^protein2genome=0/ s/protein2genome=0/protein2genome=1/" "$maker_opts_file"
sed -i "/^model_org=all/ s/model_org=all/model_org=simple/" "$maker_opts_file"

## Run Maker
### Specify number of nodes to use.
mpiexec -n 40 $maker

## Merge gffs
${gff3_merge} -d Mar.3.4.6.p1_Q30Q30A.maker.output/Mar.3.4.6.p1_Q30Q30A_master_datastore_index.log

## GFF with no fasta in footer
${gff3_merge} -n -s -d Mar.3.4.6.p1_Q30Q30A.maker.output/Mar.3.4.6.p1_Q30Q30A_master_datastore_index.log > Mar.3.4.6.p1_Q30Q30A.maker.all.noseq.gff

## Merge all fastas
${fasta_merge} -d Mar.3.4.6.p1_Q30Q30A.maker.output/Mar.3.4.6.p1_Q30Q30A_master_datastore_index.log

# transcript alignments
awk '{ if ($2 == "est2genome") print $0 }' Mar.3.4.6.p1_Q30Q30A.maker.all.noseq.gff > Mar.3.4.6.p1_Q30Q30A.maker.all.est2genome.gff
est_gff1=Mar.3.4.6.p1_Q30Q30A.maker.all.est2genome.gff
# protein alignments
awk '{ if ($2 == "protein2genome") print $0 }' Mar.3.4.6.p1_Q30Q30A.maker.all.noseq.gff > Mar.3.4.6.p1_Q30Q30A.maker.all.protein2genome.gff
protein_gff1=Mar.3.4.6.p1_Q30Q30A.maker.all.protein2genome.gff
# repeat alignments
awk '{ if ($2 ~ "repeat") print $0 }' Mar.3.4.6.p1_Q30Q30A.maker.all.noseq.gff > Mar.3.4.6.p1_Q30Q30A.maker.all.repeats.gff
repeat_gff1=Mar.3.4.6.p1_Q30Q30A.maker.all.repeats.gff


