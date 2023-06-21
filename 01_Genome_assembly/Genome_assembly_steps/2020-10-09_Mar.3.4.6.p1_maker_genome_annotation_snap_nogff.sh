#!/bin/bash
# Start log

#script for snap 1 and 2 after typo in snap line in earlier script

exec 1>2020-10-09_Mar.3.4.6.p1_maker_genome_annotation_snap.log 2>&1

# Load MAKER module 

module load exonerate/2.2.0
module load repeatmasker/4.1.0
module load maker/2.31.10

# Document programs in PATH (primarily for program version ID)

date >> system_path_snap.log
echo "" >> system_path_snap.log
printf "%0.s-" {1..10} >> system_path_snap.log
echo ${PATH} | tr : \\n >> system_path_snap.log

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

## Run SNAP training, round 1
mkdir snap01 && cd snap01
${maker2zff} ../Mar.3.4.6.p1_Q30Q30A.all.gff
${fathom} -categorize 1000 genome.ann genome.dna
${fathom} -export 1000 -plus uni.ann uni.dna
${forge} export.ann export.dna
${hmmassembler} test_snap1 . > 2020-09-11_Mar_snap01.hmm

## Initiate second Maker run.
### Copy initial maker control files and
### - change gene prediction settings to 0 (i.e. don't generate Maker gene predictions)
### - set location of snaphmm file to use for gene prediction
cp ../maker_* .
sed -i "/^est2genome=1/ s/est2genome=1/est2genome=0/" maker_opts.ctl
sed -i "/^protein2genome=1/ s/protein2genome=1/protein2genome=0/" maker_opts.ctl
sed -i "/^snaphmm=/ s% %2020-09-11_Mar_snap01.hmm %" maker_opts.ctl

## Run Maker
### Set basename of files and specify number of CPUs to use
mpiexec -n 40 $maker \
-base 2020-09-11_Mar_genome_snap01

## Merge gffs
${gff3_merge} -d 2020-09-11_Mar_genome_snap01.maker.output/2020-09-11_Mar_genome_snap01_master_datastore_index.log

### GFF with no fasta in footer
${gff3_merge} -n -s -d 2020-09-11_Mar_genome_snap01.maker.output/2020-09-11_Mar_genome_snap01_master_datastore_index.log > 2020-09-11_Mar_genome_snap01.all.noseq.gff

### Merge all fastas
${fasta_merge} -d 2020-09-11_Mar_genome_snap01.maker.output/2020-09-11_Mar_genome_snap01_master_datastore_index.log

# transcript alignments
awk '{ if ($2 == "est2genome") print $0 }' 2020-09-11_Mar_genome_snap01.noseq.gff > 2020-09-11_Mar_genome_snap01.est2genome.gff
est_gff2=~/MAKER_Mya/2020-09-11-Mar.3.4.6.p1-MAKER/snap01/2020-09-11_Mar_genome_snap01.est2genome.gff
# protein alignments
awk '{ if ($2 == "protein2genome") print $0 }' 2020-09-11_Mar_genome_snap01.noseq.gff > 2020-09-11_Mar_genome_snap01.protein2genome.gff
protein_gff2=~/MAKER_Mya/2020-09-11-Mar.3.4.6.p1-MAKER/snap01/2020-09-11_Mar_genome_snap01.protein2genome.gff
# repeat alignments
awk '{ if ($2 ~ "repeat") print $0 }' 2020-09-11_Mar_genome_snap01.noseq.gff > 2020-09-11_Mar_genome_snap01.repeats.gff
repeats_gff2=~/MAKER_Mya/2020-09-11-Mar.3.4.6.p1-MAKER/snap01/2020-09-11_Mar_genome_snap01.repeats.gff

## Run SNAP training, round 2
cd ..
mkdir snap02 && cd snap02
${maker2zff} ../snap01/2020-09-11_Mar_genome_snap01.all.gff
${fathom} -categorize 1000 genome.ann genome.dna
${fathom} -export 1000 -plus uni.ann uni.dna
${forge} export.ann export.dna
${hmmassembler} test_snap1 . > 2020-09-11_Mar_snap02.hmm

## Initiate third and final Maker run.
### Copy initial maker control files and:
### - change gene prediction settings to 0 (i.e. don't generate Maker gene predictions)
### - set location of snaphmm file to use for gene prediction
cp ../maker_* .
sed -i "/^est2genome=1/ s/est2genome=1/est2genome=0/" maker_opts.ctl
sed -i "/^protein2genome=1/ s/protein2genome=1/protein2genome=0/" maker_opts.ctl
sed -i "/^snaphmm=/ s% %2020-09-11_Mar_snap02.hmm %" maker_opts.ctl

## Run Maker
### Set basename of files and specify number of CPUs to use
mpiexec -n 40 $maker \
-base 2020-09-11_Mar_genome_snap02

## Merge gffs
${gff3_merge} \
-d 2020-09-11_Mar_genome_snap02.maker.output/2020-09-11_Mar_genome_snap02_master_datastore_index.log

### GFF with no fasta in footer
${gff3_merge} -n -s -d 2020-09-11_Mar_genome_snap02.maker.output/2020-09-11_Mar_genome_snap02_master_datastore_index.log > 2020-09-11_Mar_genome_snap02.all.noseq.gff

### Merge all fastas
${fasta_merge} -d 2020-09-11_Mar_genome_snap02.maker.output/2020-09-11_Mar_genome_snap02_master_datastore_index.log
