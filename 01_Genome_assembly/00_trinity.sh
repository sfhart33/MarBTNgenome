module load trinity 
cd /ssd2/Trinity/All_samples_together

# Merge R1 and R2 files to run Trinity on full data set
	zcat MELC-2E11-1-mantle_R1_001.fastq.gz MELC-2E11-2-foot_R1_001.fastq.gz MELC-2E11-3-siphon_R1_001.fastq.gz MELC-2E11-5-muscle_R1_001.fastq.gz MELC-2E11-6-gills_R1_001.fastq.gz MELC-2E11-7-hemocytes_R1_001.fastq.gz | gzip -c > MELC-2E11_R1_allfiles-cat.fastq.gz
	zcat MELC-2E11-1-mantle_R2_001.fastq.gz MELC-2E11-2-foot_R2_001.fastq.gz MELC-2E11-3-siphon_R2_001.fastq.gz MELC-2E11-5-muscle_R2_001.fastq.gz MELC-2E11-6-gills_R2_001.fastq.gz MELC-2E11-7-hemocytes_R2_001.fastq.gz | gzip -c > MELC-2E11_R2_allfiles-cat.fastq.gz

# Run trinity
	Trinity --seqType fq --max_memory 200G --CPU 16 --trimmomatic --full_cleanup --left MELC-2E11_R1_allfiles-cat.fastq.gz --right MELC-2E11_R2_allfiles-cat.fastq.gz

