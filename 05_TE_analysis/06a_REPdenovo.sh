# Run REPdenovo
	module load repdenovo
	REPDENOVO=/opt/pnri/modules/sw/repdenovo/2019.07.20/x86_64-Linux-ubuntu-16.04
	/ssd2/REPdenovo
	cd /ssd2/REPdenovo
	export TMPDIR=/var/tmp/
# MELC-2E11 (Healthy Reference)
	python $REPDENOVO/main.py -c Assembly -g /ssd2/testfiles/REPdenovo/MELC-2E11_config.txt -r /ssd2/testfiles/REPdenovo/MELC-2E11_files.txt
	python $REPDENOVO/main.py -c Scaffolding -g /ssd2/testfiles/REPdenovo/MELC-2E11_config.txt -r /ssd2/testfiles/REPdenovo/MELC-2E11_files.txt
# MELC-A11 (USA BTN)
	python $REPDENOVO/main.py -c Assembly -g /ssd2/REPdenovo/MELC-A11_config.txt -r /ssd2/REPdenovo/MELC-A11_files.txt 
	python $REPDENOVO/main.py -c Scaffolding -g /ssd2/REPdenovo/MELC-A11_config.txt -r /ssd2/REPdenovo/MELC-A11_files.txt
# PEI-DN08 (PEI BTN)
	python $REPDENOVO/main.py -c All -g /ssd2/REPdenovo/PEI-DN08_config.txt -r /ssd2/REPdenovo/PEI-DN08_files.txt 

# run repeatclassifier
	module load repeatmodeler
	cd /ssd2/testfiles/REPdenovo/output/
	RepeatClassifier -consensi contigs.fa
	cd /ssd2/REPdenovo/MELC_A11_output/
	RepeatClassifier -consensi contigs.fa
	cd /ssd2/REPdenovo/PEI-DN08_output/
	RepeatClassifier -consensi contigs.fa

