#!/bin/bash
#run salmon to quanity expression of transcripts

read -rp "enter directory path containing rnaseq data and references: " home_dir
cd "$home_dir"

for dir in *; do
	cd "$home_dir/$dir"
	
	mkdir -p salmon_index
	mkdir -p salmon_output
	
	salmon index -t GCF* -i salmon_index -k 31
	salmon quant -i salmon_index -l A -r *.fastq -o salmon_output --validateMappings
	
	cd "$home_dir"
done
