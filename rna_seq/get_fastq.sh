#!/bin/bash
#download rna_seq data

input="master.tsv"
output_dir="./rnaseq_data"

mkdir -p "$output_dir"

while IFS=$'\t' read -r run name; do
	sub_dir="${output_dir}/${name}"
	mkdir -p "$sub_dir"
	
	echo "downloading $run -> $name"
	fasterq-dump "$run" --split-files --outdir "$sub_dir" --threads 8
 
done <"$input"
