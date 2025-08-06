#!/bin/bash

for file in *.txt; do
	gene=$(basename "$file" .txt)
	mafft --bl 30 "$file" > "${gene}_msa.txt"
	trimal -in "${gene}_msa.txt" -out "${gene}_trimal.txt"
done

mkdir -p cleaned_alignment
mv *_trimal.txt cleaned_alignment
