#!/bin/bash

for file in *.txt; do
	gene=$(basename "$file" .txt)
	mafft --bl 30 "file" > "${gene}_msa.txt"
done