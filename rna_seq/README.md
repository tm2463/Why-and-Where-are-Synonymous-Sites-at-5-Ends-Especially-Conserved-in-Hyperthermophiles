RNA-Seq

1. get_links.py
- Get links to search sra for rna-seq data for each archaeal species in assembly summary file
2. sra_setup.py
- Create master file sra-tools can parse to download FASTQ files
3. get_fastq.sh
- Download all FASTQ files off the SRA
4. get_refs.py
- Create file with links to all reference genomes for relevant species
5. get_salmon_refs.py
- Download references for Salmon
6. salmon_commands.sh
- Run salmon to quanitfy gene expression
7. expression_analysis.py
- Perform analysis on quant.sf output
