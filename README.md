# Why-and-Where-are-Synonymous-Sites-at-5-Ends-Especially-Conserved-in-Hyperthermophiles

Required Packages:
- Python:
   - Requests
   - NumPy
   - pandas
   - bs4
   - Bio
   - tqdm
   - ViennaRNA Python
   - MatPlotLib
   - Seaborn
   - SciPy
   - statsmodels
   - Selenium
- R:
   - ape
   - nlme
   - geiger
   - phangorn
   - ggplot2
   - dpylr
- Command Line:
   - sra-tools
   - Salmon
   - IQTree
   - trimAL

Running Order:

n.b. Halophiles were removed from analysis for reasons stated in report. To remove ANY genera from analysis, delete the entry from the growth temps files generated from growth_temps.py and proceed with analysis.

Database Builder:
1. fetch_data.py
   - Downloads the required cds and genome files from the NCBI FTP server
2. growth_temps.py
   - Output optimum growth temperature predictions from genome files

GC3:
1. genome_gc.R
   - Comparison of whole genome GC content for temp groups
2. codon_gc3.py
   - Tracking changes in codon GC3 content in the 5' to 3' direction
   - Visualised with gc3.R
3. codon_usage.py
   - Looking at nucleotide bias at the third codon position
   - Outputs are stiched together and fed into nn3_bias.py
   - visualised with nn3.R

Permutation Analysis:
1. elongate_cds.py
   - Extracts genome sequences -30nt to +60nt relative to start for each gene
   - Extracts intergenome for each genus
2. synonymous_shuffle.py
   - Creates synonymous permutations of extracted genome sequences biased by intergenome nucleotide frequencies
3. mean_permutations.py (this may take days to run)
   - Makes secondary structure and folding energy predictions for permuted sequences
   - Returns average permuted folding energy for each gene
4. real_fold.py
   - Makes secondary structure predictions for real sequences for comparison against permuted sequences
5. permutation_analysis.py
   - Compares synonymous permutation free energy predictions to real sequence predictions
   - Visualised with perm_graph.R
6. start_availability.py
   - Predict which bases in start codon are free/occluded
   - Data analysis with free_start.py
   - Stats with start_ttest.R
   - Visualised with box_plot.R


Phylogenetics:
1. core_genes.py
   - List all core genes for building phylogeny
2. gene_multi_fasta.py
   - Collect all core genes into multi-fasta file for alignment
3. mafft.sh
   - Independently align and trim core genes
4. concatenate.py
   - Concatenate MSA's and create partitions file
5. iqtree.sh
   - Command used to build tree in IQTree
6. PGLS
   - PGLS metadata collected in pgls_metadata.py and pgls2_metadata.py
   - PGLS performed and associated graphs created with pgls.R and pgls2.R


RNA-Seq Analysis:
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
