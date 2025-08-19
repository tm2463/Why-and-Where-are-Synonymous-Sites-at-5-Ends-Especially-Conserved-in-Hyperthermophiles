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
