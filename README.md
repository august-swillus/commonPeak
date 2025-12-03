This package identifies ChIP-seq peak regions that show significantly similar signal intensity across replicates from different biological conditions. 

This function runs a 4-step pipeline:
1. Find intersecting peaks across all replicates 
2. Count reads in those peaks 
3. Subtract any input control reads from those counts 
4. Identify commonly strong peaks between conditions (DESeq2, Wald-test for equivalence of the read counts)

Users can set the following arguments: 

narrow: TRUE for narrow peaks (default), FALSE for broad. 
qvalue_cutoff: q-value cutoff for filtering peaks (default 0.05).
lfc_threshold: log2 fold change threshold for DESeq2 (default 1).
output_dir: directory for results (required).
paired: TRUE if BAM files are paired-end; FALSE for single-end (default).
ma_plot create an MA plot (default)
results_csv create a results file listing the most common peaks (default)
