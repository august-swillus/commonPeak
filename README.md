**Overview**

This package identifies ChIP-seq peak regions that show significantly similar signal intensity across replicates from different biological conditions. It can be used to benchmark novel ChIP-seq appraoches against standard ChIP, to identify whether the novel approach identifies peaks with a comparable sensitivity and specificity to the established protocol. Additionally, it can identify stable refernce peak regions for parralel-factor ChIP normalization. 

The main function runs a 4-step pipeline:
1. Find intersecting peaks across all replicates 
2. Count reads in those peaks 
3. Subtract any input control reads from those counts 
4. Identify commonly strong peaks between conditions (DESeq2, Wald-test for equivalence of the read counts)

---

**Users can set the following arguments:**

narrow: TRUE for narrow peaks (default), FALSE for broad. 

qvalue_cutoff: q-value cutoff for filtering peaks (default 0.05).

lfc_threshold: log2 fold change threshold for DESeq2 (default 1).

output_dir: directory for results (required).

paired: TRUE if BAM files are paired-end; FALSE for single-end (default).

ma_plot: create an MA plot (default)

results_csv: create a results file listing the most common peaks with fold changes and p-value for significant similarity (default)

---

**Installation and usage example**

```r
# install.packages("devtools")   # if not installed
devtools::install_github("august-swillus/commonPeak")

library(commonPeak)

commonPeak(
  group1_peaks1 = "path/to/group1_peakFile1.bed",
  group1_peaks2 = "path/to/group1_peakFile2.bed",
  group2_peaks1 = "path/to/group2_peakFile1.bed",
  group2_peaks2 = "path/to/group2_peakFile2.bed",
  group1_bam1  = "path/to/group1_sample1.bam",
  group1_bam2 = "path/to/group1_sample2.bam"
  group2_bam1  = "path/to/group2_sample1.bam",
  group2_bam2 = "path/to/group2_sample2.bam"
  output_dir   = "results/dir"
)




