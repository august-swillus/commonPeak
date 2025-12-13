**Overview**

This package identifies ChIP-seq peak regions that show significantly similar signal intensity across replicates from different biological conditions. It can be used to benchmark novel ChIP-seq approaches against standard ChIP, to identify whether the novel approach identifies peaks with a comparable sensitivity and specificity to the established protocol. Additionally, it can identify stable reference peak regions for parallel-factor ChIP normalization. 

The pipeline consists of 4 steps:
1. Identify peak intervals present in all samples
2. Count reads in those intervals
3. Subtract input control reads from those counts 
4. Identify peaks with similar intensities in all samples across conditions (DESeq2, Wald-test for equivalence of the read counts)

---

**Users can set the following arguments:**

narrow: TRUE for narrow peaks (default), FALSE for broad. commonPeak recenters each peak interval around the point of highest read pile-up and counts reads in a fixed interval around that point. When narrow = TRUE, reads are counted when they overlap the region from −200 bp to +200 bp around the peak center. When narrow = FALSE, reads are counted in a broader interval from −500 bp to +500 bp around the peak center.

qvalue_cutoff: q-value cutoff for filtering peaks (default 0.05).

lfc_threshold: log2 fold change threshold for DESeq2 (default 1).

output_dir: directory for results (required).

paired: TRUE if BAM files are paired-end; FALSE for single-end (default).

ma_plot: create an MA plot (default)

results_csv: create a results file listing the most common peaks with fold changes and p-value for significant similarity (default)

---

**Input expectations**

The commonPeak() wrapper function expects: 

- ChIP-seq peak files (BED format)
- Corresponding BAM files for each replicate
- Optional input control BAM files 

**Importantly, the names assigned to the input files need to follow the naming convention:** 
- groupX_bedY for BED file replicate Y from group X
- groupX_bamY for BAM file replicate Y from group X
- groupX_input_bam for the input file of group X 

Samples are automatically grouped by prefix: e.g. group1_bam1 and group1_bam2 will be grouped together, and group2_bam1 and group2_bam2 as well. 

---

**Installation and usage example**

git clone https://github.com/august-swillus/commonPeak.git 
cd commonPeak 
conda env create -f environment.yml 
conda activate commonPeak

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("august-swillus/commonPeak")

library(commonPeak)

# In this example, peaks are identified that have significantly similar strenghts in group 1 and group 2, which represent different biological conditions. Each group has two replicates and an input file. 

commonPeak(
  group1_bed1        = "path/to/group1_peakFile1.bed",
  group1_bed2        = "path/to/group1_peakFile2.bed",
  group2_bed1        = "path/to/group2_peakFile1.bed",
  group2_bed2        = "path/to/group2_peakFile2.bed",
  group1_bam1        = "path/to/group1_sample1.bam",
  group1_bam2        = "path/to/group1_sample2.bam",
  group2_bam1        = "path/to/group2_sample1.bam",
  group2_bam2        = "path/to/group2_sample2.bam",
  group1_input_bam   = "path/to/group1_input.bam",
  group2_input_bam   = "path/to/group2_input.bam",
  output_dir         = "results/dir"
)
```

---

**Output file descriptions**

beds_manifest.tsv: manifest of all peak (BED) files

bams_manifest.tsv: manifest of all BAM files and input files

Intermediate files created during the multi-intersection and counting steps

normalized_counts.csv: normalized read counts per peak

deseq2_results.csv (only created if results_csv = TRUE): DESeq2 results for commonly strong peaks

MA_plot.pdf (only created if ma_plot = TRUE): MA plot of DESeq2 results

---

**Dependencies and system requirements**

Compatible with Linux and MacOS, but not Windows systems. 

All required dependencies are installed automatically. However, depending on your system, you might need to install any of the following dependencies:

R version: R ≥ 4.3

R packages (Imports): GenomicRanges, GenomicAlignments, Rsamtools, Rsubread, DESeq2, BiocParallel

System requirements: bedtools ≥ 2.30.0 (must be available in your $PATH), python3 (with pandas and numpy installed)







