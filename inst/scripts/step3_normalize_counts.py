#!/usr/bin/env python3

import os
import sys
import pandas as pd

if len(sys.argv) < 3:
    sys.exit("Usage: python3 step3_normalize_counts.py <bam_manifest.tsv> <output_dir>")

bam_manifest = sys.argv[1]
output_dir = sys.argv[2]

bam_df = pd.read_csv(bam_manifest, sep="\t")

count_path = os.path.join(output_dir, "read_counts_matrix.csv")
counts = pd.read_csv(count_path, index_col=0)

subtracted_counts = pd.DataFrame(index=counts.index)

input_bams = bam_df[bam_df["is_input"] == True]
chip_bams  = bam_df[bam_df["is_input"] == False]

if input_bams.empty:
    print("No input controls found — returning raw counts unchanged.")
    subtracted_counts = counts.copy()
else:
    print("Found input controls for groups:", input_bams["group"].unique().tolist())
    print("Subtracting input reads from ChIP samples ...")

    for group in chip_bams["group"].unique():
        group_chips  = chip_bams[chip_bams["group"] == group]
        group_inputs = input_bams[input_bams["group"] == group]

        if group_inputs.empty:
            print(f"No input found for {group}; keeping raw counts for its ChIP samples.")
            for _, chip_row in group_chips.iterrows():
                chip_key = os.path.basename(chip_row["path"]).replace(".bam", "")
                if chip_key in counts.columns:
                    subtracted_counts[chip_key] = counts[chip_key]
            continue

        input_key = os.path.basename(group_inputs.iloc[0]["path"]).replace(".bam", "")
        print(f"Group {group}: using input {input_key}")

        for _, chip_row in group_chips.iterrows():
            chip_key = os.path.basename(chip_row["path"]).replace(".bam", "")
            if chip_key not in counts.columns or input_key not in counts.columns:
                print(f"Missing columns for {chip_key} or {input_key}, skipping.")
                continue

            adj = counts[chip_key] - counts[input_key]
            adj[adj < 0] = 0
            subtracted_counts[chip_key] = adj

out_csv = os.path.join(output_dir, "normalized_counts.csv")
subtracted_counts.to_csv(out_csv)
print(f"Input-subtracted count matrix written to: {out_csv}")
