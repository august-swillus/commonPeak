#!/usr/bin/env python3
import pandas as pd
import sys, os
import pysam
import numpy as np

# --- Read inputs ---
bam_manifest = sys.argv[1]
output_dir   = sys.argv[2]

bam_df = pd.read_csv(bam_manifest, sep="\t")

# --- Read count matrix ---
count_path = os.path.join(output_dir, "read_counts_matrix.csv")
counts = pd.read_csv(count_path, index_col=0)

# --- Prepare output ---
scaled_counts = pd.DataFrame(index=counts.index)

# --- Split inputs and ChIP samples ---
input_bams = bam_df[bam_df["is_input"] == True]
chip_bams  = bam_df[bam_df["is_input"] == False]

if input_bams.empty:
    print("⚠️  No input samples detected — skipping input subtraction.")
else:
    print("🧭 Found input controls for groups:", input_bams["group"].unique().tolist())

# --- Compute library sizes from BAM files ---
lib_sizes = {}
print("\n🧮 Computing library sizes for all BAM files ...")

for _, row in bam_df.iterrows():
    bam_path = row["path"]
    sample_name = os.path.basename(bam_path).replace(".bam", "")
    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            lib_sizes[sample_name] = bam.mapped
    except Exception as e:
        print(f"⚠️  Could not read {bam_path}: {e}")
        lib_sizes[sample_name] = 0

lib_sizes_df = pd.DataFrame.from_dict(lib_sizes, orient="index", columns=["LibrarySize"])
lib_sizes_df.to_csv(os.path.join(output_dir, "library_sizes.csv"))
print(f"✅ Library sizes saved to: {output_dir}/library_sizes.csv\n")

# --- Compute mean library size and DiffBind-style scaling factors ---
mean_lib_size = np.mean([s for s in lib_sizes.values() if s > 0])
scaling_factors = {sample: mean_lib_size / size if size > 0 else 1.0
                   for sample, size in lib_sizes.items()}

scaling_df = pd.DataFrame({
    "LibrarySize": [lib_sizes[s] for s in lib_sizes],
    "ScalingFactor": [scaling_factors[s] for s in lib_sizes],
    "LibSizeAfterScaling": [mean_lib_size if lib_sizes[s] > 0 else np.nan for s in lib_sizes]
}, index=list(lib_sizes.keys()))

scaling_df.to_csv(os.path.join(output_dir, "scaling_factors.csv"))
print(f"📊 Global mean library size: {mean_lib_size:,.0f}")
print("✅ Scaling factors written to: scaling_factors.csv\n")

# --- Print scaling factors for ChIP samples only ---
print("📈 Scaling factors (ChIP samples only):")
for _, row in chip_bams.iterrows():
    chip_key = os.path.basename(row["path"]).replace(".bam", "")
    sf = scaling_factors.get(chip_key, 1.0)
    print(f"   {chip_key}: scaling_factor={sf:.6f}")

# --- Apply subtraction and scaling ---
print("\n⚙️  Subtracting input and scaling to mean library size ...")

for group in chip_bams["group"].unique():
    group_chips  = chip_bams[chip_bams["group"] == group]
    group_inputs = input_bams[input_bams["group"] == group]

    if group_inputs.empty:
        print(f"⚠️  No input found for {group}; only scaling ChIP samples.")
        for _, chip_row in group_chips.iterrows():
            chip_key = os.path.basename(chip_row["path"]).replace(".bam", "")
            if chip_key in counts.columns:
                sf = scaling_factors.get(chip_key, 1.0)
                scaled_counts[chip_key] = (counts[chip_key] * sf).round()
        continue

    input_key = os.path.basename(group_inputs.iloc[0]["path"]).replace(".bam", "")
    print(f"🔗 Group {group}: using input {input_key}")

    for _, chip_row in group_chips.iterrows():
        chip_key = os.path.basename(chip_row["path"]).replace(".bam", "")

        if chip_key not in counts.columns or input_key not in counts.columns:
            print(f"⚠️  Missing columns for {chip_key} or {input_key}, skipping.")
            continue

        sf = scaling_factors.get(chip_key, 1.0)
        adj = counts[chip_key] - counts[input_key]
        adj[adj < 0] = 0
        scaled_counts[chip_key] = (adj * sf).round()

# --- Save normalized matrix ---
out_csv = os.path.join(output_dir, "normalized_counts.csv")
scaled_counts.to_csv(out_csv)
print(f"\n✅ Normalized (DiffBind-style) count matrix written to: {out_csv}\n")
