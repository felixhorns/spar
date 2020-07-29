# Preprocess single-cell data for SPAR

# This script performs the following preprocessing steps:
# 1) Filter for high-quality contigs using the productive flag.
# 2) Filter for singlet cells by choosing those having exactly 1 productive IGH and 1 productive IGK/L.
# 3) Select contigs that belong to those singlet cells and are productive.
# 4) Reshape contigs into a dataframe in which each row consists of a singlet cell,
# and that row contains the information about the IGH contig and the IGK/L contig associated with that cell.
# 5) Write the reshaped dataframe to output file.

# Input is the `all_contigs_annotations.csv` from 10X Genomics cellranger.
# Output is a filtered and reshaped list of contigs belonging to singlet cells.

import sys
import pandas as pd

# Specify inputs and outputs
infile = sys.argv[1]  # contig annotation file from 10X Genomics cellranger output.
sample_name = sys.argv[2]  # name of sample
outfile = sys.argv[3]  # filtered, reshaped list of productive contigs belonging to singlet cells, typically cell_contigs.csv.

# Load data
contigs_all = pd.read_csv(infile, header=0)
print "Total contigs:", contigs_all.shape[0]

# Set sample name
contigs_all["sample"] = sample_name

# Filter for high-quality contigs using productive flag
# Productive implies full-length, removes all TCR and None chain assignments, removes most "Multi" chain assignments.
# Productive flag is sufficient to filter for high quality contigs.
contigs_aggr_filtered = contigs_all.groupby(["productive", "sample", "barcode", "chain"]).size().unstack(fill_value=0).xs("True", level="productive")
print "Productive contigs:", contigs_aggr_filtered.shape[0]

# Filter for singlet cells (having exactly 1 IGH and 1 IGK/L)
singlets = contigs_aggr_filtered.loc[(contigs_aggr_filtered["IGH"] == 1) &
                                      (((contigs_aggr_filtered["IGL"] == 1) &
                                        (contigs_aggr_filtered["IGK"] == 0)) |
                                       ((contigs_aggr_filtered["IGL"] == 0) &
                                        (contigs_aggr_filtered["IGK"] == 1)))]

# Select contigs belonging to these singlets that are productive and IGH/IGK/IGL.

# Filter contigs for these cells
contigs_valid = contigs_all.set_index(["sample", "barcode"]).loc[singlets.index]

# Filter contigs for only IGH, IGL, or IGK
contigs_valid = contigs_valid.loc[contigs_valid["chain"].isin(["IGH", "IGL", "IGK"])]

# Filter contigs for only productive
contigs_valid = contigs_valid.loc[contigs_valid["productive"] == "True"]

# Make dataframe of cells and associated heavy and light chain contigs
# Each row is a singlet cell and contains the heavy and light chain contig of that cell.

# Select IGH contigs
selector = contigs_valid["chain"] == "IGH"
df_tmp_IGH = contigs_valid[selector].copy(deep=True)

# Select IGK/L contigs
selector = contigs_valid["chain"].isin(["IGK", "IGL"])
df_tmp_IGKL = contigs_valid[selector].copy(deep=True)
df_tmp_IGKL.shape

# Choose columns
columns_left = ["is_cell", "contig_id", "high_confidence", "length", "chain", "full_length", "productive", "reads", "umis"]
columns_right = ["contig_id", "high_confidence", "length", "chain", "full_length", "productive", "reads", "umis"]
suffixes = ["_IGH", "_IGKL"]

# Merge IGH and IGK/L contigs, resulting in one row per cell
df_cell_contigs = df_tmp_IGH[columns_left].merge(df_tmp_IGKL[columns_right], left_index=True, right_index=True, suffixes=suffixes, validate="one_to_one")
print "Cells after filtering (singlet cells with productive contigs):", df_cell_contigs.shape[0]

# Write to file
df_cell_contigs.to_csv(outfile)
