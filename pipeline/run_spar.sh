#!/bin/bash

# Selective PCR for Antibody Retrieval (SPAR) pipeline

# This is the main wrapper script for running the pipeline.
# It performs essential preprocessing steps, then designs primers.

# All steps are performed in the working directory ($WDIR).
# Working directory is typically the output directory of 10X Genomics
# cellranger, e.g. /path/to/cellranger/sample/outs.
# Intermediate files are written to the working directory.

# Final output cell_contigs_primers.csv, consisting of primers for each singlet cell, is written to the chosen output directory ($OUTDIR)

# Example of how to run this script:
# source run_spar.sh /path/to/cellranger/sample/outs my_sample_name /path/to/spar/output

# Specify inputs and outputs
WDIR=$1  # path to working directory (typically /path/to/cellranger/sample/outs/")
SAMPLE=$2  # name of sample (used as a unique identifier for each library or sample)
OUTDIR=$3  # path to output directory (where file indicating primers will be written)

# Specify path to resources
RESOURCES="../resources/"  # path to resources directory (contains immcantation singularity image)

echo "Welcome to the SPAR (Selective PCR for Antibody Retrieval) pipeline!"
echo

echo "Working directory is"
echo $WDIR
echo

date
echo

echo "Activating Anaconda environment..."
conda activate spar
echo

echo "Preparing directory..."
echo

echo "Filtering for productive contigs from singlet cells..."
python filter_reshape_contigs.py $WDIR/all_contig_annotations.csv $SAMPLE $WDIR/cell_contigs.csv
echo

echo "Converting annotations from JSON to CSV..."
python convert_annotations_json_to_csv.py $WDIR/all_contig_annotations.json $WDIR/all_contig_annotations.C_REGION.csv
echo

echo "Getting list of UMIs for each contig..."
# python make_contig_umis_json.py $WDIR/all_contig.bam $WDIR/contig_umis.json
echo

echo "Running IgBLAST to get annotations of V and J positions..."
# source igblast_immcantation.sh $WDIR $RESOURCES
echo

echo "Designing SPAR primers..."
python design_primers.py $WDIR $WDIR/cell_contigs.csv $OUTDIR
echo

echo "Done!!"
date
echo

echo "Outputs are in"
echo $OUTDIR
