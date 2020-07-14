#!/bin/bash

WDIR=$1 # path to working directory (which contains input files)
CELL_CONTIGS="data/cell_contigs.csv" # path to table containing cells and contigs
OUTDIR="outs" # path to output directory

echo "Welcome to the SPAR (Selective PCR for Antibody Retrieval) pipeline!"
date

echo "WDIR is"
echo $WDIR

echo "Activating Anaconda environment..."
source conda activate spar

echo "Preparing directory..."

echo "Converting annotations from JSON to CSV..."
date
python convert_annotations_json_to_csv.py $WDIR/all_contig_annotations.json $WDIR/all_contig_annotations.C_REGION.csv

echo "Getting list of UMIs for each contig..."
date
python make_contig_umis_json.py $WDIR/all_contig.bam $WDIR/contig_umis.json

echo "Running IgBLAST to get annotations of V and J positions..."
date
source igblast_immcantation.sh $WDIR

echo "Designing SPAR primers..."
date
python design_primers.py $WDIR $CELL_CONTIGS $OUTDIR

echo "Done!!"
date

echo "Outputs are in"
echo $OUTDIR
