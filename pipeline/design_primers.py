# Design primers for SPAR

import sys
import time
import json
import pandas as pd

import spar

if __name__ == "__main__":

    ##### Read arguments 

    wdir = sys.argv[1] # directory that contains files below
    infile_cell_contigs = sys.argv[2] # cell_contigs.csv which contains name of IGH and IGKL contigs for each cell
    outdir = sys.argv[3] # output directory where cell_contigs_primers.csv will be written

    infile_contigs_seqs_all = wdir + "/all_contig_annotations.C_REGION.csv"
    infile_contig_umis = wdir + "/contig_umis.json"
    infile_igblast = wdir + "/all_contig.igblast_db-pass.tab"

    outfile = outdir + "/cell_contigs_primers.csv"

    ##### Load inputs

    cell_contigs = pd.read_csv(infile_cell_contigs, header=0) # table indicating IGH, IGKL contigs for all cells
    
    contigs_seqs_all = pd.read_csv(infile_contigs_seqs_all, header=0) # sequences with annotation of C region

    contig_umis = json.load(open(infile_contig_umis)) # map from contig_id to UMI sequences

    usecols = ["SEQUENCE_ID", "V_SEQ_START", "V_SEQ_LENGTH", "J_SEQ_START", "J_SEQ_LENGTH", "SEQUENCE_VDJ"]
    contigs_igblast = pd.read_csv(infile_igblast, header=0, sep="\t", index_col=0, usecols=usecols) # igblast results with annotation of V, J positions

    ##### Design primers

    # Design primers for PCR1
    cell_contigs_primers = spar.design_primers_contigs_pcr1_IGH_IGKL(cell_contigs, contigs_seqs_all, contig_umis, verbose=False)

    # Design primers for PCR2 
    cell_contigs_primers = spar.design_primers_contigs_pcr2_IGH_IGKL(cell_contigs_primers, contigs_seqs_all, contigs_igblast, verbose=False)

    ##### Write result to file
    cell_contigs_primers.to_csv(outfile)

    print "Done!!"
