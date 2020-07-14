# SPAR (Selective PCR for Antibody Retrieval)

SPAR (Selective PCR for Antibody Retreival) is a strategy for cloning antibodies in single cells within pooled sequence libraries. The method is described in a paper published in PLOS ONE in 2020 (URL pending).

This repository contains software for designing primers for SPAR and Jupyter Notebooks that reproduce the figures shown in the paper.

This software is designed for 10X Genomics Single Cell 5' V(D)J libraries. 

## Citation

Publication pending.

## Environment

A complete Python 2.7 environment is specified by the YAML file in the repository, which can be used to create an environment using Anaconda.

Immcantation 2.6.0 is required to run igblast and can be installed using Docker or Singularity using the instructions [here](https://immcantation.readthedocs.io/en/version-2.4.0/docker/intro.html).

## How to run the primer design pipeline

The pipeline can be run using `run_spar.sh`. The inputs are:

- Path to 10X Genomics cellranger output directory, which contains `all_contig_annotations.json`, `all_contig.fasta`, and `all_contig.bam`.
- `cell_contigs.csv`, which is a table containing all cells and contigs generated by `make_cell_contigs.csv`.

## Contents

### `pipeline`

`run_spar.sh`: Runs the primer design pipeline.

`convert_annotations_json_to_csv.py`: Converts annotations of contigs from JSON to CSV.

`make_contig_umis_json.py`: Gets all UMI sequences associated with each contig.

`igblast_immcantation.sh`: Runs igblast to annotate contigs.

`design_primers.py`: Designs primers. Wrapper for core algorithm for primer design.

`spar.py`: Python module containing core algorithm for primer design.

### `analysis`

`figures.ipynb`: Analysis of primer features and retrievability of the repertoire. Generates figures shown in paper. 

## Disclaimer
This project is not maintained. Software is provided as is and requests for support may not be addressed.

## Contact
If you have questions or comments, please contact Felix Horns at <rfhorns@gmail.com>.
