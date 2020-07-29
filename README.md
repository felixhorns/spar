# SPAR (Selective PCR for Antibody Retrieval)

SPAR (Selective PCR for Antibody Retreival) is a strategy for cloning antibodies in single cells within pooled sequence libraries. The method is described in a paper published in PLOS ONE in 2020 (URL pending).

This repository contains software for designing primers for SPAR and Jupyter Notebooks that reproduce the figures shown in the paper.

This software is designed for 10X Genomics Single Cell 5' V(D)J libraries. 

## Citation

Publication pending.

## Environment

A complete Python 2.7 environment is specified by the YAML file in the repository, which can be used to create an environment using Anaconda. Key dependencies are `pandas`, `primer3`, and `pysam`.

Immcantation 2.6.0 is required to run igblast and can be installed using Docker or Singularity using the instructions [here](https://immcantation.readthedocs.io/en/version-2.4.0/docker/intro.html). Containerization such as Docker or Singularity greatly simplifies installation of igblast and is highly recommended. The pipeline is currently configured to use Singularity. Singularity image for immcantation must be located at `resources/immcantation-2.6.0.simg`.

## Running the SPAR primer design pipeline

The wrapper script that runs the pipeline is `run_spar.sh`.

Inputs:
- Path to the working directory. This must contain `all_contig_annotations.csv`, `all_contig_annotations.json`, `all_contig.fasta`, and `all_contig.bam`. This is typically the output directory of 10X Genomics cellranger (e.g. `/path/to/cellranger/sample/outs`).
- Sample name. This is used as a unique identifier for each library or sample. Typically, this is used to distinguish different 10X libraries.
- Path to the output directory. This is where the file indicating SPAR primers will be written.

Output is a table, called `cell_contigs_primers.csv`, listing the optimal primers to retrieve the variable region of each contig.

Example of running SPAR pipeline:

`./run_spar.sh /path/to/cellranger/sample/outs my_sample_name /path/to/spar/output`

## Contents

### `pipeline`

`run_spar.sh`: Runs the primer design pipeline.

`filter_reshape_contigs.py`: Pre-processes single-cell data for SPAR. Filters for productive contigs and singlet cells. Reshapes data into a table where each row is a cell, and the row contains the heavy and light chain contigs for that cell.

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
