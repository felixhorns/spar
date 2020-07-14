IMG="resources/immcantation-2.6.0.simg" # path to immcantation singularity image

WDIR=$1 # working directory containing all_contig.fasta

# results will be output into the working directory

singularity exec -B $WDIR:/data $IMG changeo-igblast -s /data/all_contig.fasta -n all_contig.igblast -o /data/ -p 4 -i
