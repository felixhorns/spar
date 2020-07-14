# Get all UMI sequences associated with each contig_id in a BAM file

import sys
import pysam
import json

infile = sys.argv[1]
outfile = sys.argv[2]

# Get UMIs associated with each contig_id in the file

i = 0

contig_umis = {}

f = pysam.AlignmentFile(infile, "rb")

for read in f.fetch():

    # check if read has a UMI
    if read.has_tag("UB"):
        umi = read.get_tag("UB")
    else:
        # if it does not have a UMI tag, skip it
        continue

    # add UMI to set for contig_id
    if read.reference_name in contig_umis:
        contig_umis[read.reference_name].add(umi)
    else:
        contig_umis[read.reference_name] = set([umi])

    # report how many reads have been processed
    i += 1
    if i % 1000000 == 0:
        print i, "reads processed"

# Write result to file

# convert to dictionary of lists
contig_umis = dict(zip(contig_umis.keys(), map(list, contig_umis.values())))

# make a JSON string out of dictionary
out_str = json.dumps(contig_umis)

# write to file
with open(outfile, 'w') as f:
    f.write(out_str)

print "Done!!"
