# Convert annotations in JSON to CSV

import sys
import json

def extract_C_REGION_contig_match_start(annotation):
    """ Find position where C region starts """
    
    # find annotation of C region (if it exists)
    annotation_C_REGION = None
    for x in annotation["annotations"]:
        if x["feature"]["region_type"] == "C-REGION":
            annotation_C_REGION = x
            
    # extract C region boundary
    if annotation_C_REGION is not None:
        annotation_C_REGION_contig_match_start = annotation_C_REGION["contig_match_start"]
    else:
        annotation_C_REGION_contig_match_start = "None"

    return annotation_C_REGION_contig_match_start

infile = sys.argv[1]
outfile = sys.argv[2]

annotations = json.load(open(infile))

out = open(outfile, 'w')

header = ["barcode", "contig_name", "is_cell", "high_confidence", "productive", "C_REGION_contig_match_start", "sequence", "aa_sequence"]
out.write(",".join(header) + "\n")

for annotation in annotations:

    annotation_C_REGION_contig_match_start = extract_C_REGION_contig_match_start(annotation)

    out_line = ",".join(map(str, [annotation["barcode"], annotation["contig_name"],
                                  annotation["is_cell"], annotation["high_confidence"],
                                  annotation["productive"], annotation_C_REGION_contig_match_start,
                                  annotation["sequence"], annotation["aa_sequence"]]))
    out.write(out_line + "\n")

out.close()
