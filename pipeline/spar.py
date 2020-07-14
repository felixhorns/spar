import time
import json
import numpy as np
import pandas as pd
import primer3

def design_primers_contigs_pcr1_IGH_IGKL(cell_contigs, contigs_seqs_all, contig_umis, verbose=False):
    """ Design primers for PCR1 for IGH and IGKL for all contigs in cell_contigs """
    
    if verbose:
        print "Starting PCR1 primer design"
        print
    
    # IGH
    
    if verbose:
        print "Starting PCR1 primer design for IGH..."
        
    start_time = time.time()
    
    contig_ids = cell_contigs["contig_id_IGH"] # Get contig_ids
    contigs_seqs_all_focal = contigs_seqs_all.set_index("contig_name").loc[contig_ids] # get sequences
    cell_contigs_primers = _design_primers_contigs_pcr1(contig_ids, contigs_seqs_all_focal, contig_umis, cell_contigs, left_on="contig_id_IGH", column_prefix="IGH_pcr1_", failure_cause_suffix="~IGH_pcr1", verbose=verbose)

    if verbose:
        print "Design finished for IGH"
        print "Elapsed time for PCR1 IGH primer design (s):", time.time() - start_time
        print
        
    # IGKL
    
    if verbose:
        print "Starting PCR1 primer design for IGKL..."
        
    start_time = time.time()
    
    contig_ids = cell_contigs["contig_id_IGKL"] # Get contig_ids
    contigs_seqs_all_focal = contigs_seqs_all.set_index("contig_name").loc[contig_ids] # get sequences
    cell_contigs_primers = _design_primers_contigs_pcr1(contig_ids, contigs_seqs_all_focal, contig_umis, cell_contigs_primers, left_on="contig_id_IGKL", column_prefix="IGKL_pcr1_", failure_cause_suffix="~IGKL_pcr1", verbose=verbose)
    
    if verbose:
        print "Design finished for IGKL"
        print "Elapsed time for PCR1 IGKL primer design (s):", time.time() - start_time
        print

    return cell_contigs_primers

def _design_primers_contigs_pcr1(contig_ids, contigs_seqs, contig_umis, cell_contigs, left_on="", column_prefix="", failure_cause_suffix="", verbose=False):
    """ Design primers for contigs specified by contig_ids, merge results into dataframe cell_contigs """
    
    ### Find and drop contigs with no C region match
    
    # Find contigs with no C region match
    contig_ids_no_C_region_match = contigs_seqs.loc[contigs_seqs["C_REGION_contig_match_start"] == "None"].index

    # Drop contigs with no C region match
    contigs_seqs = contigs_seqs.loc[contigs_seqs["C_REGION_contig_match_start"] != "None"]
    contig_ids = contigs_seqs.index # Update contig ids to remove dropped contigs

    print "Contigs with no C region match:", len(contig_ids_no_C_region_match)
    print "Remaining contigs (after dropping those with no C region match):", len(contig_ids)
    
    ### Prepare dictionaries for looking up sequences and sequence features

    seqs = dict(zip(contigs_seqs.index, contigs_seqs.sequence)) # map from contig_id to sequence
    C_region_starts = dict(zip(contigs_seqs.index, contigs_seqs.C_REGION_contig_match_start)) # map from contig_id to C region start
    cbcs_all = [x.split("-")[0] for x in contigs_seqs.index]
    cbcs = dict(zip(contigs_seqs.index, cbcs_all)) # map from contig_id to cell barcode (cbc)

    # map from contig_id to full molecule sequence for all UMIs
    partial_read1 = "CTACACGACGCTCTTCCGATCT" # partial read 1 sequence (5' of CBC)
    tso = "TTTCTTATATGGG" # TSO sequence (3' of UMI)
    seqs_full_cdna = {}

    for contig_id in contig_ids:
        seq = seqs[contig_id]
        cbc = cbcs[contig_id]
        umis = sorted(contig_umis[contig_id])
        for umi in umis:
            contig_id_umi = str(contig_id + "~" + umi)
            seqs_full_cdna[contig_id_umi] = str(partial_read1 + cbc + umi + tso + seq)

    ### Design primers

    best_primer_pairs = {}
        
    for i, contig_id in enumerate(contig_ids):
        best_primer_pairs[contig_id] = _find_primers_pcr1(contig_id, cbcs, contig_umis, C_region_starts, seqs_full_cdna, partial_read1, tso)

        if verbose:
            if i % 20 == 0:
                print i, "pairs processed"

    if verbose:
        print "Finished with primer design"

    ### Merge results into cell contig dataframe
    
    # Reform all primer results into an appropriate data frame for merging

    res = {}
    index = []
    
    contig_ids_no_primer_found = [] # keep track of contigs that had no primer pair

    for contig_id, best_primer_pair in best_primer_pairs.items():
        
        if best_primer_pair is None: 
            contig_ids_no_primer_found.append(contig_id)
            continue
            
        index.append(contig_id)
        reform = {(outerKey, innerKey): values for outerKey, innerDict in best_primer_pair.iteritems() for innerKey, values in innerDict.iteritems()}

        for key, value in reform.items():
            if key in res:
                res[key].append(value)
            else:
                res[key] = [value]

    print "Contigs with no primer pair found:", len(contig_ids_no_primer_found)
                
    # Format as dataframe
    res = pd.DataFrame(res, index=index)
    res.columns = [column_prefix + x + "_" + y for (x, y) in zip(res.columns.get_level_values(0), res.columns.get_level_values(1))] # Rename columns to flatten multilevel index
    column_failure_cause = column_prefix + "failure_cause"
    res[column_failure_cause] = ""
    
    # Append empty rows for contigs with a failure, add annotation of failure cause

    for contig_id in contig_ids_no_primer_found:
        # No primer pair was found
        res = res.append(pd.Series(name=contig_id))
        
    res[column_failure_cause].loc[contig_ids_no_primer_found] = "no_primer_found" + failure_cause_suffix
        
    for contig_id in contig_ids_no_C_region_match:
        # No C region match was found
        res = res.append(pd.Series(name=contig_id))
        
    res[column_failure_cause].loc[contig_ids_no_C_region_match] = "no_C_region_match" + failure_cause_suffix
    
    # Merge
    cell_contigs_primers = cell_contigs.copy(deep=True) # Copy dataframe (so we don't screw up original)
    cell_contigs_primers = cell_contigs_primers.merge(res, left_on=left_on, right_index=True)

    return cell_contigs_primers

def _find_primers_pcr1(contig_id, cbcs, contig_umis, C_region_starts, seqs_full_cdna, partial_read1, tso):
    
    # Choose global primer3 settings
    primer3_global_args_pcr1 = {
                                'PRIMER_OPT_SIZE': 26,
                                'PRIMER_MIN_SIZE': 22,
                                'PRIMER_MAX_SIZE': 35,
                                'PRIMER_OPT_TM': 67.0,
                                'PRIMER_MIN_TM': 53.0,
                                'PRIMER_MAX_TM': 72.0,
                                'PRIMER_MIN_GC': 30.0,
                                'PRIMER_MAX_GC': 70.0,
                                'PRIMER_MAX_POLY_X': 5,
                                'PRIMER_SALT_MONOVALENT': 50.0,
                                'PRIMER_SALT_DIVALENT': 2.5,
                                'PRIMER_DNA_CONC': 200.0,
                                'PRIMER_MAX_NS_ACCEPTED': 0,
                                'PRIMER_MAX_SELF_ANY': 8,
                                'PRIMER_MAX_SELF_END': 3,
                                'PRIMER_MAX_SELF_ANY_TH': 47.0,
                                'PRIMER_PAIR_MAX_COMPL_ANY': 8,
                                'PRIMER_PAIR_MAX_COMPL_END': 3,
                                'PRIMER_PRODUCT_SIZE_RANGE': [[250,300],[300,350],[350,400],[400,450],[450,500],[500,550],[550,600],[600,650],[650,700],[700,750]],
                                'PRIMER_PICK_INTERNAL_OLIGO': 0
                                }

    SEQUENCE_INCLUDED_REGION_PARTIAL_READ1 = 5 # number of bases that may be included within Partial Read 1 sequence
    
    # Get CBC, UMIs, and C region start
    cbc = cbcs[contig_id]
    umis = sorted(contig_umis[contig_id]) # sort UMIs to ensure reproducibility (when taking first UMI due to a tie for best primer score)
    C_region_start = int(C_region_starts[contig_id])
        
    # Find primer pair for each molecule
    
    primer3_parsed_results_all_umis = []
    
    for umi in umis:
        
        contig_id_umi = str(contig_id + "~" + umi)
        seq_cdna = seqs_full_cdna[contig_id_umi]
        start_target = len(partial_read1) + len(cbc) + len(umi) # position of beginning of TSO sequence
        len_target = len(tso) + C_region_start # length of entire target (TSO plus contig sequence up to C region)
        start_included_region = len(partial_read1) - SEQUENCE_INCLUDED_REGION_PARTIAL_READ1 # allow primers to include up to SEQUENCE_INCLUDED_REGION_PARTIAL_READ1 bases of Partial Read 1 sequence
        len_included_region = len(seq_cdna) - start_included_region
                
        seq_args = {'SEQUENCE_ID': contig_id_umi,
                    'SEQUENCE_TEMPLATE': seq_cdna,
                    'SEQUENCE_TARGET': [start_target, len_target],
                    'SEQUENCE_INCLUDED_REGION': [start_included_region, len_included_region]}

        primer3_res = primer3.bindings.designPrimers(seq_args, global_args=primer3_global_args_pcr1)
        # adding a mispriming library makes this the designPrimers() step take a very long time, so we proceed without explicitly considering mispriming
        
        primer3_parsed_res, primer3_parsed_notes = primer3_parser(primer3_res) # parse primer3 results (into a sensible dictionary)
        
        # add annotations to each pair
        for primer_info in primer3_parsed_res:
            amplicon_sequence = _find_amplicon_sequence(seq_cdna, primer_info) # determine amplicon sequence
            primer_info['pair']['umi'] = umi # add annotation of UMI
            primer_info['pair']['amplicon_sequence'] = amplicon_sequence # add annotation of amplicon sequence
            
#             if amplicon_sequence is not None:
#                 primer_info['pair']['len_amplicon_sequence'] = len(amplicon_sequence)
#             else:
#                 primer_info['pair']['len_amplicon_sequence'] = None
        
        primer3_parsed_results_all_umis.extend(primer3_parsed_res) # add to list of all possible primers

    # Sort primer pairs by score (penalty)
    primer3_parsed_results_all_umis = sorted(primer3_parsed_results_all_umis, key=lambda k: k['pair']['penalty'])

    # Return best pair if a pair is found
    if len(primer3_parsed_results_all_umis) > 0:
        return primer3_parsed_results_all_umis[0]
    else:
        return None

def _find_amplicon_sequence(sequence_template, primer_info):
    """ Find sequence of amplicon, given template and primers """
    
    index_start = None
    index_end = None
    
    if "left" in primer_info.keys():
        if "position" in primer_info["left"].keys():
            index_start = primer_info["left"]["position"]
    
    if "right" in primer_info.keys():
        if "position" in primer_info["right"].keys():
            index_end = primer_info["right"]["position"] + 1
            
    if index_start is not None and index_end is not None:    
        return sequence_template[index_start:index_end]
    else:
        return None
    
def primer3_parser(primer3_results):
    ''' Parse Primer3 designPrimers output, and sort it into a hierachical
    dictionary structure of primer pairs.

    This method return 2 outputs, the list of primer pairs and a dictionary with
    notes (the explanatory output from Primer3).

    Author: Martin CF Thomsen
    '''

    primer_pairs = {}
    notes = {}
    for k in primer3_results:
        if 'PRIMER_RIGHT' == k[:12]:
            key = 'right'
            tmp = k[13:].split('_', 1)
            if tmp[0].isdigit():
                id = int(tmp[0])
                if not id in primer_pairs:
                    primer_pairs[id] = {
                        'pair': {},
                        'right': {},
                        'left': {},
                        'internal': {},
                        }
                if len(tmp) > 1:
                    key2 = tmp[1].lower()
                    primer_pairs[id][key][key2] = primer3_results[k]
                else:
                    primer_pairs[id][key]['position'] = primer3_results[k][0]
                    primer_pairs[id][key]['length'] = primer3_results[k][1]
            elif tmp[0] == 'EXPLAIN':
                notes[key] = primer3_results[k]
            elif tmp == ['NUM', 'RETURNED']:
                pass
            else:
                print k
        elif 'PRIMER_LEFT' == k[:11]:
            key = 'left'
            tmp = k[12:].split('_', 1)
            if tmp[0].isdigit():
                id = int(tmp[0])
                if not id in primer_pairs:
                    primer_pairs[id] = {
                        'pair': {},
                        'right': {},
                        'left': {},
                        'internal': {},
                        }
                if len(tmp) > 1:
                    key2 = tmp[1].lower()
                    primer_pairs[id][key][key2] = primer3_results[k]
                else:
                    primer_pairs[id][key]['position'] = primer3_results[k][0]
                    primer_pairs[id][key]['length'] = primer3_results[k][1]
            elif tmp[0] == 'EXPLAIN':
                notes[key] = primer3_results[k]
            elif tmp == ['NUM', 'RETURNED']:
                pass
            else:
                print k
        elif 'PRIMER_PAIR' == k[:11]:
            key = 'pair'
            tmp = k[12:].split('_', 1)
            if tmp[0].isdigit():
                id = int(tmp[0])
                if not id in primer_pairs:
                    primer_pairs[id] = {
                        'pair': {},
                        'right': {},
                        'left': {},
                        'internal': {},
                        }
                if len(tmp) > 1:
                    key2 = tmp[1].lower()
                    primer_pairs[id][key][key2] = primer3_results[k]
                else:
                    print (k, primer3_results[k])
            elif tmp[0] == 'EXPLAIN':
                notes[key] = primer3_results[k]
            elif tmp == ['NUM', 'RETURNED']:
                pass
            else:
                print k
        elif 'PRIMER_INTERNAL' == k[:15]:
            key = 'internal'
            tmp = k[16:].split('_', 1)
            if tmp[0].isdigit():
                id = int(tmp[0])
                if not id in primer_pairs:
                    primer_pairs[id] = {
                        'pair': {},
                        'right': {},
                        'left': {},
                        'internal': {},
                        }
                if len(tmp) > 1:
                    key2 = tmp[1].lower()
                    primer_pairs[id][key][key2] = primer3_results[k]
                else:
                    primer_pairs[id][key]['position'] = primer3_results[k][0]
                    primer_pairs[id][key]['length'] = primer3_results[k][1]
            elif tmp[0] == 'EXPLAIN':
                notes['pair'] = primer3_results[k]
            elif tmp == ['NUM', 'RETURNED']:
                pass
            else:
                print (k, tmp[0])
        else:
            print k

    return (list(map(primer_pairs.get, sorted(primer_pairs.keys()))), notes)

def design_primers_contigs_pcr2_IGH_IGKL(cell_contigs_primers, contigs_seqs_all, contigs_igblast, verbose=False):
    """ Design primers for PCR2 for IGH and IGKL for all contigs in cell_contigs """
    
    if verbose:
        print "Starting PCR2 primer design"
        print
    
    # IGH
    
    if verbose:
        print "Starting PCR2 primer design for IGH..."
        
    start_time = time.time()
    
    contig_ids = cell_contigs_primers["contig_id_IGH"] # Get contig_ids
    contigs_seqs_all_focal = contigs_seqs_all.set_index("contig_name").loc[contig_ids] # get sequences
    cell_contigs_primers = _design_primers_contigs_pcr2(cell_contigs_primers, contigs_igblast, column_contig_id="contig_id_IGH", column_amplicon_sequence="IGH_pcr1_pair_amplicon_sequence", column_prefix="IGH_pcr2_", failure_cause_suffix="~IGH_pcr2", left_on="contig_id_IGH", verbose=False)

    if verbose:
        print "Design finished for PCR2"
        print "Elapsed time for PCR2 IGH primer design (s):", time.time() - start_time
        print
        
    # IGKL
    
    if verbose:
        print "Starting PCR2 primer design for IGKL..."
        
    start_time = time.time()
    
    contig_ids = cell_contigs_primers["contig_id_IGKL"] # Get contig_ids
    contigs_seqs_all_focal = contigs_seqs_all.set_index("contig_name").loc[contig_ids] # get sequences
    cell_contigs_primers = _design_primers_contigs_pcr2(cell_contigs_primers, contigs_igblast, column_contig_id="contig_id_IGKL", column_amplicon_sequence="IGKL_pcr1_pair_amplicon_sequence", column_prefix="IGKL_pcr2_", failure_cause_suffix="~IGKL_pcr2", left_on="contig_id_IGKL", verbose=False)

    if verbose:
        print "Design finished for PCR2"
        print "Elapsed time for PCR2 IGKL primer design (s):", time.time() - start_time
        print

    # Determine whether cell is successful (i.e., whether there was an error at any step)
    columns = [x for x in list(cell_contigs_primers.columns) if "failure" in x]
    cell_contigs_primers["success"] = (cell_contigs_primers[columns] == "").all(axis=1)
        
    return cell_contigs_primers

def _design_primers_contigs_pcr2(cell_contigs_primers, contigs_igblast, column_contig_id, column_amplicon_sequence, column_prefix, failure_cause_suffix, left_on, verbose=False):
    """ Design primers for contigs specified by contig_ids, merge results into dataframe cell_contigs """
    
    best_primer_pairs = {}
    
    contig_ids_failed_pcr1 = []
    
    for row_index, row in cell_contigs_primers.iterrows():
        
        contig_id = row[column_contig_id]
        
        if row[column_amplicon_sequence] is not np.nan:
            best_primer_pairs[contig_id] = _design_primer_contig_pcr2(contig_id, row[column_amplicon_sequence], contigs_igblast)
        else:
            contig_ids_failed_pcr1.append(contig_id)
                                    
    ##### Merge results into cell contig dataframe
    
    # Reform all primer results into an appropriate data frame for merging

    res = {}
    index = []
    
    contig_ids_no_primer_found = [] # keep track of contigs that had no primer pair

    for contig_id, best_primer_pair in best_primer_pairs.items():
        
        if best_primer_pair is None: 
            contig_ids_no_primer_found.append(contig_id)
            continue
            
        index.append(contig_id)
        reform = {(outerKey, innerKey): values for outerKey, innerDict in best_primer_pair.iteritems() for innerKey, values in innerDict.iteritems()}

        for key, value in reform.items():
            if key in res:
                res[key].append(value)
            else:
                res[key] = [value]

    print "Contigs with no primer pair found:", len(contig_ids_no_primer_found)
                
    # Format as dataframe
    res = pd.DataFrame(res, index=index)
    res.columns = [column_prefix + x + "_" + y for (x, y) in zip(res.columns.get_level_values(0), res.columns.get_level_values(1))] # Rename columns to flatten multilevel index
    column_failure_cause = column_prefix + "failure_cause"
    res[column_failure_cause] = ""
    
    # Append empty rows for contigs with a failure, add annotation of failure cause

    for contig_id in contig_ids_no_primer_found:
        # No primer pair was found
        res = res.append(pd.Series(name=contig_id))
    res[column_failure_cause].loc[contig_ids_no_primer_found] = "no_primer_found" + failure_cause_suffix
    
    for contig_id in contig_ids_failed_pcr1:
        # PCR 1 failed
        res = res.append(pd.Series(name=contig_id))
    res[column_failure_cause].loc[contig_ids_failed_pcr1] = "pcr1_failed" + failure_cause_suffix
    
    # Merge
    cell_contigs_primers_result = cell_contigs_primers.copy(deep=True) # Copy dataframe (so we don't screw up original)
    cell_contigs_primers_result = cell_contigs_primers_result.merge(res, left_on=left_on, right_index=True)

    return cell_contigs_primers_result


def _design_primer_contig_pcr2(contig_id, template_sequence, contigs_igblast):
    
    """ Design primer sequences that start on the V and end on the J """
    
    # Choose global primer3 settings
    primer3_global_args_pcr2 = {
                                'PRIMER_OPT_SIZE': 20,
                                'PRIMER_MIN_SIZE': 13,
                                'PRIMER_MAX_SIZE': 35,
                                'PRIMER_OPT_TM': 60.0,
                                'PRIMER_MIN_TM': 50.0,
                                'PRIMER_MAX_TM': 70.0,
                                'PRIMER_MIN_GC': 30.0,
                                'PRIMER_MAX_GC': 70.0,
                                'PRIMER_MAX_POLY_X': 5,
                                'PRIMER_SALT_MONOVALENT': 50.0,
                                'PRIMER_SALT_DIVALENT': 2.5,
                                'PRIMER_DNA_CONC': 200.0,
                                'PRIMER_MAX_NS_ACCEPTED': 0,
                                'PRIMER_MAX_SELF_ANY': 8,
                                'PRIMER_MAX_SELF_END': 3,
                                'PRIMER_MAX_SELF_ANY_TH': 47.0,
                                'PRIMER_PAIR_MAX_COMPL_ANY': 8,
                                'PRIMER_PAIR_MAX_COMPL_END': 3,
                                'PRIMER_PRODUCT_SIZE_RANGE': [[100,150], [150,200],[200,250],[250,300],[300,350],[350,400],[400,450],[450,500],[500,550],[550,600],[600,650],[650,700]],
                                'PRIMER_PICK_INTERNAL_OLIGO': 0
                                }
    
    # Find boundaries of V and J in template
    seq_VDJ = contigs_igblast.loc[contig_id]["SEQUENCE_VDJ"].replace("-", "")
    V_start = template_sequence.find(seq_VDJ)
    J_end = V_start + len(seq_VDJ)
    len_VDJ = len(seq_VDJ)
    
    start_target = V_start
    len_target = len_VDJ # length of entire target
    
    seq_args = {'SEQUENCE_ID': contig_id,
            'SEQUENCE_TEMPLATE': template_sequence,
            'SEQUENCE_INCLUDED_REGION': [start_target, len_target],
            'SEQUENCE_FORCE_LEFT_START': start_target,
            'SEQUENCE_FORCE_RIGHT_START': start_target + len_target - 1}
    
    # print seq_args
            
    primer3_res = primer3.bindings.designPrimers(seq_args, global_args=primer3_global_args_pcr2)
    primer3_parsed_res, primer3_parsed_notes = primer3_parser(primer3_res) # parse primer3 results (into a sensible dictionary)
    primer3_parsed_res_sorted = sorted(primer3_parsed_res, key=lambda k: k['pair']['penalty'])
        
    # Return best pair if a pair is found
    if len(primer3_parsed_res_sorted) > 0:
        return primer3_parsed_res_sorted[0]
    else:
        return None
