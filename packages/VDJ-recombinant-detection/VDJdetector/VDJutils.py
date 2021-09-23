import pandas as pd
import re

"""
Functions for identifying VDJ recombinants from V and J alignments

CDR3 boundaries (immunoDB annotation) for V genes (FR3--CDR3) and Jgenes (CDR3--FR4) are mapped from reference positions to read positions using these functions. a complete VDJ recombinant will pass 

"xxxC"

"""

def import_VJ_reference(reference_info_file):
    """
    Read VJ reference file as a dictionary of records

    Parameters
    ----------
    primer_info_file : path
        Path to immunoDB summary file for VJ genes        

    Returns
    -------
    reference : dict
        Dictionary of reference, organized by sequence ID
    """

    reference = {}
    genes = pd.read_csv(reference_info_file, sep=",")
    if not 'CDR3_target_pos' in genes:
        genes['CDR3_target_pos'] = None
    if not 'locus_type' in genes:
        genes['locus_type'] = ""
    if not 'IDB_type' in genes:
        genes['IDB_type'] = ""
    if not 'name' in genes:
        genes['name'] = ""

    for seqId, seq, refName, cdr3_pos, locus_type, idb_type in zip(genes['sequence_id'],
                                                                   genes["sequence"],
                                                                   genes["name"],
                                                                   genes['CDR3_target_pos'],
                                                                   genes['locus_type'],
                                                                   genes['IDB_type']):
        if not seqId in reference:                
            reference[seqId] = {"seq": seq,                                    
                                "cdr3_pos": cdr3_pos,
                                "locus_type": locus_type,
                                "IDB_type": idb_type}
    return reference


def parse_cdr3_boundaries(alignment, read):
    """
    Get the cdr3 boundary element for either the Vgene or Jgene.

    Parameters
    ----------
    alignment : dict
    dictionary of alignment information for some IG gene
    read : dict
    read strucuture dictionary, returned by FastqIterator.

    Returns 
    -------
    alignment : alignment dictionary, filled with cdr3 boundary information.
    
    Example:
    Suppose for a given Read we have the following Gene alignment

    Read = ACTGGTCAGT-TAGCTGATGCGCTA
                 |||| ||||||*|||||||
    Gene =       CAGTCTAGCTGATGCGCTA 
    
    And suppose for the Gene we have the following Annotation for the CDR3:
    (Gene = CAGTCTAGCTGA[TGC]GCTA, [TGC] = CDR3 start at position 13 of the Gene )
    
    The Goal of this function is to map the CDR3 codon position from the Gene (13) to the Read (18) using the alignment coordinates and cigar.
    This function works in three steps
       1) collect the start positions of the Gene alignment against the Read
       2) walk the cigar, updating the position of the Read to the corresponding Gene position.
       3) If the Read alignment spans the CDR3 codon start and end positions, update the alignment information with the read CDR3 position.  

    Finally, Since the CDR3 start and end positions are annotated at the codon level, the CDR3 read position is adjusted to the first base of the CDR3-start codon (V genes) and the last base of the CDR3-end codon (J genes).

   """
    # define boundaries to cdr3 pos
    # get alignment info needed to find cdr3    
    if alignment is None:
        return alignment
    if alignment['cdr3_pos'] is None:
        alignment['cdr3_read_pos'] = -1
        return alignment
    cigar_tups = cig2tups(alignment['cigar'])
    ref_pos = alignment['start'] + 1
    read_pos = alignment['query_start']
    read_end = alignment['query_end']
    cdr3_codon_start = alignment['cdr3_pos']
    cdr3_codon_end = alignment['cdr3_pos'] + 3
    # useful variables
    cdr3Start = 0
    cdr3End = 0
    # walk down the alignment, parse sequence of the cdr3 boundary.
    for cig_type, cig_len in cigar_tups:
        for i in range(cig_len):
            # boundary start
            if ref_pos == cdr3_codon_start:
                cdr3Start = read_pos
            if ref_pos == cdr3_codon_end:
                cdr3End = read_pos
            # move down alignment
            ref_pos, read_pos = walk_cigar(
                cig_type, ref_pos, read_pos)
    if (cdr3End - cdr3Start) == 3 and cdr3End > 0:
        if alignment['reference_name'][3] == 'V':                
            alignment['cdr3_read_pos'] = cdr3Start
        else:
            alignment['cdr3_read_pos'] = cdr3End
    else:
        alignment['cdr3_read_pos'] = -1
    return alignment

def fetch_cdr3(read, valn, jaln):
    """
    Given alignments with cdr3 boundaries, fetch the cdr3 sequence

    Parameters
    ----------
    read : dict
    read strucuture dicitionary, returned by FastqIterator.
    valn : dict
    alignment dictionary for the V gene
    jaln : dict
    alignment dictionary for the J gene

    Returns 
    -------
    str : string of the cdr3 sequence, if the V and J gene alignment identify boundary elements.
    """
    seq = read["seq"]
    qual = read["qual"]
    Vpos = valn['cdr3_read_pos']
    Jpos = jaln['cdr3_read_pos']
    if Vpos == -1 or Jpos == -1:
        return ("", "")
    else:
        return (seq[Vpos:Jpos], qual[Vpos:Jpos])

def cig2tups(cigar):
    """
    Convert alignment cigar string to list of tuples

    Given a Cigar string : 24M1X5M1D10M
    Parse cigar type and length, producing a list representative tuples : 
    [(0, 24), (-1, 1), (0, 5), (2, 1), (0, 10)]


    Parameters
    ----------
    cigar : str
    Alignment cigar string

    Returns
    -------
    cigar_tups : list of tuples
    """
    cigar_tups = []
    for x in re.findall('(\d+[MSXID])', cigar):
        cig_val = int(re.sub("[MSXID]", "", str(x)))
        cig_type = re.sub(str(cig_val), "", str(x))
        if cig_type == 'S':
            continue
        if cig_type == 'M':
            cigar_tups.append((0, cig_val))
        if cig_type == 'X':
            cigar_tups.append((-1, cig_val))
        if cig_type == 'I':
            cigar_tups.append((1, cig_val))
        if cig_type == 'D':
            cigar_tups.append((2, cig_val))
    return cigar_tups

def walk_cigar(cig_type, query_pos, read_pos):
    """
    Given Alignment Coordinates for query and subject, determine the 
    alignment of each positive relative to the cigar type.  

    Parameters
    ----------
    cig_type : int
    integer of the cigar type, indicating if the alignment is an indel or on the same position
    query_pos : int
    Current position in the query sequence
    read_pos : int
    Current position in the read sequence
    errors : int
    Number of errors tranversed, given the cig_type. 

    Returns
    -------
    tuple : query position, read position, and adjuested error calculation for moving one position in the alignment. 
    """
    if cig_type == 0:  # mapped
        query_pos += 1
        read_pos += 1
    elif cig_type == -1:  # mismatch
        query_pos += 1
        read_pos += 1
    elif cig_type == 1:  # insertion
        read_pos += 1
    elif cig_type == 2:  # deletion
        query_pos += 1
    return (query_pos, read_pos)

def search_J_motif(aln_pair):
        """
        for read with an identified vgene, but no jgene, look for the J cdr3 boundary motif: [FW]G[A-Z]G
        
        Parameters
        ----------
        aln_pair : dict
        read alignment information, returned by BamStreamer and CDR3/alignment classified
        
        Returns
        -------
        (aln_pair, jalnType) : tuple
        tuple of aln_pair object and alignment classification for J
        """
        jalnType = "missing"
        cdr3Start = aln_pair["vgene"]["cdr3_read_pos"]
        ##look for [FW]G[A-Z]G motif
        motif = "(TT[TC]|TGG)GG....GG."
        found = 0
        cdr3End = -1
        jStats = {}
        for x in re.finditer(motif, aln_pair["seq"]):            
            if x.start() > cdr3Start:
                found += 1
                cdr3End = int(x.start()) + 3
                jStats["cdr3_read_pos"] = cdr3End
                break
        if found > 0:            
            cdr3 = fetch_cdr3(aln_pair, aln_pair["vgene"], jStats)
            aln_pair["cdr3"] = cdr3[0]
            aln_pair["cdr3_qual"] = cdr3[1]
            jalnType = "motif_found"
        return (aln_pair, jalnType)
