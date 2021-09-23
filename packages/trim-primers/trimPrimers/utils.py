import pandas as pd
def import_primer_reference(reference_info_file):
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
