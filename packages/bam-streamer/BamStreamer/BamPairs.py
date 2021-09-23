import pysam
from BamStreamer.BamUtils import read_alignments
from BamStreamer.BamUtils import cmp


"""
Stream Two bam file alignments by read name. Useful for reading alignments for reads aligned against two references (e.g. Gene segements, spikeins, etc...)

"""


def set_seq_info(alnPair, alignment):
    """
    Parse read stats from bam record, record high level stats as new dict values

    Parameters
    ----------
    alnPair : dict
    dictionary of bam Pair info from read_two_bams()
    alignment : dict
    bam file record, returned by BamParser

    Returns
    -------
    alnPair : dict
    dictionary with new fields, summarizing read information from bam format
    """
    alnPair["read_name"] = alignment["query_name"]
    alnPair["seq"] = alignment["pyread"].query_sequence
    alnPair["quality_values"] = alignment["pyread"].query_qualities
    alnPair["qual"] = alignment["pyread"].qual
    alnPair["avgQual"] = sum(alnPair["quality_values"])/len(alnPair["quality_values"])
    return alnPair

def read_two_bams(bam1, bam1Name, bam2, bam2Name):
    """
    Read two bam files, combine and return read alignment information from both files.
    Keeping track of read names and read names order, ensure pairs are returned in proper order.

    Parameters
    ----------
    bam1 : str
    path to first bam file
    bam1Name : str
    The name to store for the first bam File
    bam2 : str
    path to second bam file
    bam2Name : str
    The name to store for the second bam File
    """
    reader1 = read_alignments(bam1)
    reader2 = read_alignments(bam2)        
    r1 = next(reader1, None)
    r2 = next(reader2, None)
    while True:
        alnPair = {bam1Name:None, bam2Name:None}
        if r1 is None and r2 is None:
            ##break if no more alignments
            break
        elif r1 is None:
            ##keep reading r2 if r1 is done
            alnPair[bam2Name]=r2
            alnPair = set_seq_info(alnPair, r2)
            r2 = next(reader2, None)
        elif r2 is None:
            ##keep reading r1 if r2 is done
            alnPair[bam1Name]=r1
            alnPair = set_seq_info(alnPair, r1)
            r1 = next(reader1, None)
        else:
            ##both alignments exist and are the same
            if r1["query_name"] == r2["query_name"]:
                alnPair[bam1Name]=r1
                alnPair[bam2Name]=r2
                alnPair = set_seq_info(alnPair, r1)
                r1 = next(reader1, None)
                r2 = next(reader2, None)
            ##r1 alignments are behind r2, keep reading r1
            elif cmp(r1["query_name"], r2["query_name"]) < 0:
                alnPair[bam1Name]=r1
                alnPair = set_seq_info(alnPair, r1)
                r1 = next(reader1, None)
            ##r2 alignments are behind r1, keep reading r2
            elif cmp(r1["query_name"], r2["query_name"]) > 0:
                alnPair[bam2Name]=r2
                alnPair = set_seq_info(alnPair, r2)
                r2 = next(reader2, None)
        yield alnPair
                    
