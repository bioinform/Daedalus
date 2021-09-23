import pysam

def read_alignments(fname):
    """
    Read Primary alignments from Bam File, calculating high level stats from pysam alignment object.

    Parameters
    ----------
    fname : path
        Path to swifr aligned BAM/SAM file.

    Yields
    ------
    alignment : dict
        Dictionary containing the alignment information of the read.
    """
    with pysam.AlignmentFile(fname, check_sq=False) as sam:
        prev_query_name = None
        for read in sam.fetch(until_eof=True):
            if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:
                assert prev_query_name is None or cmp(prev_query_name, read.query_name) <= 0, 'File is not sorted by query name.'
                prev_query_name = read.query_name
                strand = "+"
                if read.is_reverse:
                    strand = "-"                    
                cigar_stats = read.get_cigar_stats()[0]
                alignment = {
                    'query_name': read.query_name,
                    'reference_name': read.reference_name,
                    'start': read.reference_start + 1,
                    'end': read.reference_end,
                    'query_start': read.query_alignment_start + 1,
                    'query_end': read.query_alignment_end,
                    'score': read.get_tag('AS'),
                    'cigar': read.cigarstring,
                    'percent_identity': cigar_stats[0] * 100 / sum([cigar_stats[0], cigar_stats[1], cigar_stats[2], cigar_stats[8]]),
                    'mismatches': cigar_stats[8],
                    'insertions': cigar_stats[1],
                    'deletions': cigar_stats[2],
                    'clipped': cigar_stats[4],
                    'strand' : strand,
                    'pyread' : read
                }

                yield alignment

def cmp(q1, q2):
    """Compare read names.

    Fields in read names (separated by ':') are compared as string or integer depends on their types.

    Parameters
    ----------
    q1 : str
        read name
    q2 : str
        read name

    Returns
    -------

    int
        0 if q1 == q2; -1 if q1 < q2; 1 if q1 > q2.
    """
    if q1 == q2:
        return 0
    q1 = q1.split(':')
    q2 = q2.split(':')
    for a, b in zip(q1, q2):
        if a.isdigit():
            a = int(a)
        if b.isdigit():
            b = int(b)
        if a < b:
            return -1
        elif a > b:
            return 1
        else:
            continue
