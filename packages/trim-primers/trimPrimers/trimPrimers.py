import pandas
import os
import random
import argparse
import re
from FastqStreamer.FastqReader import FastqReader
from BamStreamer.BamPairs import read_two_bams
from trimPrimers.utils import *

class IPetePrimerTrimmer(FastqReader):
    """
    Parse primer alignments from reads, trim primers from ends
    """

    def __init__(self, primer_info, expected_v_start):
        """
        Parameters
        ----------
        reference_annot_file : Str
        reference annotation format for primer sequences

        """
        FastqReader.__init__(self)
        self.expected_v_start = expected_v_start
        self.reference = import_primer_reference(primer_info)
        
    def get_aln_stats(self, read, whichGene):
        """
        get alignment information about a read

        Parameters
        ----------
        read : dict
        read record returned by FastqReader

        whichGene : Str
        which read "read1" or "read2"

        Returns
        -------
        stats : Str
        A string containing alignment information
        """
        
        stats = "{name}\t{start}\t{end}\t{qStart}\t{qEnd}\t{score}\t{cigar}\t{mismatch}\t{insertion}\t{deletion}\t{percentId}\t{alnType}\t{primerLen}"
        if read[whichGene] is not None:
            refName = read[whichGene]["reference_name"]
            return stats.format(name=read[whichGene]["reference_name"].split("__")[0],
                                start=read[whichGene]["start"],
                                end=read[whichGene]["end"],
                                qStart=read[whichGene]["query_start"],
                                qEnd=read[whichGene]["query_end"],
                                score=read[whichGene]["score"],
                                cigar=read[whichGene]["cigar"],
                                mismatch=read[whichGene]["mismatches"],
                                insertion=read[whichGene]["insertions"],
                                deletion=read[whichGene]["deletions"],
                                percentId=read[whichGene]["percent_identity"],
                                alnType=read[whichGene]["classification"],
                                primerLen=len(self.reference[refName]["seq"]))
        else:
            return stats.format(name="",
                                start="",
                                end="",
                                qStart="",
                                qEnd="",
                                score="",
                                cigar="",
                                mismatch="",
                                insertion="",
                                deletion="",
                                percentId="",
                                alnType="",
                                primerLen="")

    def write_aln_stats(self, read, filename):
        """
        write alignment information to file, combining V and J gene alignments

        Parameters
        ----------
        read : dict
        read record returned by FastqReader

        filename: Str
        The name of the stats file to be written

        Returns
        -------
        None
        """
        fh = self.file_handles_[filename]
        stats = "{name}\t{vstats}\t{jstats}\t{avg_read_qual}\n"
        vstats = self.get_aln_stats(read, "vgene")
        jstats = self.get_aln_stats(read, "jgene")
        fh.write(
            stats.format(name=read["read_name"],
                         vstats=vstats,
                         jstats=jstats,
                         avg_read_qual=read["avgQual"]
                         )
        )

    def create_stats_handle(self, filename):
        """
        Create a file handle for alignment stats

        Parameters
        ----------
        filename : Str
        The name of the stats file to be written

        Returns
        -------
        None
        """
        header = "name\tv_gene\tv_start\tv_end\tv_qStart\tv_qEnd\tv_score\tv_cigar\tv_mismatch\tv_insertion\tv_deletion\tv_percent_id\tv_alnType\tv_primer_len\tj_gene\tj_start\tj_end\tj_qStart\tj_qEnd\tj_score\tj_cigar\tj_mismatch\tj_insertion\tj_deletion\tj_percent_id\tj_alnType\tj_primer_len\tavg_read_qual\n"
        if not filename in self.file_handles_:
            fh = open(filename, 'w')
            self.file_handles_[filename] = fh
            fh.write(header)
        fh = self.file_handles_[filename]
        return filename

    def trim_alignments(self, read, whichGene):
        """
        get alignment information about a read

        Parameters
        ----------
        read : dict
        read record returned by FastqReader

        whichGene : Str
        which read "read1" or "read2"

        Returns
        -------
        classify : Str
        alignment classification string containing alignment information
        """
        
        classify = ""
        seq = read["seq"]        
        start = read[whichGene]["query_start"]
        end = read[whichGene]["query_end"]        
        mismatch = read[whichGene]["mismatches"]
        insertion = read[whichGene]["insertions"]
        deletion = read[whichGene]["deletions"]        
        strand = read[whichGene]["strand"]
        refLen = len(self.reference[read[whichGene]["reference_name"]]["seq"])
        clipped = refLen - (end - start + 1)
        indels = insertion + deletion
        totDiff = indels + mismatch + clipped
        trimSeq = read["trim_seq"]
        trimQual = read["trim_qual"]    
        ## check V primer alignment
        if whichGene == "vgene":
            if mismatch <= 2 and indels <=1 and clipped <=1 and totDiff <=3:     
                trimSeq = trimSeq[end:]
                trimQual = trimQual[end:]
                classify = "good_quality"
            else:
                classify = "low_quality"
        ## check J primer alignment
        if whichGene == "jgene":
            if mismatch <= 2 and indels <=1 and clipped <=1 and totDiff <=3:        
                classify = "good_quality"
                toTrim = (len(seq) - start) + 1 - 5 ## +1 for zero-based, and - 5 to leave CDR3 boundary
                trimSeq = trimSeq[:-toTrim]
                trimQual = trimQual[:-toTrim]
            else:
                classify = "low_quality"
        read[whichGene]["classification"] = classify    
        read["trim_seq"] = trimSeq
        read["trim_qual"] = trimQual                
        return read

    def trim_primers(self, vbam, jbam, basename):
        """
        Identify primer alignments, trim alignments from reads
        
        Parameters
        ----------
        Vbam : string
        Read alignments against V gene primers
        Jbam : string
        Read alignments against J gene primers
        basename : string
        basename for report files

        Returns
        -------
        None
        """
        #######################
        ## create file handles
        #######################
        trim_reads = self.create_fastq_handle("{}_trim_primers.fastq".format(basename))
        too_short = self.create_fastq_handle("{}_too_short.fastq".format(basename))
        primer_hits = self.create_stats_handle("{}_trim_primers.tsv".format(basename))
        short_hits = self.create_stats_handle("{}_too_short.tsv".format(basename))
        #############################
        # filter reads by alignment
        #############################        
        for read in read_two_bams(vbam, "vgene", jbam, "jgene"):            
            read["trim_seq"] = read["seq"]
            read["trim_qual"] = read["qual"]            
            # both V and J align
            if read["vgene"] is not None and read["jgene"] is not None:
                read = self.trim_alignments(read, "vgene")
                read = self.trim_alignments(read, "jgene")                
            elif read["vgene"] is not None and read["jgene"] is None:
                read = self.trim_alignments(read, "vgene")
            elif read["vgene"] is None and read["jgene"] is not None:
                read = self.trim_alignments(read, "jgene")
            if len(read["trim_seq"]) > 50:
                read["seq"] = read["trim_seq"]
                read["qual"] = read["trim_qual"]
                self.write_fastq(read, trim_reads)
                self.write_aln_stats(read, primer_hits)
            else:
                self.write_fastq(read, too_short)
                self.write_aln_stats(read, short_hits)
        self.close_file_handles()
        

def main():
    parser = argparse.ArgumentParser(
        description="gather amplicon stats from bam file")
    parser.add_argument("-p", "--primer_info",
                        required=True,
                        help="primer information file")
    parser.add_argument("-b", "--basename",
                        required=True,
                        help="basename for output files")
    parser.add_argument("-v", "--v_aln",
                        required=True,
                        help="V probe alignment bam file")
    parser.add_argument("-j", "--j_aln",
                        required=True,
                        help="J probe alignment bam file")
    args = parser.parse_args()
    trimmer = IPetePrimerTrimmer(args.primer_info, 14)
    trimmer.trim_primers(args.v_aln, args.j_aln, args.basename)

if __name__ == '__main__':    
    main()


