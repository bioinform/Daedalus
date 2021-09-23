import pandas as pd
import argparse
import re
from FastqStreamer.FastqReader import FastqReader
from BamStreamer.BamPairs import read_two_bams
from VDJdetector.VDJutils import *

class VDJdetector(FastqReader):
    """
    Identify VDJ recombinants from V and J gene alignment files.
    Parses the CDR3 region of T-cell and Immunoglobulin gene segments: V and J genes.
    
    Attributes
    ----------
    min_identity : float
    the minimum percent identity for an alignment to be considered good.
    classifications : dict
    classification counts for read alignments
    """

    def __init__(self, reference_annot_file, identity):
        """
        Parameters
        ----------
        reference_annot_file : Str
        reference annotation file for VJ genes sequences
        """

        self.reference = import_VJ_reference(reference_annot_file)

        FastqReader.__init__(self)
        self.min_identity = identity
        self.classifications = {"on_target": 0,                                
                                "artifact": 0}

    def annotate_alignments(self, alignment):
        """
        """
        seq = ''
        cdr3_pos = None
        locus_type = None
        IDB_type = None
        seqLen = 0
        if alignment is None:
            return alignment
        refName = alignment["reference_name"]
        if refName in self.reference:
            seq = self.reference[refName]['seq']
            cdr3_pos = self.reference[refName]['cdr3_pos']
            locus_type = self.reference[refName]['locus_type']
            IDB_type = self.reference[refName]['IDB_type']
            seqLen = len(seq)
        alignment['ref_length'] =  seqLen
        alignment['cdr3_pos'] = cdr3_pos
        alignment['locus_type'] = locus_type
        alignment['IDB_type'] = IDB_type
        alignment['ref_seq'] = seq        
        return alignment
    
        
    def get_aln_stats(self, read, whichGene):
        """
        get alignment information about a read

        Parameters
        ----------
        read : dict
        read record returned by FastqIterator
        
        whichGene : Str
        which gene "vgene" or "jgene"

        Returns
        -------
        stats : Str
        A string containing alignment information
        """
        
        stats = "{name}\t{start}\t{end}\t{qStart}\t{qEnd}\t{score}\t{cigar}\t{percentId}\t{targetLen}\t{cdr3_pos}\t{locus_type}\t{IDB_type}"
        if read[whichGene] is not None:
            return stats.format(name=read[whichGene]["reference_name"].split("__")[0],
                                start=read[whichGene]["start"],
                                end=read[whichGene]["end"],
                                qStart=read[whichGene]["query_start"],
                                qEnd=read[whichGene]["query_end"],
                                score=read[whichGene]["score"],
                                cigar=read[whichGene]["cigar"],
                                percentId=read[whichGene]["percent_identity"],
                                targetLen=len(read[whichGene]["ref_seq"]),
                                cdr3_pos=read[whichGene]["cdr3_read_pos"],
                                locus_type=read[whichGene]["locus_type"],
                                IDB_type=read[whichGene]["IDB_type"])
        else:
            return stats.format(name="",
                                start="",
                                end="",
                                qStart="",
                                qEnd="",
                                score="",
                                cigar="",
                                percentId="",
                                targetLen="",
                                cdr3_pos="",
                                locus_type="",
                                IDB_type="")

    def write_aln_stats(self, read, filename):
        """
        write alignment information to file, combining V and J gene alignments

        Parameters
        ----------
        read : dict
        read record returned by FastqIterator

        filename: Str
        The name of the stats file to be written

        Returns
        -------
        None
        """
        fh = self.file_handles_[filename]
        stats = "{name}\t{vstats}\t{jstats}\t{Dregion}\t{cdr3}\t{cdr3_qual}\t{avg_read_qual}\n"
        vstats = self.get_aln_stats(read, "vgene")
        jstats = self.get_aln_stats(read, "jgene")
        fh.write(
            stats.format(name=read["read_name"],
                         vstats=vstats,
                         jstats=jstats,
                         Dregion=read["Dregion"],
                         cdr3=read["cdr3"],
                         cdr3_qual=read["cdr3_qual"],
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
        filename : Str
        the filename, which serves as a key to the filehandle
        """
        header = "name\tv_gene\tv_start\tv_end\tv_qStart\tv_qEnd\tv_score\tv_cigar\tv_percent_id\tv_target_len\tv_cdr3_pos\tv_locus_type\tv_idb_type\tj_gene\tj_start\tj_end\tj_qStart\tj_qEnd\tj_score\tj_cigar\tj_percent_id\tj_target_len\tj_cdr3_pos\tj_locus_type\tj_idb_type\td_seq\tcdr3\tcdr3_qual\tavg_read_qual\n"
        if not filename in self.file_handles_:
            fh = open(filename, 'w')
            self.file_handles_[filename] = fh
            fh.write(header)
        fh = self.file_handles_[filename]
        return filename

    def classify_alignment(self, read, whichGene):
        """
        get alignment information about a read

        Parameters
        ----------
        read : dict
        read record from BamStreamer

        whichGene : Str
        which gene to inspect "vgene" or "jgene"

        Returns
        -------
        classify : Str
        alignment classification string containing alignment information
        """
        
        classify = ""

        refStart = read[whichGene]["start"]
        refEnd = read[whichGene]["end"]        
        cdr3_pos = read[whichGene]["cdr3_pos"]
        has_cdr3 = False

        ##cdr3 position points to start of codon.
        ##on target must have the complete codon
        if cdr3_pos is not None:
            if whichGene == "vgene" and refEnd >= (cdr3_pos + 2) and refStart <= (cdr3_pos - 1):
                has_cdr3 = True
            if whichGene == "jgene" and refStart <= (cdr3_pos - 1) and refEnd >= (cdr3_pos + 1):
                has_cdr3 = True

        goodAln = False
        if read[whichGene]["percent_identity"] >= self.min_identity:
            goodAln = True        
        
        # both the primer and target region must align to be considered on-target
        # a minimum overlap of 10bases is considered a hit.
        if has_cdr3 and goodAln:
            classify = "on_target"
        else:
            classify = "artifact"
        return classify

    def parse_alignments(self, vbam, jbam, basename):
        """
        classify alignments against reads, write stats, and separate reads into categories.

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
        ######################
        # create file handles
        ######################
        on_target = self.create_fastq_handle("{}_on_target.fastq".format(basename))
        on_target_hits = self.create_stats_handle("{}_on_target.tsv".format(basename))
        artifacts = self.create_fastq_handle("{}_artifacts.fastq".format(basename))
        artifacts_hits = self.create_stats_handle("{}_artifacts.tsv".format(basename))
        ##no_aln = self.create_fastq_handle("{}_no_aln.fastq".format(basename))
        #############################
        # filter reads by alignment
        #############################
        for aln_pair in read_two_bams(vbam, "vgene", jbam, "jgene"):            
            ##fetch CDR3 annotation
            aln_pair["vgene"] = self.annotate_alignments(aln_pair["vgene"])
            aln_pair["jgene"] = self.annotate_alignments(aln_pair["jgene"])
            ## prep read info
            aln_pair["vgene"] = parse_cdr3_boundaries(aln_pair["vgene"], aln_pair)
            aln_pair["jgene"] = parse_cdr3_boundaries(aln_pair["jgene"], aln_pair)
            read_name = aln_pair["read_name"]
            vals = read_name.split(" ")[0].split(":")
            aln_pair["Dregion"] = ""
            aln_pair["cdr3"] = ""
            aln_pair["cdr3_qual"] = ""
            ##if both genes aligned
            if aln_pair["vgene"] is not None and aln_pair["jgene"] is not None:            
                valnType = self.classify_alignment(aln_pair, "vgene")
                jalnType = self.classify_alignment(aln_pair, "jgene")
                ##compare strands
                vStrand = aln_pair["vgene"]["strand"]
                jStrand = aln_pair["jgene"]["strand"]
                ## fetch cdr3 stats
                cdr3End = aln_pair["jgene"]["cdr3_read_pos"]
                cdr3Start = aln_pair["vgene"]["cdr3_read_pos"]
                ## query for read position
                DStart = aln_pair["vgene"]["query_end"]
                DEnd = aln_pair["jgene"]["query_start"]                
                if DStart < DEnd:
                    aln_pair["Dregion"] = aln_pair["seq"][DStart:(DEnd - 1)]
                else:
                    aln_pair["Dregion"] = ""
                ##on target reads must have a CDR3                
                if not cdr3Start < cdr3End:                    
                    valnType = "artifact"
                    jalnType = "artifact"
                ##genes must align to the same strand
                if not vStrand == jStrand:
                    valnType = "artifact"
                    jalnType = "artifact"                    
                # if both V and J are on-target, then the fragment is on-target
                if valnType == "on_target" and jalnType == "on_target":
                    cdr3 = fetch_cdr3(aln_pair, aln_pair["vgene"], aln_pair["jgene"])
                    aln_pair["cdr3"] = cdr3[0]
                    aln_pair["cdr3_qual"] = cdr3[1]
                    self.classifications["on_target"] += 1
                    self.write_fastq(aln_pair, on_target)
                    self.write_aln_stats(aln_pair, on_target_hits)
                elif valnType == "on_target" and jalnType == "artifact":                                        
                    ##search for motif in case of J mis-alignment
                    aln_pair, jalnType = search_J_motif(aln_pair)
                    if jalnType == "motif_found":
                        self.classifications["on_target"] += 1
                        self.write_fastq(aln_pair, on_target)
                        self.write_aln_stats(aln_pair, on_target_hits)
                    else:
                        self.classifications["artifact"] += 1
                        self.write_fastq(aln_pair, artifacts)
                        self.write_aln_stats(aln_pair, artifacts_hits)
                else:
                    ##otherwise artifacts
                    self.classifications["artifact"] += 1
                    self.write_fastq(aln_pair, artifacts)
                    self.write_aln_stats(aln_pair, artifacts_hits)
            elif aln_pair["vgene"] is not None and aln_pair["jgene"] is None:
                valnType = self.classify_alignment(aln_pair, "vgene")
                jalnType = "missing"
                if valnType == "on_target":
                    aln_pair, jalnType = search_J_motif(aln_pair)
                if jalnType == "motif_found":
                    self.classifications["on_target"] += 1
                    self.write_fastq(aln_pair, on_target)
                    self.write_aln_stats(aln_pair, on_target_hits)
                else:
                    self.classifications["artifact"] += 1
                    self.write_fastq(aln_pair, artifacts)
                    self.write_aln_stats(aln_pair, artifacts_hits)
            else:
                self.classifications["artifact"] += 1
                self.write_fastq(aln_pair, artifacts)
                self.write_aln_stats(aln_pair, artifacts_hits)
        self.close_file_handles()    
        
    def summarize(self, basename):
        """
        write summary of alignment categories to file

        Parameters
        ----------
        basename : Str
        basename of output file

        Returns
        -------
        None
        """
        keys = ["on_target", "artifact"]
        vals = []
        for key in keys:
            vals.append(str(self.classifications[key]))
        reportFile = "{}_alignment_summary.csv".format(basename)
        with open(reportFile, 'w') as Fh:
            Fh.write(",".join(keys)+"\n")
            Fh.write(",".join(vals)+"\n")


def main():
    parser = argparse.ArgumentParser(
        description="Identify VDJ recombinants from V and J gene alignments")
    parser.add_argument("-r", "--reference_info",
                        required=True,
                        help="immunoDB reference annotation")
    parser.add_argument("-b", "--basename",
                        required=True,
                        help="basename for output files")
    parser.add_argument("-v", "--vgene_aln",
                        required=True,
                        help="V gene alignment bam file")
    parser.add_argument("-j", "--jgene_aln",
                        required=True,
                        help="J gene alignment bam file")
    parser.add_argument("-i", "--identity",
                        required=False,
                        default=85,
                        type=int,
                        help="minimum gene identity to be considered a good alignment")
    args = parser.parse_args()
    parser = VDJdetector(args.reference_info, args.identity)
    parser.parse_alignments(args.vgene_aln,
                            args.jgene_aln,
                            args.basename)
    parser.summarize(args.basename)

if __name__ == '__main__':
    main()
