import pandas
import argparse
import pysam
from FastqStreamer.FastqReader import FastqReader
from BamStreamer.BamUtils import read_alignments

class IPeteSpikeIn(FastqReader):
    """
    Filter or count spike-in reads from bam files

    Attributes
    ----------
    spikes : pandas data frame
    Pandas data frame containing spike-in alignments

    spike_stats : dict
    Dictionary for counting spike-in alignments
    """

    def __init__(self, spikein, reference_annot_file, basename):
        FastqReader.__init__(self)
        self.spikes = self.load_spike_aln(spikein)
        self.spike_stats = {}
        self.basename = basename

    def load_spike_aln(self, spikeBam):
        """
        Read spike in alignment names into memory, store read and spikeID together        
        
        Parameters
        ----------
        spikeBam : str
        path to spikein alignment bam file
        
        Returns
        -------
        spikes : dict
        dictionary of read name keys and spikeId values

        """
        spikes = {}
        for read in read_alignments(spikeBam):
            read_name = read["query_name"]
            spikes[read_name] = read['reference_name']
        return spikes
    
    def split_bams(self, Vbam, Jbam):
        """
        use spikein read alignments to split V and J alignments into spike-in and native read files.
        
        Sbam : str
        spikein alignments bam file
        Vbam : str
        Vgene alignments bam file
        Jbam : str
        Jgene alignments bam file
        """

        Vsam = pysam.AlignmentFile(Vbam, "rb", check_sq=False)
        nativeV = pysam.AlignmentFile("{}_Vgene_native.bam".format(self.basename), "wb", template=Vsam)        
        spikeV = pysam.AlignmentFile("{}_Vgene_spikein.bam".format(self.basename), "wb", template=Vsam)        
        Vsam.close()
        ##split V gene alignments reads based on spikein alignments
        for read in read_alignments(Vbam):
            if read["query_name"] in self.spikes:
                spikeV.write(read['pyread'])
            else:
                nativeV.write(read['pyread'])                
        nativeV.close()
        spikeV.close()                
        Jsam = pysam.AlignmentFile(Jbam, "rb", check_sq=False)
        nativeJ = pysam.AlignmentFile("{}_Jgene_native.bam".format(self.basename), "wb", template=Jsam)        
        spikeJ = pysam.AlignmentFile("{}_Jgene_spikein.bam".format(self.basename), "wb", template=Jsam)        
        Jsam.close()
        ##split J gene alignments reads based on spikein alignments
        for read in read_alignments(Jbam):
            if read["query_name"] in self.spikes:
                spikeJ.write(read['pyread'])
            else:
                nativeJ.write(read['pyread'])                
        nativeJ.close()
        spikeJ.close()

        ###################################################################################
        ##can also process with BamLinker, but is slower with larger files, because of I/O
        ##this is bit faster with smaller files though, there are trade-offs.
        ################################################################################### 
        # ##prepare file handles for output
        # Vsam = pysam.AlignmentFile(Vbam, "rb", check_sq=False)
        # nativeV = pysam.AlignmentFile(Vbam.replace(".bam", "_native.bam"), "wb", template=Vsam)        
        # spikeV = pysam.AlignmentFile(Vbam.replace(".bam", "_spikein.bam"), "wb", template=Vsam)        
        # Vsam.close()
        # ##split V gene alignments reads based on spikein alignments
        # for read in self.read_two_bams(Vbam, "vgene", Sbam, "spikein"):
        #     if read["spikein"] is not None:                
        #         if read["vgene"] is not None:                    
        #             spikeV.write(read['vgene']['pyread'])
        #     elif read["spikein"] is None:
        #         if read["vgene"] is not None:
        #             nativeV.write(read['vgene']['pyread'])
        # nativeV.close()
        # spikeV.close()
        # ##split J gene alignments reads based on spikein alignments
        # Jsam = pysam.AlignmentFile(Jbam, "rb", check_sq=False)
        # nativeJ = pysam.AlignmentFile(Jbam.replace(".bam", "_native.bam"), "wb", template=Jsam)
        # spikeJ = pysam.AlignmentFile(Jbam.replace(".bam", "_spikein.bam"), "wb", template=Jsam)
        # Jsam.close()
        # for read in self.read_two_bams(Jbam, "jgene", Sbam, "spikein"):
        #     if read["spikein"] is not None:                
        #         if read["jgene"] is not None:
        #             spikeV.write(read['jgene']['pyread'])
        #     elif read["spikein"] is None:
        #         if read["jgene"] is not None:
        #             nativeV.write(read['jgene']['pyread'])
        # nativeJ.close()        
        # spikeJ.close()
        # Jsam.close()
        
    def report_spikein_counts(self, VJaln, outFile):
        """
        write a report file of spike-in sequences

        Parameters
        ----------
        outFile : str
        Path to output file

        Returns
        -------
        None
        """
        # summarize VJ alignments of spike-in reads
        VJspike = {}
        try:            
            aln = pandas.read_csv(VJaln, sep="\t")            
            if len(aln) > 0:
                for name, Vgene, Jgene, uid_group in zip(aln["name"], aln["v_gene"], aln["j_gene"], aln["UMI_group"]):
                    VJpair = "{}__{}".format(Vgene, Jgene)
                    if name in self.spikes:
                        spikeId = self.spikes[name]
                        if not spikeId in VJspike:
                            VJspike[spikeId] = {"match": 0, "diff": 0, "uids": []}
                        if VJpair == spikeId:
                            VJspike[spikeId]["match"] += 1
                        else:
                            VJspike[spikeId]["diff"] += 1
                        uid = "{}_{}".format(VJpair, uid_group)
                        VJspike[spikeId]["uids"].append(uid)
                # create summary for spike-in reads
                for name in self.spikes:
                    spike = self.spikes[name]
                    if not spike in self.spike_stats:
                        self.spike_stats[spike] = {"reads": 0,
                                                   "on_target": 0,
                                                   "VJ_match": 0,
                                                   "dedup_count": 0}
                    self.spike_stats[spike]["reads"] += 1
                for spike in self.spike_stats:
                    if spike in VJspike:
                        self.spike_stats[spike]["on_target"] = (
                            VJspike[spike]["match"] + VJspike[spike]["diff"])
                        self.spike_stats[spike]["VJ_match"] = VJspike[spike]["match"]
                        self.spike_stats[spike]["dedup_count"] = len(set(VJspike[spike]["uids"]))
        except Exception as e:
            print(e)        
            pass
        # spike-in report
        with open(outFile, 'w') as fhOut:
            fhOut.write("synthetic,reads,on_target,VJ_match,dedup_count\n")
            for spike in self.spike_stats:
                vals = [spike,
                        str(self.spike_stats[spike]["reads"]),
                        str(self.spike_stats[spike]["on_target"]),
                        str(self.spike_stats[spike]["VJ_match"]),
                        str(self.spike_stats[spike]["dedup_count"])]
                vals = ",".join(vals)+"\n"
                fhOut.write(vals)

def main():
    parser = argparse.ArgumentParser(
        description="Split spikein reads from native reads or count spikein molecules from a sample")
    parser.add_argument("-s", "--spike_aln",
                        required=True,
                        help="Spike in alignment bam")
    parser.add_argument("-v", "--v_aln",
                        required=False,
                        help="V gene alignment bam")
    parser.add_argument("-j", "--j_aln",
                        required=False,
                        help="J gene alignment bam")
    parser.add_argument("-VJ", "--VJaln",
                        required=False,
                        help="V & J gene alignment summary file")
    parser.add_argument("-a", "--action",
                        required=True,
                        default="split",
                        choices=["split", "count"],
                        help="fastq file")
    parser.add_argument("-r", "--reference_annot",
                        required=True,
                        help="reference annotation file")
    parser.add_argument("-b", "--basename",
                        required=True,
                        help="basename argument for report file name")
    args = parser.parse_args()
    if args.action == "count":
        if args.VJaln is None:
            parser.error("VJ alignments are required is action = count")

    spikeinParser = IPeteSpikeIn(args.spike_aln, args.reference_annot, args.basename)
    if args.action == "split":
        spikeinParser.split_bams(args.v_aln, args.j_aln)
    elif args.action == "count":
        spikeinParser.report_spikein_counts(args.VJaln, "{}_{}".format(
            args.basename, "spikein_stats.csv"))


if __name__ == '__main__':
    main()
