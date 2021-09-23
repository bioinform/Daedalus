import pandas
import os
import random
import argparse
import pysam
from FastqStreamer.FastqReader import FastqReader
from BamStreamer.BamUtils import read_alignments
##from FastqIterator import FastqIterator
##from link_bams import BamLinker


class IPeteSpikeIn(FastqReader):
    """
    Filter spike-in reads from fastq file

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
        
        ##can also process with `read_two_bams` via BamParser, but that is slower with larger files.

    def label_spikeins(self, spikeDF, spikeUMI_df):
        UMI_groups = spikeDF.UMI_group.unique()        
        spikeGroups = {}
        for group in UMI_groups:
            groupDF = spikeDF.loc[spikeDF["UMI_group"] == group]            
            family_size = groupDF.shape[0]
            Spikes = groupDF["spikein"].value_counts(dropna=False)
            spikeGroups[group] = Spikes.idxmax()            
        spikeUMI_df["spikein"] = spikeUMI_df["UMI_group"].map(lambda x : spikeGroups[x]) 
        return spikeUMI_df
            
    def cummulative_stats(self, df):
        df["cummulative_reads"] = df["UMI_family_size"].cumsum()
        df["cummulative_fraction"] = df["cummulative_reads"] / sum(df["UMI_family_size"])
        return df
    
    def split_umi_reports(self, VDJfile, UMIfams):            
        VDJreport = pandas.read_csv(VDJfile, sep="\t")
        UMIreport = pandas.read_csv(UMIfams, sep="\t")            
        ## spikein reads by name
        VDJspikein1 = VDJreport[VDJreport["name"].isin(self.spikes)]
        VDJspikein1["spikein"] = VDJspikein1["name"].map(lambda x : self.spikes[x]) 
        ## native reads by name
        VDJnative = VDJreport[~VDJreport["name"].isin(self.spikes)]
        ## spike in reads by UMI group
        VDJspikein2 = VDJnative[VDJnative["UMI_group"].isin(VDJspikein1["UMI_group"])]
        VDJspikein = pandas.concat([VDJspikein1, VDJspikein2])
        NativeReport = UMIreport[~UMIreport["UMI_group"].isin(VDJspikein["UMI_group"])]
        SpikeinReport = UMIreport[UMIreport["UMI_group"].isin(VDJspikein["UMI_group"])]        
        SpikeinReport = self.cummulative_stats(SpikeinReport)
        NativeReport = self.cummulative_stats(NativeReport)
        ## count spikeins
        SpikeinReport = self.label_spikeins(VDJspikein, SpikeinReport)
        ##label UMI report with consensus spikein label
        ##output files        
        nativeTabOut = "{}_on_target_umi_groups_native.tsv".format(self.basename)
        spikeTabOut = "{}_on_target_umi_groups_spikein.tsv".format(self.basename)
        nativeUMIOut = "{}_cdr3_dedup_report_native.tsv".format(self.basename)
        spikeUMIOut = "{}_cdr3_dedup_report_spikein.tsv".format(self.basename)
        ##write output
        VDJnative.to_csv(nativeTabOut, sep="\t", index=False)
        VDJspikein.to_csv(spikeTabOut, sep="\t", index=False)
        NativeReport.to_csv(nativeUMIOut, sep="\t", index=False)
        SpikeinReport.to_csv(spikeUMIOut, sep="\t", index=False)
        
        
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
            print("Error: ", e)        
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
        description="gather amplicon stats from bam file")
    parser.add_argument("-s", "--spike_aln",
                        required=True,
                        help="Spike in alignment bam")
    parser.add_argument("-v", "--v_aln",
                        required=True,
                        help="V gene alignment bam")
    parser.add_argument("-j", "--j_aln",
                        required=True,
                        help="J gene alignment bam")
    parser.add_argument("-VJ", "--VJ_umi",
                        required=True,
                        help="V & J gene alignment summary file")
    parser.add_argument("-d", "--dedup_report",
                        required=True,
                        help="V & J gene alignment summary file")
    parser.add_argument("-r", "--reference_annot",
                        required=True,
                        help="reference annotation file")
    parser.add_argument("-b", "--basename",
                        required=True,
                        help="basename argument for report file name")
    args = parser.parse_args()
    parser = IPeteSpikeIn(args.spike_aln, args.reference_annot, args.basename)
    parser.split_bams(args.v_aln, args.j_aln)
    parser.split_umi_reports(args.VJ_umi, args.dedup_report)
    ##parser.report_spikein_counts(args.VJaln, "{}_{}".format(
    ##args.basename, "spikein_stats.csv"))


if __name__ == '__main__':
    main()
