import pandas
import argparse
from FastqStreamer.FastqReader import FastqReader
import os


class UmiParser(FastqReader):
    """Parse UID sequences using patterns for Read1 and/or Read2 sequences.
    
    Attributes
    ----------
    read1_pattern : A pattern for UMI on Read1 (e.g. XXNNNN)
    read2_pattern : A pattern for UMI on Read2 (e.g. NNNNNNNNNNNN)
    
    """

    def __init__(self, read1_pattern, read2_pattern):
        """
        Parameters
        ----------
        uid_length : int
        Set the uid size, first X bases of read2 define the UMI
        """
        FastqReader.__init__(self)
        self.r1_umi = self.get_umi_pattern(read1_pattern)
        self.r2_umi = self.get_umi_pattern(read2_pattern)
        
    def get_umi_pattern(self, umi):
        """
        Parameters
        ----------
        umi : str
        the umi pattern to parse
        
        Returns
        -------
        (min,max) tuple of UMI positions
        """
        if not umi == "":
            pattern = []
            i = 0
            for char in umi:                
                if char == "N":
                    pattern.append(i)
                i +=1
            return (min(pattern), max(pattern))
        else:
            return (0,0)
                    

    def read_umi_label(self, read, umi_pattern, read_name, UMI_id):
        """
        Parse UMI from read sequence and quality, labeling the read header with both
        
        Parameters
        ----------
        read : dict
        A dictionary of read elements returned by FastqReader
        umi_pattern : tuple
        coordinates of the UMI in the read
        read_name : str
        The read name
        UMI_id : str
        A name for the UMI being parsed

        Returns
        -------
        newName : str
        The read header with UMI and quality information appended
        """    
        ##cut the UMI from the read seq and quality string.
        ## the second coordinate +1 for capturing the last position
        if umi_pattern == (0,0):
            return read_name             
        UID = read["seq"][umi_pattern[0]:(umi_pattern[1]+1)]
        UID_qual = read["qual"][umi_pattern[0]:(umi_pattern[1]+1)]
        qualVals = []
        for char in UID_qual:
            qualVals.append(ord(char)-33)
        qualString = ",".join([str(x) for x in qualVals])
        ##reformat read name
        newName = "{name}:{UMI_id}__{seq}__{qual}".format(name=read_name,
                                                        UMI_id = UMI_id,
                                                        seq = UID,
                                                        qual = qualString)
        return newName
        
    def label_umi(self, fq1, fq2):
        """Add UID to read1 and read2 headers. write a fastq file with "_uid.fastq" extension.

        Parameters
        ----------
        fq1 : path
        path to read1 fastq file
        fq2 : path
        path to read2 fastq file

        Returns
        -------
        None
        """
        read_pairs = []
        out1 = self.create_fastq_handle(os.path.basename(fq1).replace(".fastq", "_uid.fastq"))
        out2 = self.create_fastq_handle(os.path.basename(fq2).replace(".fastq", "_uid.fastq"))
        for read_pair in self.iter_pairs(fq1, fq2):
            # Ipete UID starts read2
            read_name = read_pair["read1"]["read_name"]
            read_name = self.read_umi_label(read_pair["read1"], self.r1_umi, read_name, "UMI_R1")
            read_name = self.read_umi_label(read_pair["read2"], self.r2_umi, read_name, "UMI_R2")
            read_pair["read1"]["read_name"] = read_name
            read_pair["read2"]["read_name"] = read_name
            self.write_fastq(read_pair["read1"], out1)
            self.write_fastq(read_pair["read2"], out2)
        self.close_file_handles()

def main():
    parser = argparse.ArgumentParser(
        description="Gather UMI information from read1 and/or read2")
    parser.add_argument("-1", "--read1",
                        required=True,
                        help="read1 fastq file")
    parser.add_argument("-2", "--read2",
                        required=True,
                        help="read2 fastq file")
    parser.add_argument("-u1", "--umi_read1",
                        type=str,
                        required=False,
                        default = "",
                        help="UMI Pattern for Read 1. for example XXNNNNNN (X = nonUMI base, N=UMI bases)")
    parser.add_argument("-u2", "--umi_read2",
                        type=str,
                        required=False,
                        default = "",
                        help="UMI Pattern for Read 2. for example NNNNNNNNNNNN (X = nonUMI base, N=UMI bases)")
    args = parser.parse_args()
    if args.umi_read1 == "" and args.umi_read2 == "":
        parser.error("A UMI pattern for read1 and/or read2 is required")
    instance = UmiParser(args.umi_read1, args.umi_read2)
    instance.label_umi(args.read1, args.read2)


if __name__ == '__main__':
    main()
