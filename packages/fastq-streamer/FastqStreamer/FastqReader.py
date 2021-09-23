class FastqReader:
    """
    Utility class for handling fastq formated files

    Attributes
    ----------
    file_handles_ : A dictionary for file handles, used for writing one or more
    Fastq file.

    """

    def __init__(self):
        self.file_handles_ = {}

    def iter_pairs(self, fastq1, fastq2):
        """
        An iterator, returning read pair info from two fastq files

        Parameters
        ----------
        fastq1 : str
        Path to Fastq file for read1
        fastq2 : str
        Path to Fastq file for read2

        Returns
        -------
        record : dictionary  conatining read pair info (name, seq, qual)
        """
        with open(fastq1) as fq1, open(fastq2) as fq2:
            count = 0
            record = {"read1": {}, "read2": {}}
            for read1, read2, in zip(fq1, fq2):
                if count % 4 == 0:
                    record = {"read1": {}, "read2": {}}
                record = self.read_pair_record(read1, read2, count, record)
                if count % 4 == 3:
                    yield record
                count += 1

    def iter_reads(self, fastq):
        """
        An iterator, returning read info from fastq

        Parameters
        ----------
        fastq : str
        Path to Fastq file

        Returns
        -------
        record : dictionary  containing read info
        """
        with open(fastq) as fq:
            count = 0
            record = {"read1": {}, "read2": {}}
            for read in fq:
                if count % 4 == 0:
                    record = {"read1": {}, "read2": {}}
                record = self.read_struct(read, count, record, "read1")
                if count % 4 == 3:
                    yield record["read1"]
                count += 1

    def read_pair_record(self, read1, read2, count, record):
        """
        combines read pair information

        Parameters
        ----------
        read1 : str
        line from read1 fastq file
        read2 : str
        line from read2 fastq file

        Returns
        -------
        record : dictionary  conatining read_pair info
        """
        record = self.read_struct(read1, count, record, "read1")
        record = self.read_struct(read2, count, record, "read2")
        return record

    def read_struct(self, line, count, record, which_read):
        """
        Store fastq file information by line

        Parameters
        ----------
        line : str
        line from fastq file
        count : int
        position in fastq file
        record : dict
        A dictionary for storing read information
        which_read : str
        Dictionary key, to keep track of read1/read2

        Returns
        -------
        record : dictionary  conatining read info
        """

        line = line.strip("\n")
        pos = count % 4
        if pos == 0:
            record[which_read]["read_name"] = line.replace("@", "").split(" ")[0]
        if pos == 1:
            record[which_read]["seq"] = line
        if pos == 3:
            record[which_read]["qual"] = line
        return record

    def create_fastq_handle(self, filename):
        """
        Maintain one or more fastq file handles

        Parameters
        ----------
        filename : str
        name of fastq file to be written

        Returns
        -------
        filename, a key for the file handle created

        """
        if not filename in self.file_handles_:
            self.file_handles_[filename] = open(filename, 'w')
        return filename

    def write_fastq(self, read, filename):
        """
        Write read record to fastq file

        Parameters
        ----------
        read : dict
        dictionary record, containing read information
        filename : str
        name of fastq file to be written        

        Returns
        -------
        None
        """
        fh = self.file_handles_[filename]
        fastq_lines = "{name}\n{seq}\n{third}\n{qual}\n"
        trimmed_read = fastq_lines.format(name="@"+read["read_name"],
                                          seq=read["seq"],
                                          third="+",
                                          qual=read["qual"])
        fh.write(trimmed_read)

    def close_file_handles(self):
        """
        close all fastq files written

        Returns
        -------
        None
        """
        for filename in self.file_handles_:
            self.file_handles_[filename].close()
