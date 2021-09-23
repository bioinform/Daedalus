

class DNA_translate:
    """
    A class for translating DNA strings
    
    Attributes
    ----------
    codon2AA : amino acid translation table

    """

    def __init__(self):
        self.codon2AA =  {
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
            'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
        }

    
    def translate(self, seq):
        """
        translate DNA string to AA string
        
        Parameters
        ----------
        seq : str
        DNA string

        Returns
        -------
        prot : str
        AA translation

        """
        AA = []
        for i in range(len(seq))[::3]:
            codon = seq[i:i+3]
            AA.append(self.translate_codon(codon))
        prot = "".join(AA)        
        return prot

    def translate_codon(self, codon):
        """
        Translate a DNA codon into an AA residue. 
        
        Parameters
        ----------
        codon : str
        A DNA string
        
        Returns
        -------
        residue : str
        The translated amino acid residue
        """
        residue = ""
        if not codon in self.codon2AA:
            if len(codon) < 3:
                residue =  "_"
            else:
                residue =  "?"
        else:
            residue = self.codon2AA[codon]
        return residue
        
    
    def translate_fuzzy(self, seq):
        """
        Translate DNA strings whose length is not divisable by three.
        Since V and J genes start the boundary. 
        Translation Starts from both ends and meets in the middle with a frameshift.
        
        Parameters
        ----------
        seq : str
        DNA string with frameshift in middle

        Returns
        -------
        prot : str
        AA translation

        """        
        ##translate from beginning to end
        forward=[]
        for i in range(len(seq))[::3]:
            codon = seq[i:i+3]
            forward.append(self.translate_codon(codon))
        ##translate from end to beginning
        reverse=[]
        for i in range(len(seq)-3, 0, -3):
            codon = seq[i:i+3]
            reverse.insert(0, self.translate_codon(codon))            
        ##place frameshift near middle
        begin = int(len(forward)/2)
        AA = []
        for i in range(len(reverse)):
            if i < begin:
                AA.append(forward[i])
            if i == begin:
                AA.append("_")
            if i >= begin:
                AA.append(reverse[i])
        prot = "".join(AA)
        return prot

    
