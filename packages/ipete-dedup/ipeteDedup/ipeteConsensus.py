import collections
import random
from numpy import median
from collections import Counter
import math
from ipeteDedup.DNA_translate import DNA_translate
##from DNA_translate import DNA_translate
import scipy.stats
import pandas

class ipeteConsensus(DNA_translate):
    """
    Define consensus for cdr3 and UMI sequences, using input from the deduping module.
        
    Attributes
    ----------
    passes : int
    count of cdr3 that pass all filters
    lowQual : quality value threshold used to define low Quality.
    """    
    def __init__(self):
        DNA_translate.__init__(self)
        
    def most_abundant_length(self, lens):
        """
        If CDR3 length differ (indels). Find the most abundant length.
        If a tie exists, choose a sequence which is divisable by three (if possible).
        
        Parameters
        ----------
        lens : list
        list of cdr3 nt lengths

        Returns
        -------
        cdr3Len : int
        The most abundant cdr3 nt length found
        """
        ##count frequency of lengths
        lenCounts = Counter(lens)
        ##sort by frequency, then modulo three
        counts = []
        for x in lenCounts:
            counts.append({"len":x, "count":lenCounts[x]})
        sortCounts = sorted(counts, key = lambda x: [-x["count"], x["len"] % 3] )
        cdr3Len = sortCounts[0]["len"]
        return cdr3Len

    def translate_cdr3(self, cdr3_nt):        
        consensusAA = "_"
        if len(cdr3_nt) % 3 == 0:
            consensusAA = self.translate(cdr3_nt)
        else:
            consensusAA = self.translate_fuzzy(cdr3_nt)
        if len(cdr3_nt) < 3:
            consensusAA = "_"
        return consensusAA
    
    def df_consensus(self, umi_fam_df, seqCol, qualCol):
        """
        Given a df with UMI family labels, define consensus sequences and qualities
        for both UMI and CDR3 sequences found in the table
        Parameters
        ----------
        umi_fam_df :pandas DF
        dataframe with UMI group labels

        Returns
        -------
        umi_fam_df :pandas DF
        dataframe with consensus information columns added        
        """
        seqs = {}
        results = []
        UMI_groups = umi_fam_df.UMI_group.unique()        
        for group in UMI_groups:
            seqs[group] = []            
            groupDF = umi_fam_df.loc[umi_fam_df["UMI_group"] == group]            
            if not seqCol in groupDF.columns:
                raise ValueError("column name {} not found".format(seqCol))
            if not qualCol in groupDF.columns:
                raise ValueError("column name {} not found".format(qualCol))            
            for seq, qual in zip(groupDF[seqCol],
                                 groupDF[qualCol]):
                seqs[group].append((seq,qual))
            consensus = self.define_consensus(seqs[group])
            consensus["UMI_group"] = group
            results.append(consensus)
        conDF = pandas.DataFrame(results)
        return conDF        
            
    def define_consensus(self, seqTups, debug=False):
        """
        Given a list of tuples with the format '(string, quality)' define consensus of the string... guided by quality scores
        
        Parameters
        ----------
        seqTups : list of tuples
        A list of tuples, where each tuple contains a string and quality scores

        Returns
        -------
        tuple with consensus information: (consensusSeq, consensusQual, family_size)
        
        consensus : dict 
        dictionary containing consensus information. Sequence, quality, and stats
        """
        baseQuals = collections.OrderedDict()
        baseCounts = collections.OrderedDict()
        family_size = len(seqTups)
        #####################
        # first filter reads
        #####################
        filtTups = self.filter_seqs(seqTups)
        compress = self.compress_seqs(filtTups)

        ####################################################################################
        # if one compressed seq represents a large majority of the reads, we have consensus
        ####################################################################################
        ## >= 9 reads, it is also reasonable to assume Q60        
        compressCounts = sorted(compress.items(), key=lambda x: -len(x[1]))
        if len(compressCounts[0][1]) >= (family_size * 0.9) and family_size >= 10:
            consensusSeq = compressCounts[0][0]
            consensusQual = "]"*len(consensusSeq)
            consensus = {"consensus_seq":consensusSeq,
                         "consensus_qual":consensusQual}
            return consensus
        ##################################
        ##if singleton, we have consensus
        ##################################
        if family_size == 1:
            consensus = {"consensus_seq":filtTups[0][0],
                         "consensus_qual":filtTups[0][1]}            
            return consensus
        ###################################################
        # otherwise, find consensus by weighing base-pairs
        ###################################################
        for seq in compress:
            seqCount = len(compress[seq])
            # count the frequency of each base
            quals = compress[seq]
            ##parse quality values from string
            pVals = {}
            for qual in quals:
                for pos, q, in enumerate(qual):
                    if not pos in pVals:
                        pVals[pos] = []
                    qVal = ord(q) - 33
                    pVal = 10**(-qVal/10)
                    pVals[pos].append(pVal)                    
            ##sum quality values
            for pos, base in enumerate(seq):
                if not pos in baseCounts:
                    baseQuals[pos] = {"A": [], "C": [], "T": [], "G": [], "N": []}
                    baseCounts[pos] = {"A": 0, "C": 0, "T": 0, "G": 0, "N": 0}    
                baseQuals[pos][base] += pVals[pos]
                baseCounts[pos][base] += len(quals)
        ##combine pvalues for each base
        basePVals = {}
        for pos in baseQuals:
            basePVals[pos] = {"A": 0, "C": 0, "T": 0, "G": 0, "N": 0}
            for base in baseQuals[pos]:
                if len(baseQuals[pos][base]) > 0:
                    pValue = scipy.stats.combine_pvalues(baseQuals[pos][base])[1]
                    if pValue > 0.63:
                        pValue = 0.63
                    basePVals[pos][base] = pValue                    
                else:
                    basePVals[pos][base] = 1

        ##################################
        # determine best base per position
        ##################################
        forConsensus = []
        forQual = []
        for pos in baseCounts:
            baseVals = sorted(basePVals[pos].items(), key=lambda x: x[1])            
            ##minVal == lowest P-value
            minVal = min([x[1] for x in baseVals])
            choices = []
            for base, val in baseVals:
                if val == minVal:
                    choices.append(base)
            maxBase = ""            
            ################################################
            # if more than one maximum, randomly choose one
            ################################################
            if len(choices) > 1:                
                ## choose highest coverage otherwise
                toChoose = []
                #####################################
                ## filter by coverage if equal pvalue
                #####################################
                maxCov = 0
                for choice in choices:
                    cov = baseCounts[pos][choice]
                    if cov > maxCov:
                        maxCov = cov
                    toChoose.append({"base": choice, "cov": cov})
                toChoose = [x for x in toChoose if x["cov"] == maxCov]
                if len(toChoose) > 1:
                    maxBase = choices[random.randint(0, len(choices)-1)]                    
                else:                    
                    maxBase = toChoose[0]["base"]
                    choices = [maxBase]
            else:
                maxBase = choices[0]
            forConsensus.append(maxBase)
            ###############################
            ## calculate consensus quality
            ###############################
            ##if more than one choice existed, the quality is set to 2 for that position
            if len(choices) > 1:
                forQual.append(2)
            else:
                basePVal = basePVals[pos][maxBase]
                if basePVal <= 0.000001:
                    baseQVal = 60
                else:
                    baseQVal = int(-10*math.log10(basePVal))
                    if baseQVal == 0:
                        raise ValueError("error {}".format(len(basePVal)))                        
                forQual.append(baseQVal)            
        #############################
        ##calculate consensus quality
        #############################
        consensusSeq = "".join(forConsensus)
        consensusQual=""
        consensusQual = "".join([chr(x+33) for x in forQual])
        consensus = {"consensus_seq":consensusSeq,
                     "consensus_qual":consensusQual}
        return consensus

    def compress_seqs(self, seqTups):
        """
        Find Exact matching sequence strings for a list of tuples.
        store qualities within exact matching strings
        
        parameters
        ----------
        seqTups : list
        a list of (seq, qual) tuples

        Returns
        -------
        compress : dict
        A dictionary of seqs with quals as values
        """
        compress = {}
        for seq, qual in seqTups:
            if not seq in compress:
                compress[seq] = []
            compress[seq].append(qual)
        return compress
        
    def filter_seqs(self, seqTups):
        """
        Filter sequences, keeping the most abundant length        

        parameters
        ----------
        seqTups : list
        a list of (seq, qual) tuples
        
        Returns
        -------
        seqs : list
        a list of filtered sequence tuples
        """
        lens = [len(x[0]) for x in seqTups]
        filtLen = self.most_abundant_length(lens)
        ## only keep sequences with same length
        seqs = [x for x in seqTups if len(x[0]) == filtLen]            
        return seqs
        
        
    
if __name__ == '__main__':
    main()
