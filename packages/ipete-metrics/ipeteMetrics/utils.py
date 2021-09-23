import pandas
import argparse
import collections
import math
import numpy as np
import re

"""
renyi_0 : Hartley/max entropy
renyi_1 : shannon entropy
renyi_2 : Collision entropy
renyi_inf : Min-entropy
"""

def calc_freq(vals):
    freqs = []
    total = float(sum(vals))
    for val in vals:
        prob = float(val/total)
        freqs.append(prob)
    return freqs

def gini(arr):
    if len(arr) < 2:
        return(1)
    ## first sort
    arr = np.array(arr)
    sorted_arr = arr.copy()
    sorted_arr.sort()
    n = arr.size
    coef_ = 2. / n
    const_ = (n + 1.) / n
    weighted_sum = sum([(i+1)*yi for i, yi in enumerate(sorted_arr)])
    return coef_*weighted_sum/(sorted_arr.sum()) - const_

def shannon_entropy(vals):
    """
    Calculate Shannon Entropy for a set of values
    """
    ## vals represent cdr3 counts
    calcs = []
    freqs = calc_freq(vals)
    for freq in freqs:
        calcs.append(freq * math.log2(freq))
    return sum(calcs) * -1

def renyi_calc(vals, alpha):
    """
    Calculate Shannon Entropy for a set of values
    """
    ## vals represent cdr3 counts
    calcs = []
    freqs = calc_freq(vals)
    for freq in freqs:
        calcs.append(math.pow(freq,alpha))    
    return (1/(1-alpha)) * math.log2(sum(calcs))

def simpsons_diversity_index(vals):
    """
    Calculate Simpsons Diveristy for a set of values
    """
    ## vals represent cdr3 counts
    total = float(sum(vals))
    if total <= 1:
        return 0
    multp = []
    for val in vals:
        multp.append(val * (val-1))
    simpsons =  sum(multp) / ( total * (total - 1) )
    return 1 - simpsons

def simpsons_dominance_index(vals):
    """
    Calculate Simpsons Dominance Index for a set of values

    Parameters
    ----------
    vals : list
    list of float values
    """

    total = float(sum(vals))
    return sum([float(x/total)**2 for x in vals])

def renyi_entropy(vals, alpha=1, normalize = False):
    """
    Calculate Renyi Entropy for a set of values

    Parameters
    ----------
    vals : list
    list of float values

    """

    alpha = float(alpha)
    ##normalization factor
    norm = 1
    if normalize == True:
        if len(vals) > 1:
            norm = math.log2(len(vals))
    ## alpha == 1: Shannon Entropy
    if alpha == 1:
        return shannon_entropy(vals) / norm
    ## boundary alpha values
    elif alpha == 'inf' or alpha == np.inf:
        return math.log2(np.max(vals)) / norm
    ## alpha != 1: Renyi Entropy
    elif alpha > 1:
        return renyi_calc(vals, alpha) / norm

def calc_D50(vals):
    """
    Calculate the D50 for clone counts

    Parameters
    ----------
    vals : list
    list of float values

    Returns
    -------
    D50 : float
    fraction of clones containing the top 50% of all UMI families        
    """
    sortVals = sorted(vals, key = lambda x: -x)
    mid = sum(sortVals)/2
    tot = 0
    D50 = 0
    for x in sortVals:
        D50 += 1
        tot += x
        if tot >= mid:
            break
    return float(D50 / len(sortVals))


def gini(arr):
    if len(arr) < 2:
        return(1)
    ## first sort
    arr = np.array(arr)
    sorted_arr = arr.copy()
    sorted_arr.sort()
    n = arr.size
    coef_ = 2. / n
    const_ = (n + 1.) / n
    weighted_sum = sum([(i+1)*yi for i, yi in enumerate(sorted_arr)])
    return coef_*weighted_sum/(sorted_arr.sum()) - const_


def import_cdr3_report(report_file):
    df = pandas.read_csv(report_file, sep="\t")
    return df

def get_functional(df):
    funDF = df.loc[df["V_IdbType"] == "gene"]
    funDF = funDF.loc[funDF["J_IdbType"] == "gene"]
    funDF = funDF.loc[~funDF.cdr3_AA.str.contains('\\*|\\?|\\_', regex=True, na=False)]
    funDF = funDF.reset_index()
    return funDF

def run_calc(sliceDF):    
    sliceDF['cdr3Id'] = sliceDF[['Vgene', 'cdr3', 'Jgene']].apply(lambda x: '_'.join(str(x)), axis=1)
    sliceCDR3 = sliceDF["cdr3Id"].value_counts().to_dict()
    vals = collections.OrderedDict({
        "total_umi_families" : sliceDF.shape[0],        
        "UMI_family_size" : sliceDF["UMI_family_size"].iloc[-1],
        "cummulative_reads" : sliceDF["cummulative_reads"].iloc[-1],
        "cummulative_fraction" : sliceDF["cummulative_fraction"].iloc[-1],
        "unique_cdr3": len(sliceCDR3), 
        "simpsons_diversity_index": simpsons_diversity_index(sliceCDR3.values()),
        "simpsons_dominance_index": simpsons_dominance_index(sliceCDR3.values()),
        "renyi_0": renyi_entropy(sliceCDR3.values(), alpha=0),
        "renyi_1": renyi_entropy(sliceCDR3.values(), alpha=1),
        "renyi_2": renyi_entropy(sliceCDR3.values(), alpha=2),
        "renyi_inf": renyi_entropy(list(sliceCDR3.values()), alpha='inf'),
        "renyi_0_norm": renyi_entropy(sliceCDR3.values(), alpha=0, normalize = True),
        "renyi_1_norm": renyi_entropy(sliceCDR3.values(), alpha=1, normalize = True),
        "renyi_2_norm": renyi_entropy(sliceCDR3.values(), alpha=2, normalize = True),
        "renyi_inf_norm": renyi_entropy(list(sliceCDR3.values()), alpha='inf', normalize = True),
        "D50": calc_D50(sliceCDR3.values())
    })
    return vals
    
    
def diversity_stats(df):
    """
    given a pandas data frame, calculate diversity/entropy stats one UMI family at a time.
    """
    stats = []
    ##no need to do all rows, slow with larger datasets
    if df.shape[0] >= 100000:
        for i in range(5000,df.shape[0],5000):
            sliceDF = df[:i]        
            vals = run_calc(sliceDF)
            stats.append(vals)            
    elif df.shape[0] >= 50000:
        for i in range(1000,df.shape[0],1000):
            sliceDF = df[:i]        
            vals = run_calc(sliceDF)
            stats.append(vals)
    elif df.shape[0] >= 10000:
        for i in range(500,df.shape[0],500):
            sliceDF = df[:i]        
            vals = run_calc(sliceDF)
            stats.append(vals)
    elif df.shape[0] >= 5000:
        for i in range(100,df.shape[0],100):
            sliceDF = df[:i]        
            vals = run_calc(sliceDF)
            stats.append(vals)
    elif df.shape[0] >= 1000:
        for i in range(20,df.shape[0],20):
            sliceDF = df[:i]        
            vals = run_calc(sliceDF)
            stats.append(vals)
    elif df.shape[0] >= 500:
        for i in range(10,df.shape[0],10):
            sliceDF = df[:i]        
            vals = run_calc(sliceDF)
            stats.append(vals)
    elif df.shape[0] >= 100:
        for i in range(2,df.shape[0],2):
            sliceDF = df[:i]        
            vals = run_calc(sliceDF)
            stats.append(vals)
    else:
        for i in range(1,df.shape[0]):
            sliceDF = df[:i]        
            vals = run_calc(sliceDF)
            stats.append(vals)
    summaryDF = pandas.DataFrame(stats)            
    return summaryDF
