import argparse
import pandas
import logging
from ipeteDedup.getLogger import getLogger
##from getLogger import getLogger
import collections
from SeqNetwork import SeqNetwork
from ipeteDedup.ipeteConsensus import ipeteConsensus
##from ipeteConsensus import ipeteConsensus
import Levenshtein
import os

def filter_reads(aln_df, params):
    """
    Filter a pandas dataframe (df) against a set of given parameters (params)

    Parameters
    ----------
    aln_df : pandas dataframe
    A pandas dataframe or dataframe chunk 
    params : dict
    A dictionary of parameters to be filtered against

    Returns
    -------
    filt : pandas dataframe
    A pandas dataframe filtered against input parameters
    """
    filt = aln_df[aln_df["cdr3"].str.len() > 0]
    filt = filt[filt["cdr3"].str.len() < params["max_cdr3_length"]]
    filt = filt[filt["avg_read_qual"] > params["min_qual"]]        
    return filt

def check_umi_columns(aln_df, params):
    """
    UMI1 and UMI2 columns are optional, if the values are not found, insert NULL values

    Parameters
    ----------
    aln_df : pandas dataframe
    A pandas dataframe or dataframe chunk 
    params : dict
    A dictionary of parameters to be filtered against

    Returns
    -------
    aln_df : pandas dataframe
    A pandas dataframe with columns values for both `umi_r1` and `umi_r2` 
    """
    if not "umi_r1" in aln_df.columns:
        if params["umi_read1"] is True:
            raise ValueError("Read1 UMI not found, please provide UMI1 args to parse-umi")
        aln_df["umi_r1"] = None
        aln_df["umi_r1_qual"] = None
    if not "umi_r2" in aln_df.columns:
        if params["umi_read2"] is True:
            raise ValueError("Read2 UMI not found, please provide UMI2 args to parse-umi")
        aln_df["umi_r2"] = None
        aln_df["umi_r2_qual"] = None    
    return aln_df

def add_dimensions(aln_df, graph, params):
    """
    Given a pandas dataframe with UMI1 and UMI2 columns, build a graph with cdr3 and UMI(s) as dimensions

    Parameters
    ----------
    aln_df : pandas dataframe
    A pandas dataframe or dataframe chunk 
    graph : SeqNetwork instance
    An instance of SeqNetwork, used for building string similarity networks
    params : dict
    A dictionary of parameters to be filtered against

    Returns
    -------
    graph : SeqNetwork instance
    A graph with dimensions added
    """

    for rid, umi1, umi1_qual, umi2, umi2_qual, \
        cdr3, cdr3_qual, read_qual in zip(aln_df['name'],
                                          aln_df['umi_r1'],
                                          aln_df['umi_r1_qual'],
                                          aln_df['umi_r2'],
                                          aln_df['umi_r2_qual'],
                                          aln_df['cdr3'],
                                          aln_df['cdr3_qual'],
                                          aln_df['avg_read_qual']):            
        if params["umi_read1"] is True:
            graph.add_dimension("UMI_R1", rid, umi1, params["umi_dist"])
        if params["umi_read2"] is True:
            graph.add_dimension("UMI_R2", rid, umi2, params["umi_dist"])
        graph.add_dimension("CDR3", rid, cdr3, params["cdr3_dist"])
    return graph
    
def ipete_dedup(alnFile, params, basename, logger):
    """
    Takes alignment file with parameters and produces dedup report files and logs progress

    Parameters
    ----------
    alnFile : tsv file
    alignment report returned by VDJdetector and with UMIs from extract-umi    
    params : dict
    A dictionary of parameters to be filtered against
    basename : string
    Basename to use for report files
    logger : logger
    Python logger for project

    Returns
    -------
    None

    """
    dedup_output = "{}_on_target_umi_groups.tsv".format(basename)
    dedup_report_file = "{}_cdr3_dedup_report.tsv".format(basename)
    ## read file in chunks
    chunksize = 100000    
    ## instatiate graph
    graph = SeqNetwork.SeqNetwork()
    consensus = ipeteConsensus()
    ## add dimensions to graph
    logger.info("add dimensions...")
    for chunk in pandas.read_csv(alnFile, chunksize=chunksize, sep="\t"):        
        filt = filter_reads(chunk, params)        
        filt = check_umi_columns(filt, params)                    
        graph = add_dimensions(filt, graph, params)
    ##build the graph
    logger.info("build graph...")
    graph.build_graph()
    ##cluster
    logger.info("partition graph...")
    graph.remove_edges(params["max_steps"], params["bidding_ratio"])
    ## define UMI family labels
    logger.info("define UMI families...")
    families = graph.collapse_graph()
    logger.info("label UMI families...")
    family_labels = graph.retreive_read_labels(families)
    ####################################
    ## read entire file, label families
    ####################################
    aln_df = pandas.read_csv(alnFile, sep="\t")
    filtered = aln_df[aln_df["name"].isin(list(family_labels.keys()))].copy(deep=True)
    filtered["UMI_group"] = filtered["name"].map(family_labels)
    filtered.to_csv(dedup_output, sep="\t", index=False)
    ##filtered = consensus.define_UMI_CDR3_consensus(filtered)
    ################################
    ## consensus and summary report
    ################################
    umi_summary = umi_group_report(filtered, params)    
    cdr3Consensus = consensus.df_consensus(filtered, "cdr3", "cdr3_qual")        
    ##cdr3Consensus.rename(columns={"consensus_seq":"cdr3_AA", "consensus_qual":"cdr3_qual"}) 
    umi_summary["cdr3"] = cdr3Consensus["consensus_seq"]
    umi_summary["cdr3_AA"] = umi_summary["cdr3"].apply(consensus.translate_cdr3)
    umi_summary["cdr3_qual"] = cdr3Consensus["consensus_qual"]
    umi_summary["cdr3_minQual"] = cdr3Consensus["consensus_qual"].apply(qualMin)
    ##rename columns    
    if params["umi_read1"] is True:
        umi1Consensus = consensus.df_consensus(filtered, "umi_r1", "umi_r1_qual")
        umi_summary["UMI_R1"] = umi1Consensus["consensus_seq"]
        umi_summary["UMI_R1_qual"] = umi1Consensus["consensus_qual"]
        umi_summary["UMI_R1_minQual"] = cdr3Consensus["consensus_qual"].apply(qualMin)
    if params["umi_read2"] is True:
        umi2Consensus = consensus.df_consensus(filtered, "umi_r2", "umi_r2_qual")
        umi_summary["UMI_R2"] = umi2Consensus["consensus_seq"]
        umi_summary["UMI_R2_qual"] = umi2Consensus["consensus_qual"]
        umi_summary["UMI_R2_minQual"] = cdr3Consensus["consensus_qual"].apply(qualMin)
    ## order by UMI family size and then min quality
    umi_summary = umi_summary.sort_values(by=['UMI_family_size', 'cdr3_minQual'], ascending=False)    
    umi_summary["cummulative_reads"] = umi_summary["UMI_family_size"].cumsum()
    umi_summary["cummulative_fraction"] = umi_summary["cummulative_reads"] / sum(umi_summary["UMI_family_size"])
    umi_summary.to_csv(dedup_report_file, sep="\t", index=False)
    summarize(aln_df, umi_summary, graph, logger)
    
def qualMin(qualString):
    minQual = 60
    qualChars = set(list(qualString))
    for qual in qualString:
        qVal = ord(qual) - 33
        if qVal < minQual:
            minQual = qVal
    return minQual

def umi_group_report(umi_fam_df, params):
    """
    given a VDJ report table with UMI_group labels, generate V-J gene stats by UMI family
    """
    family_report = []
    UMI_groups = umi_fam_df.UMI_group.unique()        
    for group in UMI_groups:
        groupDF = umi_fam_df.loc[umi_fam_df["UMI_group"] == group]            
        family_size = groupDF.shape[0]
        Vgenes = groupDF["v_gene"].value_counts(dropna=False)
        VGeneType = groupDF["v_locus_type"].value_counts(dropna=False)            
        VIdbType = groupDF["v_idb_type"].value_counts(dropna=False)            
        Jgenes = groupDF["j_gene"].value_counts(dropna=False)
        JGeneType = groupDF["j_locus_type"].value_counts(dropna=False)            
        JIdbType = groupDF["j_idb_type"].value_counts(dropna=False)                
        family_report.append(
            collections.OrderedDict({"cdr3_AA" : "",
                                     "UMI_family_size" : family_size,
                                     "UMI_R1" : "",
                                     "UMI_R1_qual" : "",
                                     "UMI_R1_minQual" : None,
                                     "UMI_R2" : "",
                                     "UMI_R2_qual" : "",                                     
                                     "UMI_R2_minQual" : None,
                                     "cdr3" : "",
                                     "cdr3_qual" : "",
                                     "cdr3_minQual" : None,
                                     "Vgene" : Vgenes.idxmax(),
                                     "Jgene" : Jgenes.idxmax(),
                                     "multiV" : len(Vgenes),
                                     "multiJ" : len(Jgenes),
                                     "V_GeneType" : VGeneType.idxmax(),
                                     "V_IdbType" : VIdbType.idxmax(),
                                     "J_GeneType" : JGeneType.idxmax(),
                                     "J_IdbType" : JIdbType.idxmax(),
                                     "UMI_group" : group
            })
        )
    return pandas.DataFrame(family_report)
    

def summarize(aln_df, umi_summary, graph, logger):
    logger.info("on target reads = {}".format( aln_df.shape[0] ))
    logger.info("filtered reads = {}".format( sum(umi_summary["UMI_family_size"]) ))        
    dimensions = list(graph.dimensions.keys())
    logger.info("number of dimensions = {}".format(len(graph.dimensions)))
    for dimName in dimensions:
        dimTotal = len(set(graph.dimensions[dimName].values()))
        logger.info("{} unique seqs = {}".format(dimName, dimTotal))
    logger.info("number of nodes = {}".format(len(graph.graph)))
    logger.info("total centroids = {}".format(graph.total_centroids))
    logger.info("UMI families = {}".format(umi_summary.shape[0]))                

            
def main():    
    parser = argparse.ArgumentParser(
        description="Deduper designed for identifying UMI families and consensus sequences from immunoPETE sequencing runs.")
    parser.add_argument("-a", "--alignment_summary",
                        required=True,
                        help="alignment summary file")
    parser.add_argument("-u1", "--umi_read1",                        
                        required=False,
                        action='store_true',
                        default = False,
                        help="dedup with read1 UMI")
    parser.add_argument("-u2", "--umi_read2",                        
                        required=False,
                        action='store_true',
                        default = False,
                        help="dedup with read2 UMI")
    parser.add_argument("-u", "--umi_dist",                        
                        required=True,
                        type=int,
                        help="maximum edit distance allowed when searching UMIs")
    parser.add_argument("-c", "--cdr3_dist",
                        required=True,
                        type=int,
                        help="maximum edit distance allowed when searching CDR3s")
    parser.add_argument("-q", "--min_qual",
                        required=True,
                        type=int,
                        help="minimum average read quality to keep a read.")
    parser.add_argument("-l", "--max_cdr3_length",
                        required=True,
                        type=int,
                        default=100,
                        help="maximum CDR3 nt sequence length allowed (otherwise read is filtered)")
    parser.add_argument("-b", "--bidding_ratio",
                        required=True,
                        type=int,
                        default=2,
                        help="The ratio when a smaller node will accept the bid of a larger node.")
    parser.add_argument("-s", "--steps",
                        required=True,
                        type=int,
                        default=2,
                        help="number of edges to tranverse")
    parser.add_argument("-p", "--processors",
                        required=False,
                        type=int,
                        default=1,
                        help="number of processors used for deduping")
    parser.add_argument("-o", "--out_basename",
                        required=True,
                        help="basename for output summary reports: read report, cdr3 consensus report and umi-cdr3 node report")

    args = parser.parse_args()        
    alnFile = args.alignment_summary
    logger = getLogger(
        os.path.basename(alnFile).replace(".tsv", ""), debug=False)        
    params = {
        "umi_dist" : args.umi_dist,
        "cdr3_dist" : args.cdr3_dist,
        "min_qual" : args.min_qual,
        "max_cdr3_length" : args.max_cdr3_length,
        "bidding_ratio" : args.bidding_ratio,
        "max_steps" : args.steps,
        "processors" : args.processors,        
        "umi_read1" : args.umi_read1,
        "umi_read2" : args.umi_read2
    }
    ipete_dedup(alnFile, params, args.out_basename, logger)
    
    
if __name__ == '__main__':
    main()
