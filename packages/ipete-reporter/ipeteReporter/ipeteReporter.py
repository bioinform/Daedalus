import pandas
import math
from SeqNetwork import SeqNetwork
from ipeteMetrics.utils import *

def ipete_reporter(dedupFile, params, basename):
    """
    Reporter for immunoPete pipeline takes the UMI report from ipete-dedup, 
    filters UMI familes for functionality and size, and clusters CDR3.
    
    Parameters
    ----------
    dedupFile : str
    path to ipete-dedup output file
    params : dict 
    command line parameters
    basename : str
    basename for output files

    """
    
    #####################################
    ## get D99, functional UMI families
    #####################################
    dedupDF = pandas.read_csv(dedupFile, sep="\t")
    dedupDF["cummulative_fraction"] = dedupDF["cummulative_reads"] / sum(dedupDF["UMI_family_size"])
    dedupDF = dedupDF.loc[dedupDF["UMI_family_size"] > 1]
    
    ##make sure we have at least a few families before cutting 99 percentile
    if dedupDF.shape[0] > 10:
        dedupDF = dedupDF[dedupDF["cummulative_fraction"] <= 0.99]
    ##fetch functional CDR3
    df = get_functional(dedupDF)
    ##High quality CDR3 sequences
    df = df.loc[df["cdr3_minQual"] >= 20]
    ##reset index values
    df = df.reset_index()

    if df.shape[0] < 10:
        df["cdr3_cluster"] = list(range(0,df.shape[0]))
        df["cdr3_network"] = list(range(0,df.shape[0]))
    else:
        ###################################
        ## cluster CDR3 based on sequences
        ###################################
        graph = SeqNetwork.SeqNetwork()
        for UID, cdr3_AA in zip(df["UMI_group"], df["cdr3_AA"]):
            graph.add_dimension("cdr3_AA", UID, cdr3_AA, params["cdr3_dist"])
        graph.build_graph()

        ######################################
        ## define CDR3 AA similarity networks
        ######################################
        networks = graph.collapse_graph()
        network_labels = graph.retreive_read_labels(networks)
        df["cdr3_network"] = df["UMI_group"].map(network_labels)

        ##remove edges from network based on dedup model, defining clusters
        ##could explore difference parameters for B-cell in future
        graph.remove_edges(params["max_steps"], params["bidding_ratio"])

        ########################
        ## define CDR3 clusters
        ########################
        families = graph.collapse_graph()
        family_labels = graph.retreive_read_labels(families)
        df["cdr3_cluster"] = df["UMI_group"].map(family_labels)

    ########################
    ## determine gene chain
    ########################
    chainVals=[]
    cdr3Ids=[]
    for Vgene, Jgene, cdr3nt in zip(df["Vgene"], df["Jgene"], df["cdr3"]):
        cdr3Id = "{}__{}__{}".format(Vgene, cdr3nt, Jgene)
        if not isinstance(Jgene, str):
            if math.isnan(Jgene):
                ## If the Jgene is missing, the FGxG motif was used to parse the site
                ## use the Vgene chain as the Jgene identity (otherwise these would always be chimera)
                Jgene = Vgene                
        chain = "chimera"            
        if re.search("TRB", Vgene) and re.search("TRB", Jgene):
            chain = "TRB"
        elif re.search("TRA", Vgene) and re.search("TRA", Jgene):
            chain = "TRA"
        elif re.search("TRD", Vgene) and re.search("TRD", Jgene):
            chain = "TRD"
        elif re.search("TRA", Vgene) and re.search("TRD", Jgene):
            ##is Delta V gene shared with alpha (DV5 for example) 
            if re.search("DV", Vgene):
                chain = "TRD"
            ##this gene is similar allele to TRAV38-2DV8 (we see in healthy controls)
            ##can choose to change these manually
            #elif re.search("TRAV38-1", Vgene):
            #    chain = "TRD"
            else:
                chain = "chimera"
        elif re.search("TRG", Vgene) and re.search("TRG", Jgene):
            chain = "TRG"
        elif re.search("IGH", Vgene) and re.search("IGH", Jgene):
            chain = "IGH"
        elif re.search("IGK", Vgene) and re.search("IGK", Jgene):
            chain = "IGK"
        elif re.search("IGL", Vgene) and re.search("IGL", Jgene):
            chain = "IGL"
        else:
            chain = "chimera"
        chainVals.append(chain)
        cdr3Ids.append(cdr3Id)        
    df["chain"] = chainVals
    df["cdr3Id"] = cdr3Ids
       
    ###########################
    ## Collect Diversity stats
    ###########################
    divStats = diversity_stats(df)

    df.to_csv("{}_functional_D99umi_CDR3_report.tsv".format(basename), sep="\t", index=False)
    divStats.to_csv("{}_CDR3_diversity_stats.tsv".format(basename), sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser(
        description="Filter CDR3 dedup report (D99 + functional), reporting high quality CDR3, CDR3 clusters, and divserity stats")
    parser.add_argument("-f", "--cdr3_report",
                        required=True,
                        help="cdr3 dedup report")
    parser.add_argument("-c", "--cdr3_dist",
                        required=True,
                        type = int,
                        help="cdr3 clustering, edit distance")
    parser.add_argument("-s", "--max_steps",
                        required=True,
                        type = int,
                        help="cdr3 clustering, scoring steps/distance")
    parser.add_argument("-r", "--bidding_ratio",
                        required=True,
                        type = int,
                        help="cdr3 clustering, bidding ratio")
    parser.add_argument("-b", "--basename",
                        required=True,
                        help="output basename for report files")
    args = parser.parse_args()    
    params = {}
    params["cdr3_dist"] = args.cdr3_dist
    params["max_steps"] = args.max_steps
    params["bidding_ratio"] = args.bidding_ratio
    ipete_reporter(args.cdr3_report, params, args.basename)

if __name__ == '__main__':
    main()
