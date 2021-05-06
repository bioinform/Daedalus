import argparse
import pandas
import os
import re
import collections
import math
from ipeteMetrics.utils import *


def match_files(path, pattern):
    """
    list all files matching a pattern

    Parameters
    ----------
    path : str
        directory location of files
    pattern : str
        A string pattern to search against files, using python re module.

    Returns
    -------
    list of file paths matching pattern

    """
    files = os.listdir(path)
    return ["{}/{}".format(path,x) for x in files if re.search(pattern, x)]

def summarizeChain(df, chain):
    chainTab = df.loc[df["chain"] == chain]
    cloneCounts = []
    for cdr3nt, df in chainTab.groupby(["cdr3Id"]):
        cloneCounts.append(df.shape[0])
    ##return deversity metrics if we have at least 50 cdr3 
    if(sum(cloneCounts) > 10):
        D50 = calc_D50(cloneCounts) ##/ sum(cloneCounts)
        stats = collections.OrderedDict({        
            "{}_renyi_1_norm".format(chain) : renyi_entropy(cloneCounts, alpha=1, normalize=True),
            "{}_renyi_2_norm".format(chain) : renyi_entropy(cloneCounts, alpha=2, normalize=True),
            "{}_simpsons_diversity_norm".format(chain) : simpsons_diversity_index(cloneCounts),
            "{}_simpsons_dominance_norm".format(chain) : simpsons_dominance_index(cloneCounts),        
            "{}_unique_cdr3".format(chain) : len(cloneCounts),
            "{}_D50".format(chain) : D50,
            "{}_gini_index".format(chain) : gini(cloneCounts),
            "{}_cdr3_clusters".format(chain) : len(set(chainTab["cdr3_cluster"]))
        })
    else:
        stats = collections.OrderedDict({        
            "{}_renyi_1_norm".format(chain) : None,
            "{}_renyi_2_norm".format(chain) : None,
            "{}_simpsons_diversity_norm".format(chain) : None,
            "{}_simpsons_dominance_norm".format(chain) : None,        
            "{}_unique_cdr3".format(chain) : len(cloneCounts),
            "{}_D50".format(chain) : None,
            "{}_gini_index".format(chain) : None,
            "{}_cdr3_clusters".format(chain) : None
        })
    return stats
           
 
def pipeline_summary(analysisDir, outFile):    
    analysisDir = os.path.realpath(analysisDir)
    sampleInfoFile = os.path.join(analysisDir, "sample_info.csv")
    #outDir = os.path.join(analysisDir, "analysis", "sample_summary")
    #if not os.path.exists(outDir):
    #    os.mkdir(outDir)
    sampleInfo = pandas.read_csv(sampleInfoFile)
    sample_stats = []
    ##FOR ROW IN TABLE...
    for sample_name in sampleInfo["sample_name"]:
        print("summarizing sample: {}".format(sample_name))
        sampleData = sampleInfo.loc[sampleInfo["sample_name"] == sample_name]            
        sampleData = sampleData.reset_index()
        ## report directories
        sampleDir = os.path.join(analysisDir, "analysis", sample_name)
        reporterDir = os.path.join(sampleDir, "ipete_reporter")
        #metricsDir = os.path.join(sampleDir, "ipete_metrics")
        dedupDir = os.path.join(sampleDir, "dedup")
        seqDir = os.path.join(sampleDir, "seq/seqtk")
        vdjDir = os.path.join(sampleDir, "VDJ_detect")    
        ## count chain types
        chains = {"TRA":0,
                  "TRB":0,
                  "TRG":0,
                  "TRD":0,
                  "IGH":0,
                  "IGK":0,
                  "IGL":0,
                  "chimera":0}
        chainVals = []
        if not os.path.exists(reporterDir):
            ##handle failed runs
            stats = collections.OrderedDict({
                "experiment_name" : sampleData["experiment_name"][0],
                "sample_name" : sample_name,            
                "status" : "fail: incomplete seg2",
                "total_reads" : None,
                "on_target_reads" : None,
                "percent_on_target" : None,
                "UMI_families" : None,
                "UMI_singletons" : None,
                "UMI_replicates" : None,
                "UMI_D99" : None,
                "UMI_D99_funct" : None,
                "UMI_D99_non_funct" : None,        
                "percent_functional" : None,  
                "cell_count" : None,
                "TRA" : None,
                "TRB" : None,
                "TRG" : None,
                "TRD" : None,
                "IGH" : None,
                "IGK" : None,
                "IGL" : None,
                "chimera" : None,
                "IGH_percent" : None,
                "TRB_percent" : None,
                "TRD_percent" : None,
                "unique_CDR3" : None,
                "cdr3_clusters" : None,
                "D50" : None,
                "gini_index" : None,
                "renyi_1_norm" : None,
                "renyi_2_norm" : None,
                "simpsons_diversity_index" : None,
                "simpsons_dominance_index" : None,
            })
            stats.update(trbStats)
            stats.update(ighStats)
            stats.update(trdStats)
            reportFiles = collections.OrderedDict({
                "cdr3_report" : None,
                "diversity_report" : None,
                "dedup_report" : None,
                "on_target_read_report" : None
            })
            stats.update(reportFiles)
            sample_stats.append(stats)
            continue
        ## count reads
        Read1 = match_files(seqDir, "R1_001_seqtk.fastq")[0]
        fqCount = int(sum(1 for line in open(Read1))/4)
        ## diversity stats
        divReport = match_files(reporterDir, "CDR3_diversity_stats.tsv")[0]
        #diversityStats = pandas.read_csv(divReport, sep="\t")
        ## cdr3 report
        cdr3Report = match_files(reporterDir, "functional_D99umi_CDR3_report.tsv")[0]
        cdr3Stats = pandas.read_csv(cdr3Report, sep="\t")
        ## dedup report
        dedupReport = match_files(dedupDir, "cdr3_dedup_report.tsv")[0]
        dedup = pandas.read_csv(dedupReport, sep="\t")
        ## read stat summary
        vdjReport = match_files(vdjDir, "alignment_summary.csv")[0]        
        vdjStats = pandas.read_csv(vdjReport)
        d99 = dedup.loc[dedup["cummulative_fraction"] < 0.99]
        d99 = d99.loc[d99["UMI_family_size"] > 1]
        d99 = d99.reset_index()                    
        ## on target reads
        otrReport = match_files(vdjDir, "on_target.tsv")[0]        
        ## non Functional
        nonFunct = (d99.shape[0] - cdr3Stats.shape[0])
        if d99.shape[0] == 0:
            percent_functional = 0
        else:
            percent_functional = (d99.shape[0] - nonFunct) * 100 / d99.shape[0]
        ##stats
        #diversityVals = diversityStats.tail(1)
        #diversityVals = diversityVals.reset_index()
        ## summarize chain
        cloneCounts = {}
        for chain, cdr3Id in zip(cdr3Stats["chain"], cdr3Stats["cdr3Id"]):
            if not cdr3Id in cloneCounts:
                cloneCounts[cdr3Id] = 0
            cloneCounts[cdr3Id] += 1 
            chains[chain]+=1
        ## diversity stats calculated by chain
        ighStats = summarizeChain(cdr3Stats, "IGH")
        trbStats = summarizeChain(cdr3Stats, "TRB")
        trdStats = summarizeChain(cdr3Stats, "TRD")
        cellCount = chains["TRB"] + chains["TRD"] + chains["IGH"]
        ## Pass or fail the sample
        status = "pass"
        ## minimum 50K sequencing reads        
        if fqCount < 50000:
            status = "fail: <50K reads"
        ## check for "N" in the umi sequence (can happen if read2 fails to be sequenced)
        ## if half or more UMI families have ambigious UMIs, fail the run
        umi1Tot = 0
        umi2Tot = 0
        rows = 0
        for umi1, umi2 in zip(cdr3Stats["UMI_R1"], cdr3Stats["UMI_R2"]):
            rows += 1
            if isinstance(umi1, str):
                if re.search("N", umi1):
                    umi1Tot += 1
            if isinstance(umi2, str):
                if re.search("N", umi2):
                    umi2Tot += 1        
        if umi1Tot >= rows * 0.5:
            status = "fail: bad UMI (>50% N's)"
        if umi2Tot >= rows * 0.5:
            status = "fail: bad UMI (>50% N's)"                        
        ##chain stats
        IGH_percent = 0
        TRB_percent = 0
        TRD_percent = 0
        if sum(chains.values()) > 0:
            IGH_percent = chains["IGH"]/sum(chains.values())
            TRB_percent = chains["TRB"]/sum(chains.values())
            TRD_percent = chains["TRD"]/sum(chains.values())        
        ## summarize sample
        CDR3clones = cdr3Stats["cdr3Id"].value_counts().to_dict()
        D50 = 0
        gini_index = None
        renyi1 = None
        renyi2 = None
        simpsonsDiv = None
        simpsonsDom = None
        if sum(CDR3clones.values()) > 0:            
            D50 = calc_D50(CDR3clones.values())
            gini_index =  gini(list(cloneCounts.values()))
            renyi1 = "{:.6f}".format(renyi_entropy(list(CDR3clones.values()), alpha=1, normalize=True))
            renyi2 = "{:.6f}".format(renyi_entropy(list(CDR3clones.values()), alpha=2, normalize=True))
            simpsonsDiv =  "{:.6f}".format(simpsons_diversity_index(CDR3clones.values()))
            simpsonsDom = "{:.6f}".format(simpsons_dominance_index(CDR3clones.values()))
        if cellCount < 10:
            status = "fail: <10 clones"
        stats = collections.OrderedDict({
            "experiment_name" : sampleData["experiment_name"][0],
            "sample_name" : sample_name,
            "status" : status,
            "total_reads" : fqCount,
            "on_target_reads" : vdjStats["on_target"][0],
            "percent_on_target" :"{:.2f}".format(vdjStats["on_target"][0]*100/fqCount),
            "UMI_families" : dedup.shape[0],
            "UMI_singletons" : dedup[dedup["UMI_family_size"] == 1].shape[0],
            "UMI_replicates" : dedup[dedup["UMI_family_size"] > 1].shape[0],
            "UMI_D99" : d99.shape[0],
            "UMI_D99_funct" : cdr3Stats.shape[0],
            "UMI_D99_non_funct" : d99.shape[0] - cdr3Stats.shape[0],        
            "percent_functional" : percent_functional,  
            "cell_count" : cellCount,
            "TRA" : chains["TRA"],
            "TRB" : chains["TRB"],
            "TRG" : chains["TRG"],
            "TRD" : chains["TRD"],
            "IGH" : chains["IGH"],
            "IGK" : chains["IGK"],
            "IGL" : chains["IGL"],
            "chimera" : chains["chimera"],
            "IGH_percent" : IGH_percent, 
            "TRB_percent" : TRB_percent, 
            "TRD_percent" : TRD_percent, 
            "unique_CDR3" : len(set(cdr3Stats["cdr3_AA"])),
            "cdr3_clusters" : len(set(cdr3Stats["cdr3_cluster"])),
            "D50" : D50,
            "gini_index" : gini_index, ##gini(list(cloneCounts.values())),
            "renyi_1_norm" : renyi1, ##"{:.6f}".format(renyi_entropy(list(CDR3clones.values()), alpha=1, normalize=True)),
            "renyi_2_norm" : renyi2, ##"{:.6f}".format(renyi_entropy(list(CDR3clones.values()), alpha=2, normalize=True)),
            "simpsons_diversity_index" : simpsonsDiv, ## "{:.6f}".format(simpsons_diversity_index(CDR3clones.values())),
            "simpsons_dominance_index" : simpsonsDom ##"{:.6f}".format(simpsons_dominance_index(CDR3clones.values()))
        })
        stats.update(trbStats)
        stats.update(ighStats)
        stats.update(trdStats)        
        reportFiles = collections.OrderedDict({
            "cdr3_report" : cdr3Report,
            "diversity_report" : divReport,
            "dedup_report" : dedupReport,
            "on_target_read_report" : otrReport
        })
        stats.update(reportFiles)
        sample_stats.append(stats)
    ##combine samples
    statsDF = pandas.DataFrame(sample_stats)
    statsDF.to_csv(outFile, index=False)        

def make_arg_parser():
    parser = argparse.ArgumentParser(
        description="Summarize Daedalus pipeline run")
    parser.add_argument("-p", "--pipeline_run",
                        required=True,
                        help="Daedalus pipeline run folder")
    parser.add_argument("-o", "--output",
                        required=False,
                        default = "sample_metrics.csv",
                        help="output file path")    
    return parser

    
def main(args):
    pipeline_summary(args.pipeline_run, args.output)

if __name__ == '__main__':
    arguments = make_arg_parser().parse_args()
    main(arguments)
