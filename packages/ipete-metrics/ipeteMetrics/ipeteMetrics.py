import pandas
import argparse
import collections
import math
import numpy as np
import re
from ipeteMetrics.utils import *

"""
Given a cdr3 dedup report from the ipete pipeline, summarize diversity stats 
for all rearrangements and all functional rearrangements 
"""

def get_ipete_metrics(report_file, basename):
    ipeteDF = import_cdr3_report(report_file)
    functDF = get_functional(ipeteDF)
    ipeteStats = diversity_stats(ipeteDF)
    functStats = diversity_stats(functDF)    
    allOut = "{}_cummulative_diversity_stats_all.csv".format(basename)
    functOut = "{}_cummulative_diversity_stats_functional.csv".format(basename)    
    ipeteStats.to_csv(allOut)
    functStats.to_csv(functOut)
            

def main():
    parser = argparse.ArgumentParser(
        description="Gather cdr3 metrics from ipete stats file")
    parser.add_argument("-c", "--cdr3_summary",
                        required=True,
                        help="alignment summary file")
    parser.add_argument("-b", "--basename",
                        required=True,
                        help="output basename for report files")
    args = parser.parse_args()    

    get_ipete_metrics(args.cdr3_summary, args.basename)

if __name__ == '__main__':
    main()
