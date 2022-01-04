import os
import argparse
import subprocess
import pandas
import glob
from multiprocessing import Pool
import shutil
import logging

def compress_file(file_path):
    logger = logging.getLogger(__name__)
    try:
        subprocess.run([
            "gzip",
            "-f",
            file_path
        ])
    except AssertionError as e:
        logger.warning(e)

def fetch_pattern(outdir, pattern, process):
    files = glob.glob(os.path.join(outdir, pattern))
    for f in files:        
        process.append([f])


def cleanup_run(analysisDir, cpus):
    logger = logging.getLogger(__name__)
    logger.info("cleaning up dir: {}".format(analysisDir))
    ###################
    ##fetch sample info
    ###################
    analysisDir = os.path.realpath(analysisDir)
    sampleInfoFile = os.path.join(analysisDir, "sample_info.csv")
    sampleInfo = pandas.read_csv(sampleInfoFile)
    #######################################
    ##remove all nextflow work directories
    #######################################
    workDir = os.path.join(analysisDir, "config")
    workDirs = []
    for f in os.listdir(workDir):
        toRm = os.path.join(workDir, f, "work")
        if os.path.exists(toRm):
            workDirs.append(toRm)
    logger.info("removing {} work dirs".format(len(workDirs)))
    for d in workDirs:
        shutil.rmtree(d)
    ########################
    ##compress output files
    ########################
    for sample_name in sampleInfo["sample_name"]:
        logger.info("Compressing files for sample: {}".format(sample_name))
        ###############
        ## key folders
        ###############
        sampleData = sampleInfo.loc[sampleInfo["sample_name"] == sample_name]            
        sampleData = sampleData.reset_index()
        ## report directory
        sampleDir = os.path.join(analysisDir, "analysis", sample_name)
        ## output folders
        reporterDir = os.path.join(sampleDir, "ipete_reporter")
        flashDir = os.path.join(sampleDir, "flash")
        dedupDir = os.path.join(sampleDir, "dedup")
        seqDir = os.path.join(sampleDir, "seq/seqtk")
        umiDir = os.path.join(sampleDir, "parse_umi")
        extractUmiDir = os.path.join(sampleDir, "extract_umi")
        trimDir = os.path.join(sampleDir, "trim_primers")
        vdjDir = os.path.join(sampleDir, "VDJ_detect")    
        ########################
        ## compress Fastq files
        ########################
        ## the raw seq file is duplicated in the parse_umi directory
        ## the only difference is that the header contains the UMI info (parsed from R1 and R2)
        toProcess = []
        fetch_pattern(seqDir, "*.fastq", toProcess)
        ## parse UMI outputs
        fetch_pattern(umiDir, "*.fastq", toProcess)
        ## flash output 
        fetch_pattern(flashDir, "*.fastq", toProcess)
        ## trim primers output 
        fetch_pattern(trimDir, "*.fastq", toProcess)
        ## vdj detect output
        fetch_pattern(vdjDir, "*.fastq", toProcess)
        #######################
        ## compress tsv files
        #######################
        ##trim primers tsv reports
        fetch_pattern(trimDir, "*.tsv", toProcess)
        ## on target read reports
        fetch_pattern(vdjDir, "*.tsv", toProcess)
        ## per read UMI report
        fetch_pattern(dedupDir, "*_umi_groups.tsv", toProcess)
        ## extract UMI file (could be removed) 
        fetch_pattern(extractUmiDir, "*.tsv", toProcess)
        if len(toProcess) == 0:
            logger.info("already compressed: {}".format(sample_name))            
        else:
            ###################################
            ##compress in parallel, per sample
            ###################################
            with Pool(cpus) as pool:
                pool.starmap(compress_file, toProcess)
            logger.info("finished compressing output files: {}".format(sample_name))
        

def parse_args():
    parser = argparse.ArgumentParser(description="""Cleanup Ipete pipeline runs.""")
    parser.add_argument('-d', '--rundir',
                        required=True,
                        help='Ipete pipeline output folder.')
    parser.add_argument('-p', '--processors', default='5',
                        type=int,
                        help='number of processors to run')
    args = parser.parse_args()
    return args


def main(args):
    logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(name)s - %(message)s', level=logging.INFO)
    logger = logging.getLogger(__name__)
    
    cleanup_run(args.rundir, args.processors)

if __name__ == '__main__':
    args = parse_args()
    main(args)
