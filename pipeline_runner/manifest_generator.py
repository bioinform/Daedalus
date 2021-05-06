import argparse
from sample_sheet import SampleSheet
import os
import collections
import pandas
import logging
import xml.etree.ElementTree as ET
from datetime import datetime
import collections

from default_params import *

def getRunData(run_folder):
    """Parse Illumina `RunInfo.xml`, set attributes of run information.
    """
    fname = os.path.join(run_folder, 'RunInfo.xml')
    tree = ET.parse(fname)
    root = tree.getroot()
    run = root.find('Run')
    sequencing_date = run.find('Date').text
    if len(sequencing_date) == 6:
        sequencing_date = datetime.strptime(sequencing_date, '%y%m%d').date()
    elif len(sequencing_date) == 8:
        sequencing_date = datetime.strptime(sequencing_date, '%Y%m%d').date()
    else:
        self.logger.warning(
            'Unrecognized sequencing date format: {}. Record raw string instead.'.format(self.sequencing_date)
        )
    runData = collections.OrderedDict({
        "run_id" : run.attrib['Id'],
        "run_number" : int(run.attrib['Number']),
        "flowcell_id" : run.find('Flowcell').text,
        "instrument" : run.find('Instrument').text,
        "sequencing_date" : sequencing_date
    })
    return runData


def create_manifest(sample_sheet_path, manifest_out, exptParams):
    """
    Take a sample sheet and generate a manifest with default config values for immunoPETE
    
    """
    sample_sheet = SampleSheet(sample_sheet_path)
    
    experiment_name = sample_sheet.header.experiment_name
    instrument_type = sample_sheet.header.instrument_type

    sample_data = []
    for sample in sample_sheet.samples:    
        sampleParams = getDefaultParams()   
        sampData = collections.OrderedDict({
            "experiment_name" : experiment_name,
            "sequencing_platform" : instrument_type,
            "sample_id" : sample.sample_id,
            "sample_name" : sample.sample_name,
            "i7_index_id" : sample.i7_index_id,
            "i7_index" : sample.index,
            "i5_index_id" : sample.i5_index_id,
            "i5_index" : sample.index2,
            "sample_sheet" : sample_sheet_path
        })
        sampleParams = setParamValues(sampleParams, sampData)
        sampleParams = setParamValues(sampleParams, exptParams)
        runData = getRunData(sampleParams["run_folder"])
        sampleParams = setParamValues(sampleParams, runData)
        sample_data.append(sampleParams)
    df = pandas.DataFrame(sample_data)
    df.to_csv(manifest_out, index=False)




    
def parse_args():
    parser = argparse.ArgumentParser(description="""Generate ImmunoPETE Manifest""")
    parser.add_argument('--pipeline_run_id',
                        type=str,
                        required=True,
                        help='provide a name for the pipeline run')
    parser.add_argument('--sequencing_run_folder',
                        type=str,
                        required=True,
                        help='the full path to the sequencing run output folder, matching the sample sheet')
    # parser.add_argument('--sequencing_platform',
    #                     type=str,
    #                     default = "NextSeq",
    #                     choices = ['NextSeq', 'HiSeq', 'NovaSeq', 'MiSeq'],
    #                     required=True,
    #                     help='the sequencing platform used: NextSeq, HiSeq, NovaSeq, MiSeq, etc...')    
    parser.add_argument('--trim_primers',
                        type=bool,
                        default=True,
                        help='Trim V and J primers for ImmunoPETE. default = True')
    parser.add_argument('--spikein',
                        type=bool,
                        default=False,
                        help='Identify ImmunoPETE spikeins. default = False')
    parser.add_argument('--umi_mode',
                        type=bool,
                        default=True,
                        help='use UMIs for dedup/consensus. default = True')
    parser.add_argument('--umi1',
                        default=None,
                        help='parse UMI from read 1. default = ""')
    parser.add_argument('--umi2',
                        type=str,
                        default="NNNNNNNNN",
                        help='parse UMI from read 2. default = "NNNNNNNNN"')
    parser.add_argument('--umi_type', type=str,
                        default='R2', choices = ["R2", "R1", "Both"]) 
    parser.add_argument('--subsample',
                        type=int,
                        default=1,
                        help='Default = 1 (all reads). Any value > 1 will subset the read pairs randomly')
    parser.add_argument('--input_type', type=str,
                        default='DNA', choices = ["DNA", "RNA"],
                        help='Sample input Type')
    parser.add_argument('--output',
                        default="SampleManifest.csv",
                        help='file name for sample manifest')    
    parser.add_argument('sample_sheet',
                        help='sample sheet file name.')
    args = parser.parse_args()
    return args


def main(args):
    exptParams = collections.OrderedDict({
        "trim_primers" : args.trim_primers,
        "umi1" : args.umi1,
        "umi2" : args.umi2,
        "umi_mode" : args.umi_mode,
        "umi_type" : args.umi_type,
        "checkSpikeIn" : args.spikein,
        "inputType" : args.input_type,
        "subsample" : args.subsample,
        "project" : args.pipeline_run_id,
        "run_folder" : args.sequencing_run_folder
    })
    
    create_manifest(args.sample_sheet, args.output, exptParams)
                        

if __name__ == '__main__':
    args = parse_args()
    main(args)
