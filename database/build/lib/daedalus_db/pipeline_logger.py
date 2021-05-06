#!/usr/bin/env python
import argparse
import glob
import gzip
import logging
import os
import sys
import datetime

import numpy as np
import pandas as pd
import collections
from collections import defaultdict

from daedalus_db.daedalus_db import DaedalusDB
from daedalus_db.run_info import RunInfo
from sample_sheet import SampleSheet

class PipelineLogger:
    """Log Daedalus manifest and pipeline run info

    Attributes
    ----------
    sample_info_fname : str
        Path to the sample_info.csv in the pipeline run folder. 
    sample_sheet_fname : str
        Path to Illumina sample sheet.
    seq_run_path : str
        Path to raw sequencing run.
    pipe_run_path : str
        Path to immunoPETE pipeline run.
    records : pandas.DataFrame
        Dataframe contains instrument, run_number, flowcell_id, sequencing_date, sample_name, lane, i7_index_id,
        index, i5_index_id, index2, cdr3_nt, cdr3_aa, cdr3_qual, v, d, j, umi_seq, umi_qual, family_size.
    """
    def __init__(self, sample_manifest_fname, url='sqlite:///daedalus.db', echo=False, newDb=False):
        """

        Parameters
        ----------
        sample_manifest_fname : str
            Path to the sample_info.csv in the pipeline run folder. 
        url : str
            URL to the database.
        echo : bool
            Whether to echo SQL statement when executing.
        """        
        self.sample_manifest_fname = sample_manifest_fname
        self.sample_sheet_fname = None        
        self.seq_run_path = None
        self.pipe_run_path = None
        self.logger = logging.getLogger(__name__)
        self.daedalus_db = DaedalusDB(url, echo, newDb)
        self.records = None

    def get_param_info(self):
        paramCols = ['subsample', 'subsampleSeed', 'primer_targets', 'inputType',
                     'trim_primers', 'checkSpikein', 'umi_mode', 'umi1', 'umi2', 'umi_type',
                     'vPrimerRef', 'jPrimerRef', 'vGeneRef', 'jGeneRef', 'primerRef', 'geneRef', 'spikeRef',
                     'vGeneScore', 'vGeneKmer', 'vGeneFrac', 'jGeneScore', 'jGeneKmer', 'jGeneFrac',
                     'vPrimScore', 'vPrimKmer', 'vPrimFrac', 'jPrimScore', 'jPrimKmer', 'jPrimFrac',
                     'spikeScore', 'spikeKmer', 'spikeFrac',
                     'min_read_qual', 'minGeneIdentity', 'max_cdr3_length',
                     'bidding_ratio', 'max_steps', 'umi_edit_dist', 'cdr3_edit_dist']        
        param_info = self.sample_manifest[paramCols]
        return param_info
    
    def get_analysis_info(self):
        analysis_info = pd.DataFrame({
            'analysis_name' : self.sample_manifest["project"][0],
            'analysis_folder' : self.pipe_run_path,
            'analysis_date' : datetime.date.today(),
            'pipeline_version' : 'omics bricks',
            'segment_1_status' : 'pending',
            'segment_2_status' : 'pending',
            'segment_3_status' : 'pending'            
        }, index=[0])
        return analysis_info

    def get_sequencing_info(self, seq_run_path):
        """Get sequencing run information from RunInfo.xml.

        Parameters
        ----------
        seq_run_path : str
            Path to the raw Illumina sequencing run.

        Returns
        -------
        run_info : pandas.DataFrame
            Dataframe contains the sequencing run information: instrument, run_number, flowcell_id, sequencing_date.
        """
        run_info = RunInfo(os.path.join(seq_run_path, 'RunInfo.xml'))
        run_info = pd.DataFrame({
            'sequencing_platform' : self.sample_manifest["sequencing_platform"][0],
            'instrument': run_info.instrument,
            'run_number': run_info.run_number,
            'flowcell_id': run_info.flowcell_id,
            'sequencing_date': run_info.sequencing_date,
            'sample_sheet' : self.sample_manifest["sample_sheet"][0],
            'run_folder' : self.sample_manifest["run_folder"][0],
            
        }, index=[0])
        return run_info

            
    def get_sample_info(self, sample_sheet_fname):
        """Get sample information from sample sheet.

        Parameters
        ----------
        sample_sheet_fname : str
            Path to Illumina sample sheet.

        Returns
        -------
        sample_info : pandas.DataFrame
            Dataframe contains sample information: name, i7_index_id, index, i5_index_id, index2, sample_date.
        """
        sample_sheet = SampleSheet(sample_sheet_fname)
        instrument_type = sample_sheet.header.instrument_type
        sample_data = []
        for sample in sample_sheet.samples:    
            sampData = collections.OrderedDict({
                "sample_name" : sample.sample_name,
                "sample_id" : sample.sample_id,                
                "i7_index_id" : sample.i7_index_id,
                "i7_index" : sample.index,
                "lane" : None,
                "i5_index_id" : sample.i5_index_id,
                "i5_index" : sample.index2
            })
            sample_data.append(sampData)
        sample_info = pd.DataFrame(sample_data)
        return sample_info

        
    def prepare_records(self):
        """Gather information and results to generate records that contain CDR3 detected in a single sequencing run.

        Gather sample information from sample sheet. Gather sequencing run information from RunInfo.xml in raw
        sequencing run folder. Gather lane information from FASTQ header. Gather CDR3 results from pipeline run.

        Returns
        -------
        records : pandas.DataFrame
            Dataframe contains CDR3 records.
        """
        
        self.parse_manifest(self.sample_manifest_fname)
        self.logger.info('Gather records.')
        
        analysis_info = self.get_analysis_info() 
        sequencing_info = self.get_sequencing_info(self.seq_run_path)
        param_info = self.get_param_info() 
        sample_info = self.get_sample_info(self.sample_sheet_fname)
        
        records = {"sequencing_info" : sequencing_info,
                   "sample_info" : sample_info,
                   "analysis_info" : analysis_info,
                   "param_info" : param_info}
        
        return records

    
    def insert_records(self, overwrite=False):
        """Insert CDR3 records to the database.

        Parameters
        ----------
        overwrite : bool
            Whether to overwrite the records in the database if the sequencing run is already stored in the database.
            This will delete the records of existing run and insert the records provided.
        """
        self.logger.info('Insert CDR3 records to the database.')
        records = self.records.copy()
        # Only record D99 families. Most of CDR3 sequences don't pass D99 are singletons.
        records = records.loc[records['cummulative_fraction'] < 0.99, :]
        records.drop(columns='cummulative_fraction', inplace=True)
        self.daedalus_db.insert(records, overwrite)
        

    def parse_manifest(self, sample_manifest_fname):
        manifest = pd.read_csv(sample_manifest_fname)
        self.pipe_run_path = os.path.abspath(os.path.expanduser(os.path.dirname(sample_manifest_fname)))   
        self.sample_manifest = manifest
        self.sample_sheet_fname = manifest['sample_sheet'][0]
        self.seq_run_path = manifest['run_folder'][0]

    def setParamValues(self, params, updateValues):
        for k,v in updateValues.items():
            if k in params:
                params[k] = v
            else: ##create new values
                params[k] = v
        return params
        

    def get_lane_info(self, sample_info, pipe_run_path):
        """Get sequencing lane from FASTQ file.

        Parameters
        ----------
        sample_info : pd.DataFrame
            Dataframe contains sample information.
        pipe_run_path : str
            Path to the immunoPETE pipeline run results. The folder should contain FASTQ files for each sample under
            `analysis/fastq`.

        Returns
        -------
        lane_info : pandas.DataFrame
            Dataframe contains lane information derived from FASTQ header: instrument, run_number, flowcell_id, name,
           ` lane. Most of these are redundant to the run_info, except  `lane`.
        """
        lane_info = defaultdict(list)
        for _, row in sample_info.iterrows():
            fastq_path_pattern = os.path.join(
                pipe_run_path,
                'analysis',
                'fastq',
                row['sample_name'] + '_S[0-9]*_R1_001.fastq.gz'
            )

            fastq_path = glob.glob(fastq_path_pattern)
            if len(fastq_path) == 1:
                fastq_path = fastq_path[0]
            else:
                if len(fastq_path) > 1:
                    sys.stderr.write('More than 1 fastq file match given sample name: {}.\n'.format(fastq_path))
                else:
                    sys.stderr.write('No fastq file match given sample name: {}.\n'.format(fastq_path_pattern))
                sys.exit(1)

            with gzip.open(fastq_path, 'rt') as f:
                header = f.readline()
                instrument, run_number, flowcell_id, lane, *_ = header.lstrip('@').split(':')
                lane_info['instrument'].append(instrument)
                lane_info['run_number'].append(int(run_number))
                lane_info['flowcell_id'].append(flowcell_id)
                lane_info['lane'].append(int(lane))
                lane_info['sample_name'].append(row['sample_name'])

        lane_info = pd.DataFrame(lane_info)
        return lane_info



def main():
    parser = argparse.ArgumentParser(description="""Given sample sheet and the pipeline run path, insert the CDR3
        results to the database. And mark possible CDR3 cross-sample contamination.""")
    parser.add_argument('--runs_to_keep', type=int, default=50,
        help='Number of recent sequencing runs to keep in the database. Older runs will be deleted from the database.')
    parser.add_argument('--url', default='sqlite:////sc1/groups/pls-redbfx/databases/daedalus_db/daedalus.db',
                    help='URL to the Daedalus sample database.')
    parser.add_argument('--new_database',
                        default=False,
                        action='store_true',
                        help='This is a new database, create the db file from scratch.')
    parser.add_argument('sample_manifest',
                        help='Path to the sample_manifest file, for initiating pipeline runs.')
    args = parser.parse_args()

    logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(name)s - %(message)s', level=logging.INFO)
    logger = logging.getLogger(__name__)

    logger.info('Detect cross-sample contamination.')
    pipelineLogger = PipelineLogger(args.sample_manifest, args.url)
    pipelineLogger.prepare_records()
    if args.mode == 'cross-run':
        pipelineLogger.detect(
            update_db=args.update_db,
            cross_run=True,
            num_runs_detect=args.num_runs_detect,
            runs_to_keep=args.runs_to_keep)
    else:
        pipelineLogger.detect(args.out_fname, update_db=args.update_db, cross_run=False)
    logger.info('Finished.')

if __name__ == '__main__':
    main()
