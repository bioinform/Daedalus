#!/usr/bin/env python
import argparse
import logging
import os
import shutil
import subprocess
import sys
import re
import collections
import time

from manifest_parser import ManifestParser
from manifest_parser import NextflowConfig
from default_params import *

def wait_for_memory(mem_required_gb=12):
    """Wait for enough memory to run the function.
    """
    logger = logging.getLogger(__name__)
    def decorator(func):
        def wrapper(*args, **kwargs):
            cmd = 'free -g'
            curr_free_mem_gb = 0
            while True:
                try:
                    curr_free_mem_gb = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True).decode('ascii').split('\n')[1].split()[6]
                except OSError as e:
                    pass
                    logger.warning('WARNING: OSError occured during resource checking: {}'.format(str(e)))
                if float(curr_free_mem_gb) >= float(mem_required_gb):
                    result = func(*args, **kwargs)
                    return result
                logger.info('Waiting for the memory. Request: {}G. Current available memory: {}G'.format(mem_required_gb, curr_free_mem_gb))
                time.sleep(60)
        return wrapper
    return decorator


class PipelineRunner:
    """A runner that runs Nextflow pipeline based on the manifest file.

    Attributes
    ----------
    repo_path : str
        Path to the immunoPETE repo.
    pipeline_path : str
        Path to the Nextflow pipeline script.
    html_generator_path : str
        Path to the html_generator.py.
    pipeline_reporter_path : str
        Path to the pipeline_summary.py.
    group : str
        Group to submit jobs on SC1.
    genes : set
        Genes to detect: {'TRA', 'TRB', 'TRD', 'TRG', 'IGH', 'IGL', 'IGK'}.
    idb_type : set
        IDB type to detect: {'gene', 'pseudogene', 'pseudo-CDR3'}. Set to None will include all.
    local : bool
        If true, run the pipeline in local machine.
    no_fairshare : bool
        If true, submit jobs without Fair Share policy.
    """

    def __init__(self, group, genes, idb_type=None, local=False, no_fairshare=False):
        self.repo_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        self.segment1_path = os.path.join(self.repo_path,
                                          'nextflow',
                                          'segment1-ipete.nf')
        self.segment2_path = os.path.join(self.repo_path,
                                          'nextflow',
                                          'segment2-ipete.nf')
        self.segment3_path = os.path.join(self.repo_path,
                                          'nextflow',
                                          'segment3-ipete.nf')
        self.group = group
        self.genes = genes
        self.idb_type = idb_type
        self.local = local
        self.no_fairshare = no_fairshare
        self.logger = logging.getLogger(__name__)

    def match_files(self, path, pattern):
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

        
    def build_seg2_configs(self, seg1Output, sample_info):
        configs = collections.OrderedDict()
        for sample in sample_info["sample_name"]:            
            R1files = self.match_files(seg1Output, "{}_S[0-9]+_R1_001.fastq.gz".format(sample))
            R2files = self.match_files(seg1Output, "{}_S[0-9]+_R2_001.fastq.gz".format(sample))            
            sampleParams = sample_info.loc[sample_info["sample_name"] == sample]            
            sampleParams = sampleParams.reset_index()
            if len(R1files) == 1 and len(R2files) == 1:                
                sampleParams["fastq1"] = R1files[0]
                sampleParams["fastq2"] = R2files[0]
                sampleParams["sample"] = sample
                configs[sample] = NextflowConfig()
                for column in sampleParams.columns:
                    configs[sample][column] = sampleParams[column][0]
            else:
                continue
                    
        return configs
            
    def run(self, fname, out_dir, resume, prefix=None, wait=False):
        """Run the Nextflow pipeline based on the manifest file.

        Parameters
        ----------
        fname : str
            Path to the manifest file.
        out_dir : str
            Parent directory of the Nextflow pipeline run. The Nextflow results will be stored at `out_dir/project`,
            where the `project` is inferred from the manifest file.
        resume : bool
            When Nextflow work directory exists, resume the pipeline.
        prefix : str
            Prefix project name with specified prefix. If None, the current date will be used.
        wait : bool
            Wait for the pipeline to finish.

        Returns
        -------
        None
        """
        out_dir = os.path.realpath(out_dir)
        fname = os.path.abspath(os.path.expanduser(fname))
        mp = ManifestParser(fname, out_dir, prefix)
        
        if mp.is_valid():
            for (project, config), (_, sample_info) in zip(mp.configs.items(), mp.sample_info.items()):                
                work_dir = os.path.dirname(os.path.abspath(config['finalDir']))
                manifest_dir = os.path.join(out_dir, 'manifest')
                manifest_path = os.path.join(manifest_dir, project + '_manifest.csv')
                config_path = os.path.join(work_dir, 'config.txt')
                config_dir = os.path.join(work_dir, 'config')
                seg1_work_dir = os.path.join(config_dir, 'seg1')
                seg3_work_dir = os.path.join(config_dir, 'seg3')
                sample_info_path = os.path.join(work_dir, 'sample_info.csv')
                nextflow_log_path = os.path.join(work_dir, 'log')

                if not os.path.exists(manifest_dir):
                    os.makedirs(manifest_dir)
                    os.chmod(manifest_dir, mode=0o777)
                
                if not os.path.exists(work_dir) or resume:
                    try:
                        os.makedirs(work_dir)
                        os.makedirs(os.path.join(work_dir, "analysis"))
                        os.makedirs(seg1_work_dir, exist_ok=True)
                        os.makedirs(seg3_work_dir, exist_ok=True)
                        os.makedirs(nextflow_log_path, exist_ok=True)
                    except FileExistsError:
                        self.logger.info(
                            'Folder {} exists, resuming the pipeline.'.format(work_dir))
                    shutil.copyfile(fname, manifest_path)
                    config['pipeline_summary'] = getDefaultParams()["pipeline_summary"]
                    config.write(config_path)
                    sample_info.to_csv(sample_info_path, index=False)
                    # Generate reference in the work directory
                    self._run_segment1(seg1_work_dir, config, config_path,
                                       sample_info_path, os.path.join(nextflow_log_path, 'seg1.log'))
                    seg1Output = os.path.join(work_dir, "analysis", "fastq")
                    seg2Configs = self.build_seg2_configs(seg1Output, sample_info)
                    seg2_procs = {}
                    seg2_log = {}
                    for sample in seg2Configs:
                        sample_dir = os.path.join(work_dir, "analysis", sample)
                        seg2_work_dir = os.path.join(config_dir, 'seg2_{}'.format(sample))
                        if not os.path.exists(sample_dir):
                            try:
                                os.makedirs(sample_dir)
                                os.makedirs(seg2_work_dir, exist_ok=True)
                            except FileExistsError:
                                self.logger.info(
                                    'Folder {} exists, resuming the pipeline.'.format(work_dir))
                        sample_config_path = os.path.join(sample_dir, 'config.txt')
                        sampleConfig = seg2Configs[sample]
                        sampleConfig['finalDir'] = config['finalDir']
                        sampleConfig.write(sample_config_path)
                        if self.local:
                            proc = self._run_segment2(seg2_work_dir, sampleConfig, sample_config_path,
                                           sample_info_path, os.path.join(nextflow_log_path, 'seg2_{}.log'.format(sample)))
                        else:
                            proc = self._submit_segment2(seg2_work_dir, sampleConfig, sample_config_path,
                                            sample_info_path, os.path.join(nextflow_log_path, 'seg2_{}.log'.format(sample)))
                        seg2_procs[sample] = proc
                        time.sleep(5)
                
                    for sample, proc in seg2_procs.items():
                        exit_code = proc.wait()
                        status = 'success' if exit_code == 0 else 'fail'
                        
                        seg2_log[sample] = status
                        self.logger.info('Finished {sample}: {status}'.format(sample=sample, status=status))

                    ##wait for files to complete transfering from seg2
                    time.sleep(10)
                    self._run_segment3(seg3_work_dir, config, config_path,
                                       os.path.join(nextflow_log_path, 'seg3.log'))

                    
                else:
                    sys.stderr.write('Directory exists: {}.\n'.format(work_dir))
                    sys.stderr.write('To resume the unfinished pipeline, use "-r/--resume".\n')
                    return

    @wait_for_memory()
    def _run_segment1(self, work_dir, config, config_path, sample_info_path, nextflow_log_path):
        """Run the Nextflow command.
        Parameters
        ----------
        work_dir : str
            Path of the working directory where the Nextflow command run from.
        config : NextflowConfig
            Object contains configs parsed from manifest.
        config_path : str
            Path of the Nextflow config file.
        sample_info_path : str
            Path of the sample_info.csv.
        nextflow_log_path : str
            Path of the log file generated by Nextflow.

        Returns
        -------
        None
        """
        os.chdir(work_dir)

        profile = 'ipete_docker' if not self.local else 'local'

        command = '''
        export NXF_OPTS="-Xmx512M";
        nextflow run {segment1_path} -resume -profile {profile} -c {config_path} -ansi-log false &> {log}
        '''.format(segment1_path=self.segment1_path,
                   profile=profile,
                   config_path=config_path,
                   log=nextflow_log_path)
        self.logger.info('Submitting Segment1.')
        subprocess.run(command, shell=True)
        self.logger.info('Pipeline is finished. Check results at {}'.format(nextflow_log_path))

    def _submit_segment2(self, work_dir, config, config_path, sample_info_path, nextflow_log_path):
        os.chdir(work_dir)
        profile = 'ipete_docker' if not self.local else 'local'

        script = '''
        export NXF_OPTS="-Xmx512M";
        nextflow run {segment2_path} -resume -profile {profile} -c {config_path} -ansi-log false &> {log}
        '''.format(segment2_path=self.segment2_path,
                   profile=profile,
                   config_path=config_path,
                   log=nextflow_log_path)
        with open('run_workflow.sh', 'w') as f:
            f.write(script)
            f.flush()
        if self.no_fairshare:
            command = f'qsub -V -S /bin/bash -cwd -l h_vmem=12G -N nf-submit -sync y run_workflow.sh'
        else:
            command = f'qsub -V -S /bin/bash -cwd -l h_vmem=12G -N nf-submit -P {self.group} -sync y run_workflow.sh'

        self.logger.info('Submitting Segment2.')
        proc = subprocess.Popen(command, shell=True)
        self.logger.info('Pipeline is submitted. Check the status at {}'.format(nextflow_log_path))
        return proc


    @wait_for_memory()
    def _run_segment2(self, work_dir, config, config_path, sample_info_path, nextflow_log_path):
        """Run the Nextflow command.
        Parameters
        ----------
        work_dir : str
            Path of the working directory where the Nextflow command run from.
        config : NextflowConfig
            Object contains configs parsed from manifest.
        config_path : str
            Path of the Nextflow config file.
        sample_info_path : str
            Path of the sample_info.csv.
        nextflow_log_path : str
            Path of the log file generated by Nextflow.

        Returns
        -------
        proc : Popen
            A Popen object returned by subprocess.Popen.
        """
        os.chdir(work_dir)

        profile = 'ipete_docker' if not self.local else 'local'

        command = '''
        export NXF_OPTS="-Xmx512M";
        nextflow run {segment2_path} -resume -profile {profile} -c {config_path} -ansi-log false &> {log}
        '''.format(segment2_path=self.segment2_path,
                   profile=profile,
                   config_path=config_path,
                   log=nextflow_log_path)

        self.logger.info('Submitting Segment2.')
        proc = subprocess.Popen(command, shell=True)
        self.logger.info('Pipeline is submitted. Check the status at {}'.format(nextflow_log_path))
        return proc

    @wait_for_memory()
    def _run_segment3(self, work_dir, config, config_path, nextflow_log_path):
        """Run the Nextflow command.
        Parameters
        ----------
        work_dir : str
            Path of the working directory where the Nextflow command run from.
        config : NextflowConfig
            Object contains configs parsed from manifest.
        config_path : str
            Path of the Nextflow config file.
        nextflow_log_path : str
            Path of the log file generated by Nextflow.

        Returns
        -------
        proc : Popen
            A Popen object returned by subprocess.Popen.
        """
        os.chdir(work_dir)

        profile = 'ipete_docker' if not self.local else 'local'

        command = '''
        export NXF_OPTS="-Xmx512M";
        nextflow run {segment3_path} -resume -profile {profile} -c {config_path} -ansi-log false &> {log}
        '''.format(segment3_path=self.segment3_path,
                   profile=profile,
                   config_path=config_path,
                   log=nextflow_log_path)

        self.logger.info('Submitting Segment3.')
        ##subprocess.run(command, shell=True)
        proc = subprocess.Popen(command, shell=True)
        proc.wait()
        self.logger.info('Pipeline is submitted. Check the status at {}'.format(nextflow_log_path))

    def get_version(self):
        """Get pipeline version using git.

        Returns
        -------
        version : str or None
            The tag string if the pipeline contains a tag. Otherwise it's the `branch`-`shorthash`.
        """
        current_dir = os.getcwd()
        pipeline_dir = os.path.dirname(self.pipeline_path)
        os.chdir(pipeline_dir)
        try:
            version = subprocess.run('git describe --tags', shell=True,
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
            version = version.stdout.strip().decode("utf-8")
        except subprocess.CalledProcessError:
            try:
                version = subprocess.run('git rev-parse --abbrev-ref HEAD', shell=True,
                                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
                version = version.stdout.strip().decode("utf-8")
                hashtag = subprocess.run('git rev-parse --short HEAD', shell=True,
                                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
                hashtag = hashtag.stdout.strip().decode("utf-8")
                version = version + '-' + hashtag
            except subprocess.CalledProcessError:
                self.logger.error('No git information at {}'.format(pipeline_dir))
                return

        os.chdir(current_dir)
        return version


def parse_args():
    parser = argparse.ArgumentParser(description="""Run Nextflow pipeline.""")
    parser.add_argument('--chain', type=str, default='A,B,G,D,H,K,L',
                        help='Chain type. A,B,G,D for TR; H,K,L for IG. Use comma to separate multiple chains. Default is A,B,G,D,H,K,L')
    parser.add_argument('-r', '--resume', action='store_true',
                        help='When Nextflow work directory exists, resume the pipeline.')
    parser.add_argument('-g', '--group', default='rssprbf',
                        help='Group to submit jobs on SC1.')
    parser.add_argument('--gene_type', type=str, default='ALL',
                        help='Gene type to detect: IG, TR, ALL. Default is ALL.')
    parser.add_argument('--idb_type', type=str,
                        help='IDB type to include: gene, noCDR3. Use comma to separate multiple types. Default will include all.')
    parser.add_argument('--local', action='store_true',
                        help='Run pipeline locally.')
    parser.add_argument('--no_fairshare', action='store_true',
                        help='Submit jobs without Fair Share policy.')
    parser.add_argument('--wait', action='store_true',
                        help='Wait for the pipeline to finish. Without this option, nextflow will run in background.')
    parser.add_argument('-o', '--out_dir', default='/sc1/groups/pls-redbfx/pipeline_runs/iPETE',
                        help='Nextflow output directory.')
    parser.add_argument('--prefix',
                        help='Prefix project name with specified prefix. If None, the current date will be used.')
    parser.add_argument('manifest',
                        help='Manifest file name.')
    args = parser.parse_args()
    chains = set(args.chain.strip(',').upper().split(','))
    gene_types = set(args.gene_type.strip(',').upper().split(','))
    tr_chains = set(['A', 'B', 'D', 'G'])
    ig_chains = set(['H', 'L', 'K'])
    assert chains.issubset(tr_chains.union(ig_chains))
    assert gene_types.issubset(['IG', 'TR', 'ALL'])
    genes = set()
    if 'IG' in gene_types:
        chains_to_detect = chains.intersection(ig_chains)
        if len(chains_to_detect) == 0:
            parser.error('Please specifiy chains {} for IG.'.format(ig_chains))
        for c in chains_to_detect:
            genes.add('IG'+c)
    if 'TR' in gene_types:
        chains_to_detect = chains.intersection(tr_chains)
        if len(chains_to_detect) == 0:
            parser.error('Please specifiy chains {} for TR.'.format(tr_chains))
        for c in chains_to_detect:
            genes.add('TR'+c)
    args.genes = genes
    if args.idb_type:
        args.idb_type = set(args.idb_type.strip(',').split(','))
    return args


def main(args):
    logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(name)s - %(message)s', level=logging.INFO)
    logger = logging.getLogger(__name__)

    logger.info('Start.')
    pr = PipelineRunner(
        group=args.group,
        genes=args.genes,
        idb_type=args.idb_type,
        local=args.local,
        no_fairshare=args.no_fairshare)
    pr.run(args.manifest, args.out_dir, args.resume, args.prefix, args.wait)


if __name__ == '__main__':
    args = parse_args()
    main(args)
