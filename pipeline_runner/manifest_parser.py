#!/usr/bin/env python
import datetime
import logging
import os
import pandas as pd
import numpy as np
import math

from collections import OrderedDict
from itertools import chain

class ManifestParser:
    """A parser that read csv/txt manifest and convert it to Nextflow config file.

    Attributes
    ----------
    out_dir : str
        Parent directory of the Nextflow pipeline run.
    manifest : pandas.DataFame
        Dataframe contains at least three columns: `read1`, `read2`, `project`.
    configs : collections.OrderedDict
        Dict of key-value pairs where each key is a project name, each value is another dict contains key-value pairs of
        configurations converted from the manifest file.
    sample_info : collections.OrderedDict
        Dict of key-value pairs where each key is a project name, each value is a pandas.DataFrame of project specific
        manifest.
    prefix : str
        Prefix project name with specified prefix. If None, the current date will be used.
    """

    def __init__(self, fname, out_dir, prefix=None):
        self.logger = logging.getLogger(__name__)
        self.required_cols = ('project',)
        self.out_dir = out_dir
        self.prefix = prefix
        self.manifest = self._read(fname)
        self._is_manifest_valid = self._validate_manifest()
        self.configs, self.sample_info = self._convert()
        self._is_config_valid = self._validate_config()
        
    def _convert(self):
        """Convert manifest to Nextflow config file.

        Returns
        -------
        configs : collections.OrderedDict
            Dict of key-value pairs where each key is a project name, each value is another dict contains key-value pairs of
            configurations converted from the manifest file.
        sample_info : collections.OrderedDict
            Dict of key-value pairs where each key is a project name, each value is a pandas.DataFrame of project specific
            manifest.
        """
        if not self._is_manifest_valid:
            return
        self.logger.info('Parsing manifest.')

        configs = OrderedDict()
        sample_info = OrderedDict()

        for project, manifest in self.manifest.groupby('project'):
            manifest = manifest.copy()
            
            if not self.prefix:
                prefix = str(datetime.date.today())
            else:
                prefix = self.prefix
            project = prefix + '-' + project.strip('/').replace(' ', '_')

            configs[project] = NextflowConfig()
            configs[project]['project'] = project
            configs[project]['finalDir'] = os.path.join(self.out_dir, project, 'analysis')
            configs[project]['analysis_dir'] = os.path.join(self.out_dir, project)
            configs[project]['run_dir'] = manifest["run_folder"][0]
            configs[project]['sample_sheet'] = manifest["sample_sheet"][0]            
            configs[project]['runBcl2Fastq'] = 'true'
            manifest['project'] = project
            sample_info[project] = manifest
            
        return configs, sample_info

    def _read(self, fname):
        """Read manifest file.

        The manifest should contain at least three columns: `read1`, `read2`, `project`. It can contain comments and/or
        descriptions start with `#`.

        Parameters
        ----------
        fname : str
            Path to the manifest file.

        Returns
        -------
        manifest : pandas.DataFrame or None
            DataFrame contains manifest information without comments.
        """
        self.logger.info('Read manifest file {}.'.format(fname))

        sep = None
        _, extension = os.path.splitext(fname)
        
        if extension == '.csv':
            sep = ','
        elif extension == '.txt' or extension == '.tsv':
            sep = '\t'
        else:
            self.logger.error('Unrecognized file extension: {}. Please provide a comma-delimited ".csv" or a tab-delimited ".txt/.tsv" manifest file.'.format(fname))
            return

        manifest = None
        try:
            manifest = pd.read_csv(fname, sep=sep, comment='#')
        except FileNotFoundError as e:
            self.logger.error(e)
            return
    
         # for column, colType in zip(manifest.columns, manifest.dtypes):
         #    #if colType == "bool":
         #    #    manifest[column] = manifest[column].map({True: 'true', False: 'false'})
         #    vals = []
         #    hasNa = False
         #    for val in manifest[column]:
         #        try:
         #            if math.isnan(val):
         #                hasNa = True
         #                vals.append("")
         #            else:
         #                vals.append(val)
         #        except:
         #            vals.append(val)
         #    if hasNa == True:
         #        manifest[column] = vals
            
        return manifest

    def is_valid(self):
        """Validate the manifest and config files.

        Returns
        -------
        is_valid : bool
            Boolean indicates whether the stored manifest and config are valid or not.
        """
        return self._is_manifest_valid and self._is_config_valid

    def _validate_manifest(self):
        """Validate stored manifest file.
        """
        # Check manifest is loaded and contains required columns.
        if self.manifest is None:
            return False
        else:
            for col in self.required_cols:
                if col not in self.manifest.columns:
                    self.logger.error('Missing required column: {}.'.format(col))
                    return False
        return True

    def _validate_config(self):
        """Validate the stored configs.
        """
        if not self.configs or not self.sample_info:
            return False
        return True

class NextflowConfig(OrderedDict):
    """A collections.OrderedDict contains key-value pairs of Nextflow config.
    """
    def __init__self(self, *args, **kwargs):
        super(NextflowConfig, self).__init__(*args, **kwargs)
        self.logger = logging.getLogger(__name__)

    def write(self, fname):
        """Write Nextflow config file parsed from manifest.

        Parameters
        ----------
        fname : str
            Path to the output Nextflow config file.

        Returns
        -------
        None
        """
        if self.__len__() != 0:
            with open(fname, 'w') as f:
                f.write(self.__str__())
        else:
            self.logger.error('Config is empty.')

    def __str__(self):
        config_str = ''
        for k, v in self.items():
            if isinstance(v, str):
                v = '"{}"'.format(v)
            if isinstance(v, bool):
                v = 'true' if v else 'false'
            try:
                if np.issubdtype(v, np.bool_):
                    v = 'true' if v else 'false'
            except:
                pass
            try:
                if math.isnan(v):
                    v = ""
            except:
                pass            
            config_str += 'params.{} = {}\n'.format(k, v)
        return config_str
