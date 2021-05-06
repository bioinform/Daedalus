import os
import logging
import pandas as pd

class ReferenceGenerator:
    """Generate references based on immunoDB.

    Attributes
    ----------
    v : pandas.DataFrame
        A dataframe contains annotated V sequence records from immunoDB.
    j : pandas.DataFrame
        A dataframe contains annotated J sequence records from immunoDB.
    """
    def __init__(self, v_fname, j_fname, vprim_fname, jprim_fname):
        self.v = self._load(v_fname)
        self.j = self._load(j_fname)
        self.vprim = self._load(vprim_fname)
        self.jprim = self._load(jprim_fname)
        self.logger = logging.getLogger(__name__)

    def export_fasta(self, geneDF, out_file):
        """
        Write fasta file from immunoDB, pandas dataframe.
        
        Parameters
        ----------
        geneDF : pandas DF
        data frame of gene targets
        out_file : str
        path to output results
        """
        fa_out = open(out_file, 'w')
        for seqId, seq in zip(geneDF["sequence_id"],
                              geneDF["sequence"]):
            fa_out.writelines(">{}\n{}\n".format(seqId, seq))
        fa_out.close()
                              
        
    def extract(self, genes, out_dir, idb_type=None):
        """Extract annotated sequence records from immunoDB based on given gene types.

        Parameters
        ----------
        genes : set
            Genes to detect: {'TRA', 'TRB', 'TRD', 'TRA', 'IGH', 'IGL', 'IGK'}.
        out_dir : str
            Path to the output directory.
        idb_type : set
            IDB type to include.
        """
        gene_pattern = '^' + '|^'.join(genes)
        v = self.v.loc[self.v['symbol'].str.contains(gene_pattern),:]
        j = self.j.loc[self.j['symbol'].str.contains(gene_pattern),:]
        vprim = self.vprim.loc[self.vprim['sequence_id'].str.contains(gene_pattern),:]
        jprim = self.jprim.loc[self.jprim['sequence_id'].str.contains(gene_pattern),:]
        
        if idb_type:
            v = v.loc[v['IDB_type'].isin(idb_type),:]
            j = j.loc[j['IDB_type'].isin(idb_type),:]

        try:
            os.makedirs(out_dir)
        except FileExistsError:
            self.logger.info(
                'Folder {} exists, overwrite.'.format(out_dir))
            
        self.export_fasta(v, os.path.join(out_dir, 'Vgene_targets.fasta'))
        self.export_fasta(j, os.path.join(out_dir, 'Jgene_targets.fasta'))
        self.export_fasta(vprim, os.path.join(out_dir, 'Vgene_primers.fasta'))
        self.export_fasta(jprim, os.path.join(out_dir, 'Jgene_primers.fasta'))

        vj = pd.concat([v, j], axis=0)
        vj.to_csv(os.path.join(out_dir, 'immunoDB_VJtargets.csv'), index=False)
        vjprim = pd.concat([vprim, jprim], axis=0)
        vjprim.to_csv(os.path.join(out_dir, 'immunoDB_VJprimers.csv'), index=False)

        
    def _load(self, fname):
        """Load annotated sequences from immunoDB.

        Parameters
        ----------
        fname : str
            Path to the file that contains annotated sequences.

        Returns
        -------
        reference : pandas.DataFrame
            A dataframe contains annotated sequence records from immunoDB.
        """
        reference = pd.read_csv(fname)
        return reference
