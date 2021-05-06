import logging
import warnings

import pandas as pd

from contextlib import contextmanager
from sqlite3 import Connection as SQLite3Connection

from sqlalchemy import create_engine, event, and_
from sqlalchemy import exc as sa_exc
from sqlalchemy.engine import Engine
from sqlalchemy.orm import sessionmaker

from daedalus_db.models import Base, Analysis, SequencingRun, Sample, IpeteParams #, SampleSummary

class DaedalusDB:
    """Connector to Daedalus Pipeline database.

    Attributes
    ----------
    url : str
        URL to the database.
    engine : object
        SQLAlchemy database engine.
    Session : object
        SQLAlchemy Session.
    """
    def __init__(self, url='sqlite:///daedalus.db', echo=False, newDb=False, timeout=60):
        self.url = url
        self.engine = create_engine(self.url, echo=echo, connect_args={'timeout': timeout})
        self.Session = sessionmaker(bind=self.engine)
        self.logger = logging.getLogger(__name__)
        if newDb == True:
            self.logger.info("creating new database")
            self.create_db()
            
    def checkRunInfo(self, run_info, overwrite=False):
        with self.session_scope() as session:
            self.logger.info('Inserting records to the database.')
            instrument, run_number, flowcell_id, sequencing_date = (run_info["instrument"][0],
                                                                    run_info["run_number"][0],
                                                                    run_info["flowcell_id"][0],
                                                                    run_info["sequencing_date"][0])
                                                                       
            existing_run = session.query(
                DaedalusDB
            ).filter(and_(
                DaedalusDB.instrument==instrument,
                DaedalusDB.run_number==run_number,
                DaedalusDB.flowcell_id==flowcell_id,
                DaedalusDB.sequencing_date==sequencing_date)
            ).first()
            if not existing_run or overwrite:
                session.query(
                    DaedalusDB
                ).filter(and_(
                    DaedalusDB.instrument==instrument,
                    DaedalusDB.run_number==run_number,
                    DaedalusDB.flowcell_id==flowcell_id,
                    DaedalusDB.sequencing_date==sequencing_date)
                ).delete()
                session.bulk_insert_mappings(DaedalusDB, run_info.to_dict('records'))
                
        
    def insert_analysis(self, records, overwrite=False):
        """Insert CDR3 records from a sequencing run to the database.

        Parameters
        ----------
        records : pandas.DataFrame
            Dataframe contains instrument, run_number, flowcell_id, sequencing_date, sample_name, lane, i7_index_id,
             index, i5_index_id, index2, cdr3_nt, cdr3_aa, cdr3_qual, v, d, j, umi_seq, umi_qual, family_size.
        overwrite : bool
            Whether to overwrite the records in the database if the sequencing run is already stored in the database.
            This will delete the records of existing run and insert the records provided.
        """
        
        with self.session_scope() as session:
            sequencing_info = records['sequencing_info']
            sample_info = records['sample_info']
            analysis_info = records['analysis_info']
            param_info = records['param_info']
            self.logger.info('Inserting records to the database.')

            assert sequencing_info.shape == (1, 7)
            (sequencing_platform, instrument, run_number, flowcell_id,
                sequencing_date, sample_sheet, run_folder) = sequencing_info.iloc[0]
            sequencing_run_query = session.query(
                                SequencingRun
                            ).filter(and_(
                                SequencingRun.instrument==instrument,
                                SequencingRun.run_number==run_number,
                                SequencingRun.flowcell_id==flowcell_id,
                                SequencingRun.sequencing_date==sequencing_date)
                            ).one_or_none()
            if not sequencing_run_query or overwrite:
                session.query(
                    SequencingRun
                ).filter(and_(
                    SequencingRun.instrument==instrument,
                    SequencingRun.run_number==run_number,
                    SequencingRun.flowcell_id==flowcell_id,
                    SequencingRun.sequencing_date==sequencing_date)
                ).delete()
                self.logger.info('Overwrite sequencing run.')
                sequencing_run_query = SequencingRun(**sequencing_info.iloc[0].to_dict())
                session.add(sequencing_run_query)

            for p in param_info.drop_duplicates().to_dict('records'):
                param_query = session.query(
                    IpeteParams
                ).filter(and_(
                    IpeteParams.subsample == p['subsample'],
                    IpeteParams.subsampleSeed == p['subsampleSeed'],
                    IpeteParams.primer_targets == p['primer_targets'],
                    IpeteParams.inputType == p['inputType'],
                    IpeteParams.trim_primers == p['trim_primers'],
                    IpeteParams.checkSpikein == p['checkSpikein'],
                    IpeteParams.umi_mode == p['umi_mode'],
                    IpeteParams.umi1 == (p['umi1'] if not pd.isna(p['umi1']) else None),
                    IpeteParams.umi2 == (p['umi2'] if not pd.isna(p['umi2']) else None),
                    IpeteParams.umi_type == p['umi_type'],
                    IpeteParams.vPrimerRef == p['vPrimerRef'],
                    IpeteParams.jPrimerRef == p['jPrimerRef'],
                    IpeteParams.vGeneRef == p['vGeneRef'],
                    IpeteParams.jGeneRef == p['jGeneRef'],
                    IpeteParams.primerRef == p['primerRef'],
                    IpeteParams.geneRef == p['geneRef'],
                    IpeteParams.spikeRef == p['spikeRef'],
                    IpeteParams.vGeneScore == p['vGeneScore'],
                    IpeteParams.vGeneKmer == p['vGeneKmer'],
                    IpeteParams.vGeneFrac == p['vGeneFrac'],
                    IpeteParams.jGeneScore == p['jGeneScore'],
                    IpeteParams.jGeneKmer == p['jGeneKmer'],
                    IpeteParams.jGeneFrac == p['jGeneFrac'],
                    IpeteParams.vPrimScore == p['vPrimScore'],
                    IpeteParams.vPrimKmer == p['vPrimKmer'],
                    IpeteParams.vPrimFrac == p['vPrimFrac'],
                    IpeteParams.jPrimScore == p['jPrimScore'],
                    IpeteParams.jPrimKmer == p['jPrimKmer'],
                    IpeteParams.jPrimFrac == p['jPrimFrac'],
                    IpeteParams.spikeScore == p['spikeScore'],
                    IpeteParams.spikeKmer == p['spikeKmer'],
                    IpeteParams.spikeFrac == p['spikeFrac'],
                    IpeteParams.min_read_qual == p['min_read_qual'],
                    IpeteParams.minGeneIdentity == p['minGeneIdentity'],
                    IpeteParams.max_cdr3_length == p['max_cdr3_length'],
                    IpeteParams.bidding_ratio == p['bidding_ratio'],
                    IpeteParams.max_steps == p['max_steps'],
                    IpeteParams.umi_edit_dist == p['umi_edit_dist'],
                    IpeteParams.cdr3_edit_dist == p['cdr3_edit_dist']
                )).one_or_none()
                if not param_query or overwrite:
                    if param_query:
                        self.logger.info('Overwrite params.')
                        session.delete(param_query)
                    param_query = IpeteParams(**p)
                    session.add(param_query)

            for s in records['sample_info'].to_dict('records'):
                sample_query = session.query(
                    Sample
                ).filter(and_(
                    Sample.sample_name == s['sample_name'],
                    Sample.run_id == sequencing_run_query.id,
                    Sample.lane == s['lane']
                )).one_or_none()
                if not sample_query or overwrite:
                    if sample_query:
                        self.logger.info('Overwrite sample.')
                        session.delete(sample_query)
                    sample_query = Sample(**s)
                    sample_query.sequencing_run = sequencing_run_query

            (analysis_name, analysis_folder, analysis_date, pipeline_version,
                segment_1_status, segment_2_status, segment_3_status) = records['analysis_info'].iloc[0]
            analysis_query = session.query(
                Analysis
            ).filter(and_(
                Analysis.analysis_name == analysis_name,
                Analysis.pipeline_version == pipeline_version,
                Analysis.seq_id == sequencing_run_query.id,
                Analysis.param_id == param_query.id
            )).one_or_none()
            if not analysis_query or overwrite:
                if analysis_query:
                    self.logger.info('Overwrite analysis.')
                    session.delete(analysis_query)
                analysis_query = Analysis(**records['analysis_info'].iloc[0].to_dict())
                analysis_query.seq_id = sequencing_run_query.id
                analysis_query.param_id = param_query.id
                self.logger.info('Insert analysis.')
                session.add(analysis_query)
                
            return analysis_query.id

    def insert_summary(self):
        return None
        
    def update_segment_status(self, analysis_id, segment, status):
        with self.session_scope() as session:
            analysis_query = session.query(
                Analysis
            ).filter(
                Analysis.id == analysis_id
            ).one_or_none()
            if analysis_query:
                self.logger.info('Update analysis status.')
                if segment == 1:
                    analysis_query.segment_1_status = status
                elif segment == 2:
                    analysis_query.segment_2_status = status
                elif segment == 3:
                    analysis_query.segment_3_status = status


    def query_recent(self, nums=5):
        """Get CDR3 records from recent sequencing runs.

        Parameters
        ----------
        nums : int
            The number of recent sequencing runs to retrieve.
        
        Returns
        -------
        query_results : pandas.DataFrame
            Dataframe of CDR3 records in recent sequencing runs.
        """
        session = self.Session()
        # Get recent `nums` runs.
        self.logger.info('Retrieving records of the past {} runs from the database.'.format(nums))
        recent_dates = session.query(DaedalusDB.sequencing_date).distinct().order_by(DaedalusDB.sequencing_date.desc())[:nums]
        recent_dates = [x.sequencing_date for x in recent_dates]

        query_results = pd.read_sql(
            session.query(
                DaedalusDB
            ).filter(
                DaedalusDB.sequencing_date.in_(recent_dates)
            ).statement,
            session.bind
        )
        query_results.drop('id', axis=1, inplace=True)
        return query_results

    def delete_old(self, keep_n=10):
        """Delete old sequencing runs in the database.

        Parameters
        ----------
        keep_n : int
            Number of recent sequencing runs to keep in the database.
        """
        with self.session_scope() as session:
            # Keep specific number of recent sequencing runs, delete all the old ones.
            # Records in Sample and SampleCdr3Xref will be cascade deleted.
            self.logger.info('Deleting records from old sequencing runs. Keep only past {} runs.'.format(keep_n))
            recent_dates = session.query(DaedalusDB.sequencing_date).distinct().order_by(DaedalusDB.sequencing_date.desc())[:keep_n]
            recent_dates = [x.sequencing_date for x in recent_dates]
            session.query(DaedalusDB).filter(~DaedalusDB.sequencing_date.in_(recent_dates)).delete(synchronize_session=False)

    @contextmanager
    def session_scope(self):
        """Create session to execute transaction.
        """
        session = self.Session()
        try:
            yield session
            session.commit()
            self.logger.info('Change to the database committed.')
        except:
            session.rollback()
            self.logger.info('Rollback change to the database.')
            raise
        finally:
            session.close()

    def create_db(self):
        """Create the tables in the database.
        """
        Base.metadata.create_all(self.engine)
        self.logger.info('Database created.')

@event.listens_for(Engine, "connect")
def set_sqlite_pragma(dbapi_connection, connection_record):
    """Enable foreign key constraints when using sqlite3 as backend. 
    """
    if isinstance(dbapi_connection, SQLite3Connection):
        cursor = dbapi_connection.cursor()
        cursor.execute("PRAGMA foreign_keys=ON")
        cursor.close()
