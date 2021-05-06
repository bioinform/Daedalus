import logging
import warnings

import pandas as pd

from contextlib import contextmanager
from sqlite3 import Connection as SQLite3Connection

from sqlalchemy import create_engine, event, and_
from sqlalchemy import exc as sa_exc
from sqlalchemy.engine import Engine
from sqlalchemy.orm import sessionmaker

# from daedalus_db.models import Base, Analysis, SequencingRun, Sample, IpeteParams #, SampleSummary
from models import (
    Base,
    Analysis,
    SequencingRun,
    Sample,
    IpeteParams,
    SampleSummary,
    SequencingRunSampleXref,
    SampleAnalysisXref,
)


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

    def __init__(
        self, url='sqlite:///daedalus.db', echo=False, newDb=False, timeout=60
    ):
        self.url = url
        self.engine = create_engine(
            self.url, echo=echo, connect_args={'timeout': timeout}
        )
        self.Session = sessionmaker(bind=self.engine)
        self.logger = logging.getLogger(__name__)
        if newDb == True:
            self.logger.info("creating new database")
            self.create_db()

    def checkRunInfo(self, run_info, overwrite=False):
        with self.session_scope() as session:
            self.logger.info('Inserting records to the database.')
            instrument, run_number, flowcell_id, sequencing_date = (
                run_info["instrument"][0],
                run_info["run_number"][0],
                run_info["flowcell_id"][0],
                run_info["sequencing_date"][0],
            )

            existing_run = (
                session.query(DaedalusDB)
                .filter(
                    and_(
                        DaedalusDB.instrument == instrument,
                        DaedalusDB.run_number == run_number,
                        DaedalusDB.flowcell_id == flowcell_id,
                        DaedalusDB.sequencing_date == sequencing_date,
                    )
                )
                .first()
            )
            if not existing_run or overwrite:
                session.query(DaedalusDB).filter(
                    and_(
                        DaedalusDB.instrument == instrument,
                        DaedalusDB.run_number == run_number,
                        DaedalusDB.flowcell_id == flowcell_id,
                        DaedalusDB.sequencing_date == sequencing_date,
                    )
                ).delete()
                session.bulk_insert_mappings(DaedalusDB, run_info.to_dict('records'))

    def insert_sequencing_run(self, session, sequencing_info, overwrite):
        """Insert Sequencing Run info to the database
            
            Parameters
            ----------
            session : SqlAlchemy Session object
            sequencing_info : pandas.DataFrame
            overwrite : bool
                Whether to overwrite the records in the database if the sequencing run is already stored in the database.
                This will delete the records of existing run and insert the records provided.

            Returns : sequencing_rrun_query
        """
        assert sequencing_info.shape == (1, 7)
        (
            sequencing_platform,
            instrument,
            run_number,
            flowcell_id,
            sequencing_date,
            sample_sheet,
            run_folder,
        ) = sequencing_info.iloc[0]
        sequencing_run_query = (
            session.query(SequencingRun)
            .filter(
                and_(
                    SequencingRun.instrument == instrument,
                    SequencingRun.run_number == run_number,
                    SequencingRun.flowcell_id == flowcell_id,
                    SequencingRun.sequencing_date == sequencing_date,
                )
            )
            .one_or_none()
        )
        if not sequencing_run_query or overwrite:
            session.query(SequencingRun).filter(
                and_(
                    SequencingRun.instrument == instrument,
                    SequencingRun.run_number == run_number,
                    SequencingRun.flowcell_id == flowcell_id,
                    SequencingRun.sequencing_date == sequencing_date,
                )
            ).delete()
            self.logger.info('Overwrite sequencing run.')
            sequencing_run_query = SequencingRun(**sequencing_info.iloc[0].to_dict())
            session.add(sequencing_run_query)
        return sequencing_run_query

    def insert_analysis_info(self, session, analysis_info, overwrite):
        """Insert SampleInfo to the database
            
            Parameters
            ----------
            session : SqlAlchemy Session object
            analysis_info : pandas.DataFrame
            overwrite : bool
                Whether to overwrite the records in the database if the sequencing run is already stored in the database.
                This will delete the records of existing run and insert the records provided.

            Returns : a dictionary of SqlAlchemy sample_query for each sample.
        """
        a = analysis_info.iloc[0].to_dict()
        analysis_query = (
            session.query(Analysis)
            .filter(
                and_(
                    Analysis.analysis_name == a['analysis_name'],
                    Analysis.pipeline_version == a['pipeline_version'],
                )
            )
            .one_or_none()
        )
        if not analysis_query or overwrite:
            if analysis_query:
                self.logger.info('Overwrite analysis.')
                session.delete(analysis_query)
            analysis_query = Analysis(**a)
            self.logger.info('Insert analysis.')
            session.add(analysis_query)
        return analysis_query

    def insert_ipete_param(self, session, param_info, overwrite):
        """Insert IPete params to the database
            
            Parameters
            ----------
            session : SqlAlchemy Session object
            param_info : pandas.DataFrame
            overwrite : bool
                Whether to overwrite the records in the database if the sequencing run is already stored in the database.
                This will delete the records of existing run and insert the records provided.

            Returns : a dictionary of param_query for each sample.
        """
        param_queries = {}
        for p in param_info.to_dict('records'):
            param_query = (
                session.query(IpeteParams)
                .filter(
                    and_(
                        IpeteParams.subsample == p['subsample'],
                        IpeteParams.subsampleSeed == p['subsampleSeed'],
                        IpeteParams.primer_targets == p['primer_targets'],
                        IpeteParams.inputType == p['inputType'],
                        IpeteParams.trim_primers == p['trim_primers'],
                        IpeteParams.checkSpikein == p['checkSpikein'],
                        IpeteParams.umi_mode == p['umi_mode'],
                        IpeteParams.umi1
                        == (p['umi1'] if not pd.isna(p['umi1']) else None),
                        IpeteParams.umi2
                        == (p['umi2'] if not pd.isna(p['umi2']) else None),
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
                        IpeteParams.cdr3_edit_dist == p['cdr3_edit_dist'],
                    )
                )
                .one_or_none()
            )
            if not param_query or overwrite:
                if param_query:
                    self.logger.info('Overwrite params.')
                    # bad logic here.
                    # session.delete(param_queries[p['sample_name']])
                else:
                    # Make a dictionary from params data that contains only columns in the table
                    p_ = {
                        k: p[k] for k in p if k in IpeteParams.__table__.columns.keys()
                    }
                    param_query = IpeteParams(**p_)
                    session.add(param_query)

            param_queries[p['sample_name']] = param_query

        return param_queries

    def insert_sample_summary(self, session, pipeline_results, overwrite):
        """Insert IPete params to the database
            
            Parameters
            ----------
            session : SqlAlchemy Session object
            pipeline_results : pandas.DataFrame
            overwrite : bool
                Whether to overwrite the records in the database if the sequencing run is already stored in the database.
                This will delete the records of existing run and insert the records provided.

            Returns : a dictionary of param_query for each sample.
        """
        sample_summary_queries = {}
        for ssum in pipeline_results.to_dict('records'):
            sample_summary_query = (
                session.query(SampleSummary)
                .filter(
                    and_(
                        SampleSummary.cdr3_report == ssum['cdr3_report'],
                        SampleSummary.diversity_report == ssum['diversity_report'],
                        SampleSummary.dedup_report == ssum['dedup_report'],
                        SampleSummary.on_target_read_report == ssum['on_target_read_report']
                    )
                )
                .one_or_none()
            )
            if not sample_summary_query or overwrite:
                if sample_summary_query:
                    self.logger.info('Overwite sample summary')
                    session.delete(sample_summary_query)
                # Make a dictionary from sample_summary data that contains only columns in the table
                ssum_ = {
                    k: ssum[k] for k in ssum if k in SampleSummary.__table__.columns.keys()
                }
                sample_summary_query = SampleSummary(**ssum_)
                session.add(sample_summary_query)
            sample_summary_queries[ssum['sample_name']] = sample_summary_query

        return sample_summary_queries

    def insert_samples(
        self,
        session,
        sample_info,
        analysis_query,
        sequencing_run_query,
        ipete_param_queries,
        sample_summary_queries,
        overwrite,
    ):
        """Insert SampleInfo to the database
            
            Parameters
            ----------
            session : SqlAlchemy Session object
            sample_info : pandas.DataFrame
            sequencing_run_query: SqlAlchemy query for sequencing run
            overwrite : bool
                Whether to overwrite the records in the database if the sequencing run is already stored in the database.
                This will delete the records of existing run and insert the records provided.

            Returns : a dict of SqlAlchemy seqrun_sample_query for each sample (key=sample_name).
        """
        for s in sample_info.to_dict('records'):
            # Inserting the sample, if it's new
            sample_query = (
                session.query(Sample)
                .filter(
                    and_(
                        Sample.sample_name == s['sample_name'],
                        Sample.sample_identifier == s['sample_id'],
                    )
                )
                .one_or_none()
            )
            if not sample_query:
                sample_query = Sample(
                    sample_name=s['sample_name'], sample_identifier=s['sample_id']
                )
                session.add(sample_query)

            # Inserting the seq_run_sample_xref
            seqrun_sample_xref_query = (
                session.query(SequencingRunSampleXref)
                .filter(
                    and_(
                        SequencingRunSampleXref.sample_id == sample_query.id,
                        SequencingRunSampleXref.sequencing_run_id
                        == sequencing_run_query.id,
                    )
                )
                .one_or_none()
            )
            if not seqrun_sample_xref_query or overwrite:
                if seqrun_sample_xref_query:
                    self.logger.info('Overwite seqrun_sample_query_xref')
                    session.delete(seqrun_sample_xref_query)

                seqrun_sample_xref_query = SequencingRunSampleXref(
                    sample_id=sample_query.id,
                    sequencing_run_id=sequencing_run_query.id,
                    lane=(s['lane'] if not pd.isna(s['lane']) else None),
                    i7_index_id=s['i7_index_id'],
                    i7_index=s['i7_index'],
                    i5_index_id=(
                        s['i5_index_id'] if not pd.isna(s['i5_index_id']) else None
                    ),
                    i5_index=(s['i5_index'] if not pd.isna(s['i5_index']) else None),
                )
                session.add(seqrun_sample_xref_query)

            # insert seq_run_sample_xref
            analysis_sample_xref_query = (
                session.query(SampleAnalysisXref)
                .filter(
                    and_(
                        SampleAnalysisXref.sample_id == sample_query.id,
                        SampleAnalysisXref.sequencing_run_id == sequencing_run_query.id,
                        SampleAnalysisXref.analysis_id == analysis_query.id,
                        SampleAnalysisXref.params_id
                        == ipete_param_queries[s['sample_name']].id,
                        SampleAnalysisXref.sample_summary_id
                        == sample_summary_queries[s['sample_name']].id,
                    )
                )
                .one_or_none()
            )
            if not analysis_sample_xref_query or overwrite:
                if analysis_sample_xref_query:
                    self.logger.info('Overwrite analysis_sample_xref_query')
                    session.delete(analysis_sample_xref_query)
                analysis_sample_xref_query = SampleAnalysisXref(
                    sample_id=sample_query.id,
                    sequencing_run_id=sequencing_run_query.id,
                    analysis_id=analysis_query.id,
                    params_id=ipete_param_queries[s['sample_name']].id,
                    sample_summary_id=sample_summary_queries[s['sample_name']].id,
                )
                session.add(analysis_sample_xref_query)

    def insert(self, records, overwrite=False):
        """Insert records from a sequencing run to the database.

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
            self.logger.info('Inserting records to the database.')
            # insert sequencing run
            sequencing_run_query = self.insert_sequencing_run(
                session, records['sequencing_info'], overwrite
            )

            # insert analysis info
            analysis_query = self.insert_analysis_info(
                session, records['analysis_info'], overwrite
            )

            # insert sample params (returns a dic of qparam queries, one per sample)
            ipete_param_queries = self.insert_ipete_param(
                session, records['param_info'], overwrite
            )

            # insert sample summaries / pipeline_results
            sample_summary_queries = self.insert_sample_summary(
                session, records['pipeline_results'], overwrite
            )

            # insert all samples : this will insert all the xrefs too!
            #     session,
            # sample_info,
            # analysis_query,
            # sequencing_run_query,
            # ipete_param_queries,
            # sample_summary_queries,
            seqrun_sample_xref_queries = self.insert_samples(
                session,
                records['sample_info'],
                analysis_query,
                sequencing_run_query,
                ipete_param_queries,
                sample_summary_queries,
                overwrite,
            )

    def update_segment_status(self, analysis_id, segment, status):
        with self.session_scope() as session:
            analysis_query = (
                session.query(Analysis).filter(Analysis.id == analysis_id).one_or_none()
            )
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
        self.logger.info(
            'Retrieving records of the past {} runs from the database.'.format(nums)
        )
        recent_dates = (
            session.query(DaedalusDB.sequencing_date)
            .distinct()
            .order_by(DaedalusDB.sequencing_date.desc())[:nums]
        )
        recent_dates = [x.sequencing_date for x in recent_dates]

        query_results = pd.read_sql(
            session.query(DaedalusDB)
            .filter(DaedalusDB.sequencing_date.in_(recent_dates))
            .statement,
            session.bind,
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
            self.logger.info(
                'Deleting records from old sequencing runs. Keep only past {} runs.'.format(
                    keep_n
                )
            )
            recent_dates = (
                session.query(DaedalusDB.sequencing_date)
                .distinct()
                .order_by(DaedalusDB.sequencing_date.desc())[:keep_n]
            )
            recent_dates = [x.sequencing_date for x in recent_dates]
            session.query(DaedalusDB).filter(
                ~DaedalusDB.sequencing_date.in_(recent_dates)
            ).delete(synchronize_session=False)

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
