from sqlalchemy import Column, Integer, Float, String, Boolean, Date, ForeignKey, UniqueConstraint, Index
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship

Base = declarative_base()

class Analysis(Base):
    """Pipeline Analysis, for tracking individual pipeline runs separately.

    This allows for the a single experiment to be analyzed with multiple parameters.
    All input parameters are also stored in the database.
    """
    __tablename__ = 'analysis'

    id = Column(Integer, primary_key=True)
    analysis_name = Column(String, nullable=False)
    analysis_folder = Column(String, nullable=False)
    analysis_date = Column(Date, nullable=False)
    pipeline_version = Column(String, nullable=False)
    segment_1_status = Column(String, nullable=False)
    segment_2_status = Column(String, nullable=False)
    segment_3_status = Column(String, nullable=False)

    sequencing_run_samples = relationship(
        'SampleAnalysisXref',
        back_populates='analysis',
        cascade='save-update, merge, delete, delete-orphan'
    )

    __table_args__ = (
        UniqueConstraint('analysis_name', 'analysis_folder', 'analysis_date', 'pipeline_version'),
        Index('pipe_info', 'analysis_name', 'pipeline_version')
    )


class SequencingRun(Base):
    """Table contains sequencing run information.
    """
    __tablename__ = 'sequencing_run'

    id = Column(Integer, primary_key=True)
    sequencing_platform = Column(String, nullable=False)
    instrument = Column(String, nullable=False)
    run_number = Column(Integer, nullable=False)
    flowcell_id = Column(String, nullable=False)
    sequencing_date = Column(Date, nullable=False)
    run_folder = Column(String, nullable=False)
    sample_sheet = Column(String, nullable=False)

    samples = relationship(
        'SequencingRunSampleXref',
        back_populates='sequencing_run',
        cascade='save-update, merge, delete, delete-orphan'
    )

    __table_args__ = (
        UniqueConstraint('instrument', 'run_number', 'flowcell_id', 'sequencing_date'),
    )


class Sample(Base):
    """Table contains sample information and the sequencing run their belongs to.
    """
    __tablename__ = 'sample'

    id = Column(Integer, primary_key=True)
    sample_name = Column(String, nullable=False)
    sample_identifier = Column(String, nullable=False)

    sequencing_runs = relationship(
        'SequencingRunSampleXref',
        back_populates='sample',
        cascade='save-update, merge, delete, delete-orphan'
    )

    __table_args__ = (
        UniqueConstraint('sample_name', 'sample_identifier'),
    )


class SequencingRunSampleXref(Base):
    """Table matchs samples to sequencing runs.
    """
    __tablename__ = 'sequencing_run_sample_xref'

    sample_id = Column(Integer, ForeignKey('sample.id', onupdate='CASCADE', ondelete='CASCADE'), primary_key=True)
    sequencing_run_id = Column(
        Integer, ForeignKey('sequencing_run.id', onupdate='CASCADE', ondelete='CASCADE'), primary_key=True)
    lane = Column(Integer)
    i7_index_id = Column(String, nullable=False)
    i7_index = Column(String, nullable=False)
    i5_index_id = Column(String)
    i5_index = Column(String)

    sample = relationship('Sample', back_populates='sequencing_runs')
    sequencing_run = relationship('SequencingRun', back_populates='samples')

    # analyses = relationship(
    #     'SampleAnalysisXref',
    #     back_populates='sequencing_run_sample_xref',
    #     cascade='save-update, merge, delete, delete-orphan'
    # )


class SampleAnalysisXref(Base):
    """Table matchs samples to analysis.
    """
    __tablename__ = 'sample_analysis_xref'

    # sample_id = Column(
    #     Integer, ForeignKey('sequencing_run_sample_xref.sample_id', onupdate='CASCADE', ondelete='CASCADE'), primary_key=True)
    # sequencing_run_id = Column(
    #     Integer, ForeignKey('sequencing_run_sample_xref.sequencing_run_id', onupdate='CASCADE', ondelete='CASCADE'), primary_key=True)

    sample_id = Column(
        Integer, ForeignKey('sample.id', onupdate='CASCADE', ondelete='CASCADE'), primary_key=True)
    sequencing_run_id = Column(
        Integer, ForeignKey('sequencing_run.id', onupdate='CASCADE', ondelete='CASCADE'), primary_key=True)

    analysis_id = Column(Integer, ForeignKey('analysis.id', onupdate='CASCADE', ondelete='CASCADE'), primary_key=True)
    params_id = Column(Integer, ForeignKey('ipete_params.id', onupdate='CASCADE', ondelete='CASCADE'), primary_key=True)

    sample_summary_id = Column(Integer, ForeignKey('sample_summary.id', onupdate='CASCADE', ondelete='CASCADE'), primary_key=True)

    analysis = relationship('Analysis', back_populates='sequencing_run_samples')
    # sequencing_run_sample = relationship('SequencingRunSampleXref', back_populates='analyses')
    params = relationship('IpeteParams', back_populates='sample_analyses')
    sample_summaries = relationship('SampleSummary', back_populates='sample_analyses')


# #get info from Florian, KCL, etc.

# class SampleMetadata(Base):
#     """Information at the pool level:experiment name,
#     run folder, sample sheet, platform, flowcell_id, etc..
#     """
#     __tablename__ = 'sampleMeta'

#     id = Column(Integer, primary_key=True)

#     ## sequencing info
#     input_ng = Column(String, nullable=False)


class IpeteParams(Base):
    """Table contains CDR3 and umi information.
    """
    __tablename__ = 'ipete_params'

    id = Column(Integer, primary_key=True)

    # pipeline parameters
    subsample = Column(Integer, nullable=False)
    subsampleSeed = Column(Integer, nullable=False)
    primer_targets = Column(String, nullable=False)
    inputType = Column(String, nullable=False)
    trim_primers = Column(Boolean, nullable=False)
    checkSpikein = Column(Boolean, nullable=False)
    umi_mode = Column(Boolean, nullable=False)
    umi1 = Column(String)
    umi2 = Column(String)
    umi_type = Column(String, nullable=False)
    vPrimerRef = Column(String, nullable=False)
    jPrimerRef = Column(String, nullable=False)
    vGeneRef = Column(String, nullable=False)
    jGeneRef = Column(String, nullable=False)
    primerRef = Column(String, nullable=False)
    geneRef = Column(String, nullable=False)
    spikeRef = Column(String, nullable=False)
    vGeneScore = Column(Integer, nullable=False)
    vGeneKmer = Column(Integer, nullable=False)
    vGeneFrac = Column(Float, nullable=False)
    jGeneScore = Column(Integer, nullable=False)
    jGeneKmer = Column(Integer, nullable=False)
    jGeneFrac = Column(Float, nullable=False)
    vPrimScore = Column(Integer, nullable=False)
    vPrimKmer = Column(Integer, nullable=False)
    vPrimFrac = Column(Float, nullable=False)
    jPrimScore = Column(Integer, nullable=False)
    jPrimKmer = Column(Integer, nullable=False)
    jPrimFrac = Column(Float, nullable=False)
    spikeScore = Column(Integer, nullable=False)
    spikeKmer = Column(Integer, nullable=False)
    spikeFrac = Column(Float, nullable=False)
    min_read_qual = Column(Integer, nullable=False)
    minGeneIdentity = Column(Integer, nullable=False)
    max_cdr3_length = Column(Integer, nullable=False)
    bidding_ratio = Column(Float, nullable=False)
    max_steps = Column(Integer, nullable=False)
    umi_edit_dist = Column(Integer, nullable=False)
    cdr3_edit_dist = Column(Integer, nullable=False)

    sample_analyses = relationship(
        'SampleAnalysisXref',
        back_populates='params',
        cascade='save-update, merge, delete-orphan'
    )

    __table_args__ = (
        UniqueConstraint(
            'subsample',
            'subsampleSeed',
            'primer_targets',
            'inputType',
            'trim_primers',
            'checkSpikein',
            'umi_mode',
            'umi1',
            'umi2',
            'umi_type',
            'vPrimerRef',
            'jPrimerRef',
            'vGeneRef',
            'jGeneRef',
            'primerRef',
            'geneRef',
            'spikeRef',
            'vGeneScore',
            'vGeneKmer',
            'vGeneFrac',
            'jGeneScore',
            'jGeneKmer',
            'jGeneFrac',
            'vPrimScore',
            'vPrimKmer',
            'vPrimFrac',
            'jPrimScore',
            'jPrimKmer',
            'jPrimFrac',
            'spikeScore',
            'spikeKmer',
            'spikeFrac',
            'min_read_qual',
            'minGeneIdentity',
            'max_cdr3_length',
            'bidding_ratio',
            'max_steps',
            'umi_edit_dist',
            'cdr3_edit_dist'
        ),
    )

    

class SampleSummary(Base):
    """Table contains pipeline summary stats.
    """
    __tablename__ = 'sample_summary'

    id = Column(Integer, primary_key=True)
    
    ## pipeline summary columns
    status = Column(String, nullable = False)
    total_reads = Column(Integer, nullable=True)
    on_target_reads = Column(Integer, nullable=True)
    percent_on_target = Column(Integer, nullable=True)
    UMI_families = Column(Integer, nullable=True)
    UMI_singletons = Column(Integer, nullable=True)
    UMI_replicates = Column(Integer, nullable=True)
    UMI_D99 = Column(Integer, nullable=True)
    UMI_D99_funct = Column(Integer, nullable=True)
    UMI_D99_non_funct = Column(Integer, nullable=True)
    percent_functional = Column(Float, nullable=True)
    cell_count = Column(Integer, nullable=True)
    TRA = Column(Integer, nullable=True)
    TRB = Column(Integer, nullable=True)
    TRG = Column(Integer, nullable=True)
    TRD = Column(Integer, nullable=True)
    IGH = Column(Integer, nullable=True)
    IGK = Column(Integer, nullable=True)
    IGL = Column(Integer, nullable=True)
    chimera = Column(Integer, nullable=True)

    IGH_percent = Column(Float, nullable=True)
    TRB_percent = Column(Float, nullable=True)
    TRD_percent = Column(Float, nullable=True)
    
    unique_CDR3 = Column(Integer, nullable=True)
    cdr3_clusters = Column(Integer, nullable=True)
    D50 = Column(Float, nullable=True)
    gini_index = Column(Float, nullable=True)
    
    renyi_1_norm = Column(Float, nullable=True)
    renyi_2_norm = Column(Float, nullable=True)
    simpsons_diversity_index = Column(Float, nullable=True)
    simpsons_dominance_index = Column(Float, nullable=True)
    
    TRB_renyi_1_norm = Column(Float, nullable=True)
    TRB_renyi_2_norm = Column(Float, nullable=True)
    TRB_simpsons_diversity_norm = Column(Float, nullable=True)
    TRB_simpsons_dominance_norm = Column(Float, nullable=True)
    TRB_D50 = Column(Float, nullable=True)
    TRB_gini_index = Column(Float, nullable=True)
    TRB_cdr3_clusters = Column(Integer, nullable=True)

    IGH_renyi_1_norm = Column(Float, nullable=True)
    IGH_renyi_2_norm = Column(Float, nullable=True)
    IGH_simpsons_diversity_norm = Column(Float, nullable=True)
    IGH_simpsons_dominance_norm = Column(Float, nullable=True)
    IGH_D50 = Column(Float, nullable=True)
    IGH_gini_index = Column(Float, nullable=True)
    IGH_cdr3_clusters = Column(Integer, nullable=True)

    TRD_renyi_1_norm = Column(Float, nullable=True)
    TRD_renyi_2_norm = Column(Float, nullable=True)
    TRD_simpsons_diversity_norm = Column(Float, nullable=True)
    TRD_simpsons_dominance_norm = Column(Float, nullable=True)
    TRD_D50 = Column(Float, nullable=True)
    TRD_gini_index = Column(Float, nullable=True)
    TRD_cdr3_clusters = Column(Integer, nullable=True)
    
    ## paths to report files
    cdr3_report = Column(String, nullable=True)
    diversity_report = Column(String, nullable=True)
    dedup_report = Column(String, nullable=True)
    on_target_read_report = Column(String, nullable=True)
    
    ##perhaps also include a pdf of plots.
    sample_analyses = relationship(
        'SampleAnalysisXref',
        back_populates='sample_summaries',
        cascade='save-update, merge, delete-orphan'
    )
