from sqlalchemy import create_engine
from sqlalchemy import Column, Integer, Float, String, Boolean, Date, ForeignKey, UniqueConstraint, Index
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship

Base = declarative_base()

class Analysis(Base):
    """Pipeline Analysis, for tracking individual pipeline runs separately.
    This allows for the a single experiment to be analyzed with multiple parameters. All input parameters are also stored in the database.
        
    """
    __tablename__ = 'analysis'

    ##primary id
    id = Column(Integer, primary_key=True)
    ##seqinfo = relationship('sequencing_run', back_populates='analysis', cascade='save-update, merge, delete')
    seq_id = Column(Integer, ForeignKey('sequencing_run.id', onupdate='CASCADE', ondelete='CASCADE'), nullable=False)
    param_id = Column(Integer, ForeignKey('ipete_params.id', onupdate='CASCADE', ondelete='CASCADE'), nullable=False)

    ##pipeline run info
    analysis_name = Column(String, nullable=False)
    analysis_folder = Column(String, nullable=False)
    analysis_date = Column(Date, nullable=False)
    pipeline_version = Column(String, nullable=False)
    ##what other fields are useful for us to store in the analysis?
    ##Owner: KCL
    ##
    ##pipeline status
    segment_1_status = Column(String, nullable=False)
    segment_2_status = Column(String, nullable=False)
    segment_3_status = Column(String, nullable=False)

    params = relationship('IpeteParams', back_populates='analyses')
    
        
    __table_args__ = (        
        UniqueConstraint('analysis_name', 'pipeline_version', 'seq_id', 'param_id'),
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

    samples = relationship('Sample',
                back_populates='sequencing_run',
                cascade='save-update, merge, delete, delete-orphan')

    __table_args__ = (
        UniqueConstraint('instrument', 'run_number', 'flowcell_id', 'sequencing_date'),
    )

        

class Sample(Base):
    """Table contains sample information and the sequencing run their belongs to.
    """
    __tablename__ = 'sample'

    id = Column(Integer, primary_key=True)
    sample_name = Column(String, nullable=False)
    sample_id = Column(String, nullable=False)
    run_id = Column(Integer, ForeignKey('sequencing_run.id', onupdate='CASCADE', ondelete='CASCADE'), nullable=False)
    lane = Column(Integer)
    i7_index_id = Column(String, nullable=False)
    i7_index = Column(String, nullable=False)
    i5_index_id = Column(String)
    i5_index = Column(String)

    sequencing_run = relationship('SequencingRun', back_populates='samples')

    __table_args__ = (
        UniqueConstraint('sample_name', 'run_id', 'lane'),
        UniqueConstraint('i7_index_id', 'i7_index', 'i5_index_id', 'i5_index')
    )



##get info from Florian, KCL, etc.  
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
        
    ## pipeline parameters
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

    analyses = relationship('Analysis',
        back_populates='params',
        cascade='save-update, merge, delete, delete-orphan')

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


# class SampleSummary(Base):
#     """Table contains pipeline summary stats.
#     """
#     __tablename__ = 'sample_summary'
    
#     id = Column(Integer, primary_key=True)
#     ## which pipeline run
#     analysis_id = Column(Integer, primary_key=True)
    
#     ## pipeline steps
#     pipeline_status = Column(String, nullable=False)

#     ## sample info
#     sample_name = Column(String, nullable=False)
#     sample_id = Column(String, nullable=False)
#     i7_index_id = Column(String, nullable=False)
#     i7_index = Column(String, nullable=False)
#     i5_index_id = Column(String)
#     i5_index = Column(String)
#     ## pipeline steps
#     pipeline_status = Column(String, nullable=False)
#     ## pipeline summary
#     total_reads = Column(Integer, nullable=False)
#     on_target_reads = Column(Integer, nullable=False)
#     percent_on_target = Column(Integer, nullable=False)
#     UMI_families = Column(Integer, nullable=False)
#     UMI_singletons = Column(Integer, nullable=False)
#     UMI_replicates = Column(Integer, nullable=False)
#     UMI_D99 = Column(Integer, nullable=False)
#     UMI_D99_funct = Column(Integer, nullable=False)
#     UMI_D99_non_funct = Column(Integer, nullable=False)
#     percent_functional = Column(Float, nullable=False)
#     cell_count = Column(Integer, nullable=False)
#     TRA = Column(Integer, nullable=False)
#     TRB = Column(Integer, nullable=False)
#     TRG = Column(Integer, nullable=False)
#     TRD = Column(Integer, nullable=False)
#     IGH = Column(Integer, nullable=False)
#     IGK = Column(Integer, nullable=False)
#     IGL = Column(Integer, nullable=False)
#     chimera = Column(Integer, nullable=False)
#     unique_CDR3 = Column(Integer, nullable=False)
#     cdr3_clusters = Column(Integer, nullable=False)
#     D50 = Column(Integer, nullable=False)
#     renyi_0_norm = Column(Float, nullable=False)
#     renyi_1_norm = Column(Float, nullable=False)
#     renyi_2_norm = Column(Float, nullable=False)
#     renyi_inf_norm = Column(Float, nullable=False)
#     simpsons_diversity_index = Column(Float, nullable=False)
#     simpsons_dominance_index = Column(Float, nullable=False)

#     ## report files
#     Vgene_Bam = Column(String, nullable=False)
#     Jgene_Bam = Column(String, nullable=False)
#     Spike_Bam = Column(String, nullable=True)
#     cdr3_report = Column(String, nullable=False)
#     on_target_report = Column(String, nullable=False)
#     dedup_report = Column(String, nullable=False)
#     spike_report = Column(String, nullable=False)
#     diversity_report = Column(String, nullable=False)

#     ##perhaps also include a pdf of plots.

    
    
   
    

##will be too big for sqlite, most likely

# class Cdr3Stats(Base):
#     """Table contains cdr3 information.
#     """
#     __tablename__ = 'pipesummary'
    
#     id = Column(Integer, primary_key=True)

#     cdr3 record
#     cdr3_nt = Column(String, nullable=False)
#     cdr3_aa = Column(String, nullable=False)
#     cdr3_qual = Column(String, nullable=False)
#     cdr3_min_qual = Column(Integer, nullable=False)
#     cdr3_low_qual_count = Column(Integer, nullable=False)
#     v = Column(String, nullable=False)
#     d = Column(String)
#     j = Column(String)
#     umi_seq = Column(String, nullable=False)
#     umi_qual = Column(String, nullable=False)
#     umi_min_qual = Column(Integer, nullable=False)
#     umi_low_qual_count = Column(Integer, nullable=False)
#     family_size = Column(Integer, nullable=False)
#     read_cummulative_fraction = Column(Integer, nullable=False)
# __table_args__ = (
#     Index('cdr3_seq', 'cdr3_nt', 'cdr3_aa', 'umi_seq', 'v', 'j', 'd'),
#     Index('seq_date', 'sequencing_date')
# )
