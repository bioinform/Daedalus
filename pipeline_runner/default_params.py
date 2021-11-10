import collections
import os

def setParamValues(params, updateValues):
    for k,v in updateValues.items():
        if k in params:
            params[k] = v
        else: ##create new values
            params[k] = v
    return params


def get_pipeline_path():
    pipe_path = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
    return pipe_path

def getDefaultParams():
    pipe_path = get_pipeline_path()
    ipeteParams = collections.OrderedDict({
        "experiment_name" : None, ##the analysis name for the run, defined by the assay team... something meaninful
        "sequencing_platform" : None, ##the illumina platform, parsed from the sample_sheet
        "sample_id" : None, ## the sample id, parsed from the first column of sample_sheet
        "sample_name" : None, ## the name defined in the sample sheet
        "i7_index_id" : None, ## which index1 was used
        "i7_index" : None, ## the index1 sequence
        "i5_index_id" : None, ## which index2 was used
        "i5_index" : None, ## the index2 sequence
        "sample_sheet" : None, ##the path to the sample_sheet
        "project": None, ##the name used to generate the pipeline output parent folder
        "run_folder": None, ##the path to the illumina run folder
        "subsample": 1, ##subsample parameter, 1 for all reads, >1 specifies read amount
        "subsampleSeed" : 100, ##random number seed, used by seqtk
        "inputType" : "DNA", ##the input type, specified by manifest generator, provided as a log for DNA vs RNA
        "trim_primers" : True, ##whether or not to run the primer trimming steps in the nextflow pipeline
        "checkSpikein" : False, ## whether or not to run the spikein steps in the nextflow pipeline
        "umi_mode" : True, ##whether or not to run ipete_dedup in the nextflow pipeline
        "umi1" : '""', ##if there is a UMI on read1, specify the pattern
        "umi2" : "NNNNNNNNN", ##if there is a UMI on read2, specify the pattern
        "umi_type" : "R2", ##either R1, R2, or both. Specify if a UMI starts read1, read2, or both
        "vPrimerRef" : os.path.join(pipe_path, "data/Vprimers.fasta"), ##path to vgene primers
        "jPrimerRef" : os.path.join(pipe_path,"data/Jprimers.fasta"), ##path to jgene primers
        "vGeneRef" : os.path.join(pipe_path,"data/Vgenes_cdna_updated.fasta"), #path to vgene reference
        "jGeneRef" : os.path.join(pipe_path,"data/Jgenes_cdna.fasta"), #path to jgene reference
        "primerRef" : os.path.join(pipe_path,"data/V2.2_primers.csv"), #path to primer reference csv file
        "geneRef" : os.path.join(pipe_path,"data/immunoDB_cdna_updated.csv"), #path to gene reference csv file
        "spikeRef" : os.path.join(pipe_path,"data/synthetic_seqs_fixed.fasta"), #path to spikein reference
        "vGeneScore" : 30, ##Vgene score, minimum alignment score to keep a Vgene alignment
        "vGeneKmer" : 14, ##Vgene kmer size, the kmer length used to index Vgene reference seqs and search reads against
        "vGeneFrac" : 0.3, ##Vgene fraction, any reference with kmer coverage greater than this fraction is aligned against
        "jGeneScore" : 20, ##jgene score, minimum alignment score to keep a Jgene alignment
        "jGeneKmer" : 12, ##Jgene kmer size, the kmer length used to index Jgene reference seqs and search reads against
        "jGeneFrac" : 0.3, ##Jgene fraction, any reference with kmer coverage greater than this fraction is aligned against
        "vPrimScore" : 16, ##V primer score, minimum alignment score to keep a V primer alignment
        "vPrimKmer" : 12, ##V primer kmer size, the kmer length used to index V primer reference seqs and search reads against
        "vPrimFrac" : 0.6,##V primer fraction, any reference with kmer coverage greater than this fraction is aligned against
        "jPrimScore" : 16, ##J primer score, minimum alignment score to keep a J primer alignment
        "jPrimKmer" : 12, ##J primer kmer size, the kmer length used to index J primer reference seqs and search reads against
        "jPrimFrac" : 0.6, ##J primer fraction, any reference with kmer coverage greater than this fraction is aligned against
        "spikeScore" : 40, ##spikein score, minimum alignment score to keep a spikein alignment
        "spikeKmer" : 12, ##spikein kmer size, the kmer length used to index spikein reference seqs and search reads against
        "spikeFrac" : 0.4, ##spikein fraction, any reference with kmer coverage greater than this fraction is aligned against
        "min_read_qual" : 30, ##the minimum average read quality filter, since CDR3 have no reference, Quality scores are useful to avoid artifacts
        "minGeneIdentity" : 90, ##For V and J gene alignments, a read is considered on target, only if the alignments are better than this percentage identity
        "max_cdr3_length" : 100, ##some V and Jgenes are nearby in the reference. This is the maximum allowable CDR3 length to be reported
        "bidding_ratio" : 2, ## For the dedup step, the ratio between two nodes determing when one node is an error from the other
        "max_steps" : 2, ##for each node in the graph, edges will be scored up to this many steps away
        "umi_edit_dist" : 1, ##the levenshtein edit distance used to search for similar UMIs
        "cdr3_edit_dist" : 1, ##the levenshtein edit distance used to search for similar CDR3s
        "swifr": os.path.join(pipe_path,"bin/swifr2"), ## the path to the swifr binary
        "pipeline_summary": os.path.join(pipe_path,"pipeline_runner/pipeline_summary.py") ## the path to the pipeline summary script
    })
    return ipeteParams
    

