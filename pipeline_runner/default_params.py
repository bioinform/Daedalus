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
        "experiment_name" : None,
        "sequencing_platform" : None,
        "sample_id" : None,
        "sample_name" : None,
        "i7_index_id" : None,
        "i7_index" : None,
        "i5_index_id" : None,
        "i5_index" : None,
        "sample_sheet" : None,
        "project": None,
        "run_folder": None,
        "subsample": 1,
        "subsampleSeed" : 100,
        "inputType" : "DNA",
        "trim_primers" : True,
        "checkSpikein" : False,
        "umi_mode" : True,
        "umi1" : '""',
        "umi2" : "NNNNNNNNN",
        "umi_type" : "R2",
        "vPrimerRef" : os.path.join(pipe_path, "data/Vprimers.fasta"),
        "jPrimerRef" : os.path.join(pipe_path,"data/Jprimers.fasta"),
        "vGeneRef" : os.path.join(pipe_path,"data/Vgenes_cdna_updated.fasta"),
        "jGeneRef" : os.path.join(pipe_path,"data/Jgenes_cdna.fasta"),
        "primerRef" : os.path.join(pipe_path,"data/V2.2_primers.csv"),
        "geneRef" : os.path.join(pipe_path,"data/immunoDB_cdna_updated.csv"),
        "spikeRef" : os.path.join(pipe_path,"data/synthetic_seqs_fixed.fasta"),
        "vGeneScore" : 30,
        "vGeneKmer" : 14,
        "vGeneFrac" : 0.3,
        "jGeneScore" : 20,
        "jGeneKmer" : 12,
        "jGeneFrac" : 0.3,
        "vPrimScore" : 16,
        "vPrimKmer" : 12,
        "vPrimFrac" : 0.6,
        "jPrimScore" : 16,
        "jPrimKmer" : 12,
        "jPrimFrac" : 0.6,
        "spikeScore" : 40,
        "spikeKmer" : 12,
        "spikeFrac" : 0.4,
        "min_read_qual" : 30 ,
        "minGeneIdentity" : 90,
        "max_cdr3_length" : 100,
        "bidding_ratio" : 2,
        "max_steps" : 2,
        "umi_edit_dist" : 1,
        "cdr3_edit_dist" : 1,
        "swifr": os.path.join(pipe_path,"bin/swifr2")
    })
    return ipeteParams
    

