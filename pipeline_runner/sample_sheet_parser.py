from sample_sheet import SampleSheet
import os
import collections
import pandas

#from default_params import *

sample_sheet_path = '/sc1/groups/pls-redbfx/immunoPETE/develop/Daedalus/data/iPETEV2_Expt56_LibPool11_NexSeq_Run7_031320_MANIFEST.csv'

sample_sheet = SampleSheet(sample_sheet_path)

experiment_name = sample_sheet.header.experiment_name
instrument_type = sample_sheet.header.instrument_type

sample_data = []
for sample in sample_sheet.samples:    
    #sampleParams = getDefaultParams()   
    sampData = collections.OrderedDict({
        "experiment_name" : experiment_name,
        "instrument_type" : instrument_type,
        "sample_id" : sample.sample_id,
        "sample_name" : sample.sample_name,
        "i7_index_id" : sample.i7_index_id,
        "i7_index" : sample.index,
        "i5_index_id" : sample.i5_index_id,
        "i5_index" : sample.index2
    })
    sample_data.append(sampData)
    
df = pandas.DataFrame(sample_data)

df.to_csv("Expt56_samples.csv", index=False)





