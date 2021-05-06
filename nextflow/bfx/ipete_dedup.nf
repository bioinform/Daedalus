nextflow.preview.dsl=2

include '../modules/ipete_dedup' params(params)

def ipeteDedup(sample, onTargetUMI, umi_mode) {
    (sample, dedupReport, dedupAll) =  ipete_dedup(sample, onTargetUMI, umi_mode, 
	params.max_cdr3_length, params.max_steps, params.bidding_ratio, params.umi_edit_dist,
	params.cdr3_edit_dist, params.min_read_qual)
    return [sample, dedupReport, dedupAll]
}
