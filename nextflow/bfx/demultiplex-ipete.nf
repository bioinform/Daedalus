nextflow.preview.dsl=2

include '../modules/demultiplex-ipete' params(params)

def doBcl2Fastq(raw_dir, sample_sheet) {
   (fastqs, reports, stats) = bcl2fastq(raw_dir, sample_sheet)

//     if(params.runBcl2Fastq) {
//     bcl2fastq_out.collate(2)
//         .map {fastq_pair -> [sampleShort(fastq_pair[0]), fastq_pair[0], fastq_pair[1]]}
//         .filter { sample -> !(sample[0] =~ /^Undetermined.*/) }
//         .set { groupedDemux }
// }

    
   return fastqs.flatten()
}
