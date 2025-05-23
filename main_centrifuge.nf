#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl=2

// params
params.dlCentrifugeindex = 'https://genome-idx.s3.amazonaws.com/centrifuge/p_compressed%2Bh%2Bv.tar.gz'
params.runName='Not set'

params.species='S.aureus'
params.taxid=1280

// modules
include {basecall} from './modules/basic.nf'                                     
include {download_centrifuge} from './modules/basic.nf'
include {classify} from './modules/basic.nf'
include {centrifuge_report} from './modules/basic.nf'
include {merge_data} from './modules/basic.nf'
include {results_proc} from './modules/basic.nf'



workflow {
//############################ channels ######################################//
Channel
	.fromPath( "${params.f5s}/*.fast5" )
        .map{ file -> tuple(file.baseName, file) }
	.set{ f5s } 

Channel
	.fromPath( "${params.unblocked}" )
	.set{ unblocked } 

Channel
	.fromPath( "${params.channel_toml}" )
	.set{ channel_toml } 

// ########################### workfow ########################################//
    main:

	basecall(f5s)

	download_centrifuge()

	download_centrifuge.out.db
                    .map{row -> tuple(row[0].simpleName, row[1], row[1])}
                    .set{centdb1}
        
	classify(basecall.out.fastq.combine(centdb1))

	centrifuge_report(classify.out.centFiles.combine(centdb1))

	merge_data(classify.out.centFiles.combine(basecall.out.seqsum, by:0).combine(channel_toml))

	merge_data.out.mixdata.collectFile(name:'mixdata.tab',sort:true, newLine: false, keepHeader:true) 
		.map{ file -> tuple(file.baseName, file) }
	    .set{ mixdata_ch }

	results_proc( mixdata_ch.combine(unblocked))

}