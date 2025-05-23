#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl=2

// params
params.kraken2server_port = 8080
params.kraken2server_host = '163.1.213.232'
params.krakdb = '/mnt/nanostore/dbs/kraken2_db'
params.ktaxonomy = '/mnt/nanostore/dbs/kraken2_db/ktaxonomy.tsv'
params.map_reads_min = 10
params.taxList=''

// modules
include {GET_GROUP_READS} from './modules/basic.nf'
include {SPLIT_READS} from './modules/basic.nf'
include {KRAKEN2SERVER} from './modules/basic.nf'
include {KREPORT_BASES} from './modules/basic.nf'


// workflow
include {map} from './workflows/map.nf'

workflow {
//############################ channels ######################################//
Channel
	.fromPath( "${params.fqs}/*" )
    .map{ file -> tuple(file.simpleName, file) }
	.set{ fqs } 

Channel
	.fromPath( "${params.unblocked}/*" )
	.map{ file -> tuple(file.baseName, file) }
	.set{ unblocked } 

Channel
	.fromPath( "${params.channel_toml}/*" )
	.map{ file -> tuple(file.baseName, file) }
	.set{ channel_toml } 

Channel.fromPath("${params.krakdb}")
    .map{row -> tuple( row.simpleName, file(row) ) }
    .set{classifier_db}

Channel.fromPath("${params.ktaxonomy}")
	.set{ktaxonomy}

// ########################### workfow ########################################//
    main:

	GET_GROUP_READS(fqs.combine(unblocked, by:0).combine(channel_toml, by:0))

	groups=GET_GROUP_READS.out.AS_unblocked.mix(GET_GROUP_READS.out.AS_sequenced, GET_GROUP_READS.out.control_sequenced)

	SPLIT_READS(fqs.combine(groups, by:0))

	KRAKEN2SERVER(SPLIT_READS.out.fastq.combine(classifier_db))

	KREPORT_BASES(KRAKEN2SERVER.out.txt.combine(ktaxonomy))

	map(SPLIT_READS.out.fastq, KRAKEN2SERVER.out.read_tax)
}