process basecall {
    tag {name}
    label 'gpu'

    publishDir 'fastqs/', mode: 'copy'//, pattern: "*.gz", saveAs: { filename -> "${name}.fastq.gz"}
	
	input:
	tuple val(name), path('file.fast5')

	output:
	tuple val(name), path("${name}/pass/*.fastq.gz"), emit:fastq
	tuple val(name), path("${name}/sequencing_summary.txt"), emit:seqsum


    script:
    """
	guppy_basecaller -i ./ -s ${name} \
        -c dna_r9.4.1_450bps_hac.cfg \
        --device 'auto' \
        --compress_fastq
    """
    stub:
    """
    mkdir -p ${name}/pass
    touch ${name}/sequencing_summary.txt
    touch ${name}/pass/example1.fastq.gz
    """
}


process download_centrifuge {
    output:
    tuple path("centdb/*.1.cf"), path('centdb'), emit: db

    script:
    centdb=params.dlCentrifugeindex
    """
    wget $centdb -O centdb.tar.gz
    tar -xvzf centdb.tar.gz
    mkdir centdb
    mv *.cf centdb
    """

    stub:
    """
    mkdir centdb
    touch centdb/example1.1.cf
    """
}

process classify {
    label 'classifiers'                       
    tag { run }                                                                 
    errorStrategy 'ignore'                                                      
                                                                                
    input:                                                                      
    tuple val(run),file('file.fq.gz'), val(centdbName), file('centdbfol'), file('centdb')
                                                                                
    output:                                                                     
    tuple val(run), file("${run}.cent.tsv.gz"), emit: centFiles          
                                                                                
    script:                                                                     
    """                                                                     
    data=`ls -m *fq.gz | tr -d '\\n' | tr -d ' '`                           
    centrifuge -f -x centdbfol/${centdbName} \
       --mm -q -U \$data \
       -S ${run}.cent.tsv --min-hitlen 16 -k 1                              
    gzip ${run}.cent.tsv                                                    
    """       
    stub:
    """
    touch ${run}.cent.tsv.gz
    """                                                              
}     

process centrifuge_report {
    label 'classifiers'                       
    tag { run }                                                                 
    errorStrategy 'ignore'                     
    publishDir 'centrifuge_reports', mode: 'copy'                                 

    input:
    tuple val(run), file("${run}.cent.tsv.gz"), val(centdbName), file('centdbfol'), file('centdb')

    output:
    tuple val(run), file("${run}_report.txt"), emit: reports

    script:
    """
    zcat ${run}.cent.tsv.gz | \
        centrifuge-kreport -x centdbfol/${centdbName} \
        > ${run}_report.txt
    """
    stub:
    """
    touch ${run}_report.txt
    """
}

process merge_data {
    label 'analysis'                       
    tag { run }                                                                 
    //errorStrategy 'ignore'                     
    publishDir 'mixdata', mode: 'copy'                                 

    input:
    tuple val(run), path("${run}.cent.tsv.gz"), file('seqsum.txt'), path('channels.toml')

	output:
	path("${run}_mixdata.tab"), emit: mixdata

    script:
    """
	merge_data.py --centrifuge ${run}.cent.tsv.gz \
		--seqsum seqsum.txt \
        --toml channels.toml \
		--output ${run}_mixdata.tab
    """
    stub:
    """
    touch ${run}_mixdata.tab
    """
}

process results_proc {
	label 'analysis'
	publishDir 'results', mode: 'copy'

	input:
	tuple val(name), path('mixdata.tab'), path('unblocked.txt')

	output:
	path("*.csv"), emit: csvs
	path("*.pdf"), emit: pdfs

	script:
    runName=params.runName
    species=params.species
    taxid = params.taxid
	"""
	results.py mixdata.tab unblocked.txt $runName $species $taxid
	"""
    stub:
    """
    touch read_length_boxplot.pdf 
    touch read_lengths.csv 
    touch Yields.csv 
    touch yield_time.pdf
    touch Group_counts.csv
    """
}

process KRAKEN2SERVER {
    tag {barcode + ' ' + batchID}
    publishDir "classifications/${barcode}/", mode: 'copy'

    label 'kraken2server'

    input:
    tuple val(barcode), val(batchID), file('in.fastq.gz'), val(krak_name), path('kraken_db')

    output:
    tuple val(barcode), val(batchID), path('read_tax.txt'), emit: read_tax
    tuple val(barcode), val(batchID), path("${barcode}_${batchID}_report.txt"), emit: reports
    tuple val(barcode), val(batchID), path("${barcode}_${batchID}.txt"), emit: txt

    script:
    port=params.kraken2server_port
    host=params.kraken2server_host
    """
    kraken2_client -s in.fastq.gz \
        --port $port --host-ip $host \
        -r ${barcode}_${batchID}_report.txt \
        > ${barcode}_${batchID}.txt

    echo "readID taxID" > read_tax.txt
    awk '{print \$2,\$3}' ${barcode}_${batchID}.txt >> read_tax.txt
    """
    stub:
    """
    touch kraken.txt.gz
    """

}

process KREPORT_BASES {
    tag {barcode + ' ' + batchID}

    publishDir "classifications_bases/${barcode}/", mode: 'copy'

    input:
    tuple val(barcode), val(batchID), file('cent.txt'), path('ktaxonomy.txt')

    output:
    tuple val(barcode), val(batchID), file("${barcode}_${batchID}_bases.txt"), emit: reports

    script:
    """
    create_kreport.py -i cent.txt -t ktaxonomy.txt \
        -o ${barcode}_${batchID}_bases.txt -c kraken2server
    """

    stub:
    """
    touch ${barcode}_${batchID}_bases.txt
    """
}

process GET_GROUP_READS {
    tag {name}
    label 'analysis'
    publishDir "group_reads/${name}/", mode: 'copy'

    input:
    tuple val(name), path('file.fastq.gz'), path('unblocked.txt'), path('channels.toml')

    output:
    tuple val(name), val('AS_unblocked'), path('AS_unblocked.txt'), emit: AS_unblocked
    tuple val(name), val('AS_sequenced'), path('AS_sequenced.txt'), emit: AS_sequenced
    tuple val(name), val('control_sequenced'), path('control_sequenced.txt'), emit: control_sequenced

    script:
    """
    get_groups -f file.fastq.gz -u unblocked.txt -c channels.toml
    """
    stub:
    """
    touch AS_unblocked.txt AS_sequenced.txt control_sequenced.txt
    """
}

process SPLIT_READS {
    tag {name + ' ' + group}
    label 'analysis'
    publishDir 'split_reads', saveAs: { filename -> "${name}_${group}.fastq.gz" }, mode: 'copy'

    input:
    tuple val(name),  path('file.fastq.gz'), val(group), path('reads.txt')

    output:
    tuple val(name), val(group), path('out.fastq.gz'), emit: fastq

    script:
    """
    seqtk subseq file.fastq.gz reads.txt | gzip > out.fastq.gz
    """
    stub:
    """
    touch file.fastq
    """
}