process BATCHREADS {
    tag {barcode + ' ' + batchID}

    input:
    tuple val(barcode), val(batchID), file('read_tax.csv')

    output:
    tuple val(barcode), val(batchID), path('taxids/*.txt'), emit: taxids, optional: true

    script:
    """
    mkdir taxids
    batchReads.py read_tax.csv
    touch taxids/fake1.txt taxids/fake2.txt
    """
    stub:
    """
    mkdir taxids
    touch taxids/1280.txt taxids/480.txt taxids/567.txt
    """
}

process GET_ASSEMBLY_META{
    label 'online'

    output:
    path('viral_assembly_summary.txt'), emit: viral
    path('bacteria_assembly_summary.txt'), emit: bacteria


    script:
    """
    wget ftp://ftp.ncbi.nih.gov/genomes/refseq/viral/assembly_summary.txt
    mv assembly_summary.txt viral_assembly_summary.txt
    wget ftp://ftp.ncbi.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
    mv assembly_summary.txt bacteria_assembly_summary.txt
    """
    stub:
    """
    touch viral_assembly_summary.txt bacteria_assembly_summary.txt
    """
}

process GETREFS {
    tag {taxid}
    errorStrategy 'retry'
    maxRetries 5

    input:
    tuple val(taxid), path('viral_assembly_summary.txt'), path('bacteria_assembly_summary.txt')

    output:
    tuple val(taxid), path('ref_path.txt'), emit: ref, optional: true

    script:
    """
    getRef.py -t $taxid -b bacteria_assembly_summary.txt -v viral_assembly_summary.txt
    """
    stub:
    """
    touch ref.fa.gz
    """
}

process DOWNLOADREF {
    tag {taxid}
    label 'online'
    errorStrategy 'retry'
    maxRetries 5
    publishDir "refs/${taxid}/" 

    input:
    tuple val(taxid), path('ref_path.txt')

    output:
    tuple val(taxid), path('ref.fa.gz'), emit: ref

    script:
    """
    ref_path="\$(cat ref_path.txt)"
    rsync -avP "\${ref_path}" ref.fa.gz
    """
    stub:
    """
    touch ref.fa.gz
    """
}

process SEQTK {
    tag {barcode + ' ' + batchID + ' ' + taxid}

    input:
    tuple val(barcode), val(batchID), val(taxid), path('reads.txt'), path('in.fastq.gz')

    output:
    tuple val(barcode), val(batchID), val(taxid), path('out.fastq.gz'), emit: fqs

    script:
    """
    seqtk subseq in.fastq.gz reads.txt | gzip > out.fastq.gz
    """
    stub:
    """
    touch out.fastq.gz
    """
}

process MINIMAP2 {
    tag {barcode + ' ' + batchID + ' ' + taxid}

    publishDir "bams/individual/${barcode}/", mode: 'copy'

    input:
    tuple val(barcode), val(batchID), val(taxid), path('in.fastq.gz'), path('ref.fa')

    output:
    tuple val(barcode), val(batchID), val(taxid), path("${batchID}_${taxid}.sorted.bam"), emit: bam
    tuple val(barcode), val(batchID), val(taxid), path("${batchID}_${taxid}.sorted.bam"), path("${batchID}_${taxid}.sorted.bam.bai"),path('ref.fa'),  emit: bam_index


    script:
    """
    minimap2 -t $task.cpus -ax map-ont ref.fa in.fastq.gz |\
        samtools view -F 4 -q 50 -bS - |\
        samtools sort -o ${batchID}_${taxid}.sorted.bam
    samtools index ${batchID}_${taxid}.sorted.bam
    """
    stub:
    """
    touch "${batchID}_${taxid}.sorted.bam"
    """
}

process MERGEBAMS{
    tag {barcode + ' ' + taxid}

    publishDir "${params.outdir}/bams/merged/${barcode}/"

    input:
    tuple val(barcode), val(taxid), path("in_????.sorted.bam")

    output:
    tuple val(barcode), val('final'), val(taxid), path("${barcode}_${taxid}.sorted.bam"), path("${barcode}_${taxid}.sorted.bam.bai"), emit: bams

    script:
    """
    samtools merge ${barcode}_${taxid}.sorted.bam in*.sorted.bam
    samtools index ${barcode}_${taxid}.sorted.bam
    """
    stub:
    """
    ls > files.txt
    touch "${barcode}_${taxid}.sorted.bam"
    """

}

process SAMTOOLSDEPTH {
    tag {barcode + ' ' + batchID + ' ' + taxid}
    errorStrategy 'ignore'
    maxForks 10

    publishDir "depths/${task.process.replaceAll(":","_")}/${barcode}/${taxid}", mode: 'copy'

    input:
    tuple val(barcode), val(batchID), val(taxid), path("${batchID}_${taxid}.sorted.bam"), path("${batchID}_${taxid}.sorted.bam.bai"), path('ref.fa')

    output:
    tuple val(barcode), val(batchID), val(taxid), path("${batchID}_${taxid}.txt"), emit: txt
    tuple val(barcode), val(taxid), path("${batchID}_${taxid}.txt"), emit: merged_txt

    script:
    """
    samtools depth -aa ${batchID}_${taxid}.sorted.bam  --reference ref.fa > ${batchID}_${taxid}.txt
    """
    stub:
    """
    touch "${batchID}_${taxid}.sorted.bam"
    """
}

process MERGEDEPTHS {
    tag {barcode + ' ' + batchID + ' ' + taxid}
    maxForks 1
    errorStrategy 'retry'
    maxRetries 2

    // keep this as symlink, otherwise next process might not have the file in time
    publishDir "depths/${task.process.replaceAll(":","_")}"

    input:
    tuple val(barcode), val(batchID), val(taxid), path("${batchID}_${taxid}.txt")

    output:
    tuple val(barcode), val(taxid), path("${barcode}_${taxid}.txt"), emit: txt

    script:
    bartax=barcode + '_' + taxid
    if (!(bartax in params.depths)) {
        params.depths[bartax] = 1
    }
    else if (params.depths[bartax] == 1) {
        params.depths[bartax] = 2
    }
    if (params.depths[bartax] == 1)
        """
        echo $params.depths > maps.txt
        cp ${batchID}_${taxid}.txt ${barcode}_${taxid}.txt
        """
    else 
        """
        echo $params.depths > maps.txt
        mergeDepths.py -d "${params.outdir}/depths/${task.process.replaceAll(":","_")}/${barcode}_${taxid}.txt" ${batchID}_${taxid}.txt \
            -o ${barcode}_${taxid}.txt
        """
    stub:
    """
    touch "${taxid}.txt"
    """
}

process MERGEDEPTHS_TWO {
    tag {barcode + ' ' + taxid}

    // keep this as symlink, otherwise next process might not have the file in time
    publishDir "${params.outdir}/depths/${task.process.replaceAll(":","_")}"

    input:
    tuple val(barcode), val(taxid), path("depth1.txt"), path("depth2.txt")

    output:
    tuple val(barcode), val(taxid), path("${barcode}_${taxid}.txt"), emit: txt

    script:
    """
    mergeDepths.py -d depth1.txt depth2.txt \
        -o ${barcode}_${taxid}.txt
    """
    stub:
    """
    touch "${barcode}_${taxid}.txt"
    """
}

process DEPTHSTATS {
    tag {barcode + ' ' + taxid}

    publishDir "depth_stats/${task.process.replaceAll(":","_")}/${barcode}/", mode: 'copy'

    input:
    tuple val(barcode), val(batchID), val(taxid), path("${barcode}_${taxid}.txt")

    output:
    tuple val(barcode), val(batchID), val(taxid), path("${barcode}_${batchID}_${taxid}.csv")

    script:
    """
    coverageStats.py -i ${barcode}_${taxid}.txt -o ${barcode}_${batchID}_${taxid}.csv \
        -t ${taxid} -b ${barcode} -a ${batchID}
    """
    stub:
    """
    touch "${taxid}.csv"
    """
}

process PLOTDEPTH {
    tag {barcode + ' ' + taxid}

    publishDir "depth_plots/${task.process.replaceAll(":","_")}/${barcode}/", mode: 'copy'

    input:
    tuple val(barcode), val(batchID), val(taxid), path("${barcode}_${taxid}.txt")

    output:
    tuple val(barcode), val(batchID), val(taxid), path("${barcode}_${batchID}_${taxid}.pdf")

    script:
    """
    plot_depth.py -i ${barcode}_${taxid}.txt -o ${barcode}_${batchID}_${taxid}.pdf \
        -t ${taxid} -b ${barcode} -a ${batchID}
    """
    stub:
    """
    touch "${taxid}.csv"
    """
}