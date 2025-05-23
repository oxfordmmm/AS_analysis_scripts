// include modules           
include {BATCHREADS} from '../modules/map.nf'
include {GET_ASSEMBLY_META} from '../modules/map.nf'
include {GETREFS} from '../modules/map.nf'
include {DOWNLOADREF} from '../modules/map.nf'
include {SEQTK} from '../modules/map.nf'
include {MINIMAP2} from '../modules/map.nf'
include {SAMTOOLSDEPTH} from '../modules/map.nf'
include {MERGEBAMS} from '../modules/map.nf'
include {MERGEDEPTHS} from '../modules/map.nf'
include {MERGEDEPTHS_TWO} from '../modules/map.nf'
include {DEPTHSTATS} from '../modules/map.nf'
include {PLOTDEPTH} from '../modules/map.nf'



workflow map {
    take:
        fqs
        classifications
    main:
    BATCHREADS(classifications)

    taxids=BATCHREADS.out.taxids
        .transpose(by:2)
        .map(row -> tuple(row[0], row[1], row[2].simpleName, file(row[2])))
        .filter{it[2] != 'fake1'}
        .filter{it[2] != 'fake2'}
        

    uniq_taxids=taxids
        .map(row -> tuple(row[2],row[3]))
        .filter { taxid, taxid_tsv ->
        // loads all lines in the file into memory, then counts
        long count = taxid_tsv.readLines().size()
        count > params.map_reads_min }
        .map(row -> row[0])
        .unique()
        .view()
    
    if (params.taxList != '') {
    //println params.taxList
    uniq_taxids=uniq_taxids
        .filter{it in params.taxList}

    taxids=taxids
        .filter{it[2] in params.taxList}
    }
    
    GET_ASSEMBLY_META()

    GETREFS(uniq_taxids.combine(GET_ASSEMBLY_META.out.viral)
                        .combine(GET_ASSEMBLY_META.out.bacteria))

    DOWNLOADREF(GETREFS.out.ref)

    SEQTK(taxids.combine(fqs, by:[0,1]))

    premap = SEQTK.out.fqs
                .map(row -> tuple(row[2],row[0],row[1], row[3]))
                .combine(DOWNLOADREF.out.ref, by:0)
                .map(row -> tuple(row[1],row[2],row[0], row[3], row[4]))

    MINIMAP2(premap)

    preDepth = MINIMAP2.out.bam_index
                .map(row -> tuple(row[2],row[0],row[1], row[3], row[4]))
                .combine(DOWNLOADREF.out.ref, by:0)
                .map(row -> tuple(row[1],row[2],row[0], row[3], row[4], row[5]))

    SAMTOOLSDEPTH(preDepth)

    DEPTHSTATS(SAMTOOLSDEPTH.out.txt)

    PLOTDEPTH(SAMTOOLSDEPTH.out.txt)

    emit:
    depths = SAMTOOLSDEPTH.out.txt

}

