# AS_analysis_scripts
Scripts used to analyse adaptive sampling data. This workflow was written as a one off analysis workflow for our adaptive sampling evaluation. It is here to provide information on what analysis was conducted, as a resource of information. It is not actively being developed nor designed as a standalone workflow for use with other peoples data. You are welcome to take and adapt the code, however. 

## Installing
Clone this repository to a location on your computer. Conda is required and the environment will automatically be built when running with the conda profile.

## Usage

### inputs
This workflow requires the channel tomls and lists of unblocked reads from the gridion

-Fastq files(`--fqs`), concattenated togther, one fastq file per sample/barcode.

-unblocked reads (`--unblocks`), text file of read ids unblocked by the sequencer, provided by the gridion/minion.

-channel groups (`--channel_toml`), a toml file with the groups of channels used for adaptive sampling (AS) or control (non-AS).

-kraken taxomomy (`-ktaxonomy`), a TSV file containing the kraken taxonomy used for the classification.


### Running the workflow
The workflow can be run with the following commands, make sure to change the paths accodingly.

```bash
nextflow run /path/you/downloaded/to/adaptive_sampling_scripts/main.nf \
	--fqs /your/data/fastqs \
	--unblocked unblocks \
	--channel_toml /path/you/downloaded/to/adaptive_sampling_scripts/tomls \
	--ktaxonomy /mnt/nanostore/dbs/kraken2_db/ktaxonomy.tsv \
	-with-trace \
	-resume -profile conda
```

### Outputs

-classifications, kraken2 reports for each group (AS, Non-AS) and sample (barcode).

-classifications_bases, the same, but with bases instead of reads.

-depth_stats, csv file with mapping coverage stats. One file per AS group, sample and speices TaxID. Example below for E.coli, chromasome + 3 plasmids. 

| Barcode   | batch        | taxid | chrom       |    length |       bases |    avDepth | position cov1 | position cov10 | covBreadth1x | covBreadth10x |
| --------- | ------------ | ----- | ----------- | --------- | ----------- | ---------- | ------------- | -------------- | ------------ | ------------- |
| sample_1a | AS_sequenced |   562 | NC_002127.1 |     3,306 |   1,665,828 | 2,478.911… |           672 |            577 |       0.203… |        0.175… |
| sample_1a | AS_sequenced |   562 | NC_002128.1 |    92,721 |   1,501,670 |    59.233… |        25,352 |         23,905 |       0.273… |        0.258… |
| sample_1a | AS_sequenced |   562 | NC_002695.2 | 5,498,578 | 220,371,387 |    56.022… |     3,933,629 |      3,868,375 |       0.715… |        0.704… |

-depths, the output from samtools depth