# ConFilter
Filter contamination from filial tgs long reads or stlfr co-barcode-reads by parental kmer sets and unsupervised clustering.


## USAGE

### full pipeline 

type ``` ./ConFilter.sh -h ``` and your get below :

```
Usage   :
    ./ConFilter.sh [options]

Options :

        --paternal          paternal NGS reads file in FASTQ format.
                            ( note : gzip format is NOT supported. )
                            [ optional, needed when --use_existing_libs off ]

        --maternal          maternal NGS reads file in FASTQ format.
                            ( note : gzip format is NOT supported. )
                            [ optional, needed when --use_existing_libs off ]

        --offspring         Offspring sequence file.
                            gzip format file is supported but should end by '.gz'

        --offspring_format  fasta/fastq (default fasta)

        --threshold1        minimum density of parental kmer density (default 0.001)

        --threshold2        minimum density of shared kmer density (default 0.1)

        --kmer               kmer-size (default 21. ]

        --nmer              nmer for gc_nmer(default 2)

        --sequence_platform tgs/stlfr (default tgs)

        --loop              loop number of BGMM (default 30)

        --thread            thread num.
                            [ optional, default 8 threads. ]

        --memory            x (GB) of memory to used by meryl.
                            [ optional, default 50GB. ]

        --python3           PATH to python3 file from anaconda3 ( default python3 )

        --help              print this usage message.
```

### step by step

* 00.BuildTrioLib.sh
* 11.GetTrioMatrix_tgs.sh or 12.GetTrioMatrix_stlfr.sh
* 20.GetGCNmer.sh
* bgm/main_logic.py
* merge result by tools/merge_result.sh

## Get reads from result & raw reads

try tools/ExtractReads_tgs.sh or tools/ExtractReads_stlfr.sh 


## Format of final.result.txt

```
column1 :   read-name/barcode-name.
column2 :   priori. (prior-formula result) 1 means host and 0 means contamination.
column3 :   host.   (final culster result) 1 means host and 0 means contamination.
column4 :   hit-count   how many number that this read occurs in the best cluster.
```


