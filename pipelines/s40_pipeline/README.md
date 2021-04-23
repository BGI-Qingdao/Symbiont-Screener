# Symbiont-Screener
Filter contamination from filial tgs long reads parental strobemers sets and unsupervised clustering.


## USAGE

### full pipeline 

type ``` ./Symbiont-Screener.sh -h ``` and your get below :

```
Usage   :
    ./Symbiont-Screener.sh [options]

Options :

        --paternal          paternal NGS reads file in FASTQ format.

        --maternal          maternal NGS reads file in FASTQ format.

        --offspring         Offspring sequence file.
                            gzip format file is supported but should end by '.gz'

        --offspring_format  fasta/fastq (default fasta)

        --threshold1        minimum density of parental kmer density (default 0.02)

        --threshold2        minimum density of shared kmer density (default 0.1)

        --loop              loop number of BGMM (default 30)

        --thread            thread num.
                            [ optional, default 8 threads. ]

        --memory            x (GB) of memory to used by jellyfish.
                            [ optional, default 100GB. ]

        --python3           PATH to python3 file from anaconda3 ( default python3 )

        --seed              random seed ( default 42 )

        --help              print this usage message.
```

### step by step

* 00.BuildTrioLib.sh
* 11.GetTrioMatrix_tgs.sh
* 20.GetGCNmer.sh
* bgm/main_logic.py
* merge result by tools/merge_result.sh

## Get reads from result & raw reads

try tools/ExtractReads_tgs.sh


## Format of final.result.txt

```
column1 :   read-name/barcode-name.
column2 :   priori.             (prior-formula result) 1 means host and 0 means contamination.
column3 :   host.               (suggested culster result by hit-count>5) 1 means host and 0 means contamination; normally we recommand you to ignore this.
column4 :   hit-count           how many number that this read occurs in the best cluster.
column5 :   second-hit-count    how many number that this read occurs in the second-best cluster.
```
