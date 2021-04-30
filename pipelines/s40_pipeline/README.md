# Symbiont-Screener

Filter contamination from filial tgs long reads parental strobemers sets and unsupervised clustering.

## USAGE

### full pipeline 

type ``` ./Symbiont-Screener.sh -h ``` and your get below :

```
Usage   :
    ./Symbiont-Screener.sh [options]

Options :

   Basic parameters:

        --paternal          paternal NGS reads file in FASTQ format.

        --maternal          maternal NGS reads file in FASTQ format.

        --offspring         Offspring sequence file.
                            gzip format file is supported but should end by '.gz'

        --offspring_format  fasta/fastq (default fasta)

        --thread            thread num.
                            [ optional, default 8 threads. ]

        --memory            x (GB) of memory to used by jellyfish.
                            [ optional, default 100GB. ]

  For marker generation:

        -low_depth          predict low_depth            ( default 0 )

        --high_depth        predict high_depth           ( default 0 )
                            program will automatic choose low and high depth threashold when both --low_depth and --high_depth are not setted.
                            if user can predict the depth of host sequences as x , please set like low_depth=x/4 and  high_depth=x*[3 or 5]

  For Trio-formula detection:

        --threshold1        minimum density of parental kmer density (default 0.001)
                            for ONT reads(error rate~=15%), we recommand to use 0.001.
                            for Pacbio reads(error rate<5%), we recommand to use 0.002.

        --threshold2        minimum density of shared kmer density (default 0.01)
                            this default value is very steady.

  For BGMM cluster:

        --cluster           (1/0) use cluster or not. default(0)

        --shortest          length threshold for cluster ( default 5000 )
                            only reads with lengh>shortest can be used for cluster.
                            short reads ( <=5k ) normally create noise points and hamper the cluster result.

        --loop              loop number of BGMM (default 10)

        --python3           PATH to python3 file from anaconda3 ( default python3 )

        --seed              random seed ( default 42 )

        --help              print this usage message.
```

### run this pipeline step by step.

#### first build the markers

first use ```00.BuildTrioLib/build_trio_kmer.sh``` to generate K40 markers.

then convert to randstrobe by ```tools/k40_to_randstrobe_2_10_30.sh```.

#### then get marker density of each long read and filter them

Use 11.GetTrioMatrix_tgs.sh

Here we got trio-only result in trio_only.result.txt and trio_only.filtered_readname.txt

To get the final filtered reads of trio-only, try ```seqkit grep -f trio_only.filtered_readname.txt input.fasta >output.fasta```

This step also genereate trio.4r.matrix for clustering.

#### then cluster to find more reads

[Optional : cluster them if the error rate of long reads is high and the depth of trio-only reads are low.]

##### prepare GC & 3mer matrix 

Use 20.GetGCNmer.sh

This step genereate gc.matrix for clustering.

##### run BGMM

run ```python3 bgm/main_logic.py -t trio.4r.matrix -m gc.matrix   1>cluter.result.txt 2>log ```

Here we got the cluster result in cluter.result.txt
We can get final.result.txt by tools/merge_result.sh
To get the final reads of cluster, try tools/fishSq.sh

## Results

### trio.result.txt

```
column1 :   read-name.
column2 :   trio-base-result    trio-based-formula result 1 means host and 0 means contamination.
column3 :   density of paternal specific markers
column4 :   density of maternal specific markers
column5 :   density of shared  markers
```

all read names filtered by trio-only results are in trio_only.filtered_readname.txt

use ```seqkit grep -f trio_only.filtered_readname.txt input.fasta >output.fasta``` to get host sequences.

### final.result.txt

```
column1 :   read-name.
column2 :   trio-only.          (trio-based-formula result) 1 means host and 0 means contamination.
column3 :   host.               (suggested culster result by hit-count>(loop-time/5)) 1 means host and 0 means contamination)
                                (normally we recommand you to ignore this and choose your own threshold)
column4 :   hit-count           how many number that this read occurs in the best cluster.
column5 :   second-hit-count    how many number that this read occurs in the second-best cluster.
```

Warning: reads shorter than 5000 will be discarded in final.result.txt. However you still can extract short host reads from trio_only.results.txt

example of use fishSq to get final reads

```
./tools/fishSq.sh -f final.result.txt -t 4 -i mixture.fasta -o output.fasta
```



