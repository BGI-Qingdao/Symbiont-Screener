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
	--low_depth         predict low_depth            ( default 0 )
        --high_depth        predict high_depth           ( default 0 )  
	                    program will automatic choose low and high depth threashold when --low_depth and --high_depth are not setted.
			    if user can predict the depth of host sequences as x , please set like low_depth=x/4 and  high_depth=x*[3 or 5]
  For Trio-formula detection:
  
        --threshold1        minimum density of parental kmer density (default 0.001)
	                    for Pacbio reads, we recommand to use 0.002.
        --threshold2        minimum density of shared kmer density (default 0.01)
	                    this default value is very steady.
  
  For BGMM cluster:
  
        --cluster           (1/0) use cluster or not. default(0)
	
	--shortest          length threshold for cluster ( default 5000 )
	                    only reads with lengh>shortest can be used for cluster.
			    short reads ( 1k-5k ) normally create noise points and hamper the cluster result.
        
	--loop              loop number of BGMM (default 10) 
        --python3           PATH to python3 file from anaconda3 ( default python3 )
      
        --seed              random seed ( default 42 )
        --help              print this usage message.
```

### step by step

* 00.BuildTrioLib.sh
* 11.GetTrioMatrix_tgs.sh
    *  here we got trio-only result in trio.result.txt 
* 20.GetGCNmer.sh
* bgm/main_logic.py <your paramters>    1>cluter.result.txt 2>log 
    *  here we got the cluster result in cluter.result.txt
    *  we can get final.result.txt by tools/merge_result.sh

## Results 

## trio.result.txt
### fish final target fastq
## final.result.txt

```
column1 :   read-name/barcode-name.
column2 :   priori.             (prior-formula result) 1 means host and 0 means contamination.
column3 :   host.               (suggested culster result by hit-count>5) 1 means host and 0 means contamination; normally we recommand you to ignore this.
column4 :   hit-count           how many number that this read occurs in the best cluster.
column5 :   second-hit-count    how many number that this read occurs in the second-best cluster.
```

### fish final target fastq

use tools/fishSeq.sh
