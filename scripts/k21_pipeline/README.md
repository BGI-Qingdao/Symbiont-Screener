# Symbiont-Screener (kmer version)

Separate host's sequences from contamination for TGS long reads by trio-binning markers and unsupervised clustering.

## USAGE of full pipeline 

type ``` ./Symbiont-Screener.sh -h ``` and you will get below information:

```
Usage   :
    ./Symbiont-Screener.sh [options]

Options :

   Basic parameters:

        --paternal          contaminated paternal NGS read file in FASTQ format.

        --maternal          contaminated maternal NGS read file in FASTQ format.

        --offspring         contaminated offspring's TGS read file.
                            gzip format file is supported but should end by '.gz'

        --offspring_format  fasta/fastq (default fasta)

        --thread            thread num.
                            [ optional, default 8 threads. ]
                            
        --memory            x (GB) memory used for meryl
                            [ optional, default 50GB. ]

  For marker generation:

        --low_depth         estimated lower depth for k-mer histogram (default 0)

        --high_depth        estimated higher depth for k-mer histogram (default 0)
                            this pipeline will automatically choose lower and higher depth thresholds when both --low_depth and --high_depth are not set.
                            if the user estimates that sequencing coverage or depth of the host is around x , then please set low_depth=x/4 and high_depth=x*[3 or 5]

  For trio-binning-based detection:

        --threshold1        minimum of parental-specific kmer density (default 0.001)
                            for ONT reads(error rate~=15%), we recommand 0.001.
                            for PacBio reads(error rate<5%), we recommand 0.002-0.005.

        --threshold2        minimum of shared kmer density (default 0.1)

  For BGMM-clustering-based detection:

        --cluster           (1/0) use clustering or not. default(0)

        --shortest          length threshold for clustering ( default 5000 )
                            only reads with lenghs > the shortest will be used for clustering.
                            shorter reads ( <=5k ) often create noisy points and hamper accurate clustering.

        --loop              number of BGMM clustering for consensus (default 10)

        --python3           python3 path ( default python3 )

        --seed              random seed ( default 42 )

        --help              print this usage message.
```
## Results

### trio-binning result is recorded in trio.result.txt

```
column1 :   read-name.
column2 :   trio-base-result    trio-based-formula result 1 means host and 0 means contamination.
column3 :   density of paternal specific markers
column4 :   density of maternal specific markers
column5 :   density of shared  markers
```

All read names filtered by trio-binning use pipeline's threshold1 and threshold2 parameteres are recorded in ```trio_only.filtered_readname.txt```

Use "./tools/tune_t1_t2.sh to test new parameters."¬

**There is no need to re-run the pipeline if you want to tune the threshold1 and threshold2 parameters**

I use ```seqkit grep -f trio_only.filtered_readname.txt input.fasta >output.fasta``` to extract corresponding host's sequences.
```seqkit``` can be installed by ```conda install -c bioconda seqkit```.·

Fell free to use any sequence fishing tools that you preferred to do the same thing.

### clustering result is recorded in final.result.txt

```
column1 :   read-name.
column2 :   trio-only.          (trio-based-formula result) 1 refers to host and 0 refers to contamination.
column3 :   host.               (suggested culstering result by hit-count>(loop-time/5)) 1 refers to host and 0 refers to contamination)
                                (recommand to use your own threshold based on the profile)
column4 :   hit-count           how often this read occurs in the best cluster.
column5 :   second-hit-count    how often this read occurs in the second-best cluster.
```

Warning: reads shorter than 5000bp will be discarded in ```final.result.txt```. However, you can still extract short host reads from ```trio_only.results.txt```

example of using fishSq to get final host's reads

```
./tools/fishSq.sh -f final.result.txt -t 4 -i mixture.fasta -o output.fasta
```

## How to run this pipeline step by step.

### first build markers

first use ```00.BuildTrioLib.sh``` to generate K21 k-mers.

### then calculate marker densities of each long read and filter them

Use ```11.GetTrioMatrix_tgs.sh```

Here we get trio-binning result in ```trio_only.result.txt``` and ```trio_only.filtered_readname.txt```

To get filtered reads based on trio binning, try ```seqkit grep -f trio_only.filtered_readname.txt input.fasta >output.fasta```

This step also genereates ```trio.4r.matrix``` for the following clustering.

### then cluster to find more reads [Optional]

#### prepare GC & 3mer matrix 

Use ```20.GetGCNmer.sh```

This step genereates ```gc.matrix``` for clustering.

#### run BGMM

run ```python3 bgm/main_logic.py -t trio.4r.matrix -m gc.matrix   1>cluter.result.txt 2>log ```

Here we get clustering result in ```cluter.result.txt```
We can get final.result.txt by ```tools/merge_result.sh```
To get the final filtered reads based on clustering, try ```tools/fishSq.sh```

## Clustering or not

Run clustering if the error rates of long reads are high (>5%), or the filtered data by trio-binning markers are not sufficient to support assembly.

The clustering step will detect more host reads (especially for those without any trio-binning markers ) but can introduce some sequences from closely-related symbionts because their genome signatures are similar to the host.

## How to tune parameters

### parameters for trio-binning
There are 4 parameters for trio-binning step : ```low_depth high_depth threshold1 threshold2```

-----------------------------------------

The ```low_depth``` and ```high_depth``` are used for generating parental and shared markers. Only markers with depths of [low_depth,high_depth] are used.

The pipeline supports automatic detection of those two thresholds by analysing the "Kmer depth"-"Count" graph.

2 situations that you cannot use automatic detection

1. If "Kmer depth"-"Count" graph is very complex, for instance, containing several well-matched peaks

2. If "Kmer depth"-"Count" graph contains no peaks.
   * maybe the depth of your host genome is very low because the majority of reads belong to contaminants.
   * maybe the genome sizes are huge for some contaminants.

-----------------------------------------
The ```threshold1``` and ```threshold2``` are used for trio-binning-based filtration.

The density graph ( either parent-specific marker densities (threshold1) or shared marker densities (threshold2) ) consisted by 2 part : 

1. Some (maybe sharp) noise peaks in the left ( close to zero )
2. Some gentle and wide peaks after the noise peaks.

We need to set the threshold1/2 after the just after those shark noise peaks.

To achive that, I usually check the density velocity information by :

```
# to get threshold2
>awk '{print int($10)}' trio_density.data.txt | sort -n | uniq -c | awk '{if(NR>1) print t prev/$1; t=$2;prev=$1;}'

0 0.00114416
1 0.0752022
2 0.654466
3 2.14107
4 2.1881
5 1.88301
6 1.5913
7 1.68106
8 1.94696
9 1.59053
10 0.947368 <------ the velocity trend to close to 1 now , so we set threshold2 as 10
11 0.747813
12 0.745652
13 0.788346
14 0.88142
15 0.937013  
16 1.05134
17 1.07006
18 0.900358
19 0.93
20 0.934579
21 0.957066
22 0.983001
23 0.965478
24 0.975704
25 1.09957
26 0.993366
27 1.12636
28 1.06744
29 1.13779
30 1.10584
31 1.17345
32 1.22251
33 1.10564
34 1.22954
35 1.2086
36 1.04027
37 1.24167
38 1.17264
39 1.00987
40 1.09353
41 1.02583
42 1
43 0.860317
44 0.813953
45 0.830472
46 0.889313
47 0.86755
48 0.963317
49 0.880618
50 0.87362
51 0.934633
52 0.924708
# to get threshod1
>awk '{if($8>$9) print int($8);else print int($9)}' trio_density.data.txt | sort -n | uniq -c | awk '{if(NR>1) print t prev/$1; t=$2;prev=$1;}'
```


### parameters for generating clustering result

There is only 1 parameter for clustering : 

The ```-t``` frequency threshold of reads found in best and second-best clusters for consensus.

For totally 10 loops of clustering for consensus,  ```-t 8``` or ```-t 9```  will give a high precision rate.

To obtain more host reads ( outliers far from the center of pre-detected host points ), set a lower t based on the histgram of read hits:

Try ```awk '{if(NR>1) print $4+$5;}' final.result.txt | sort -k2n | uniq -c``` and you got the data of read-hits histogram.

Assume you get
```
>awk '{if(NR>1) print $4+$5;}' final.result.txt | sort -k2n | uniq -c
10000 0
 8000 1
 1000 2
   50 3
  100 4
  300 5
  200 6
  500 7
 3000 8
 5000 9
 7000 10
>
```
data with hits=0,1,2 are the majority of contamination, while data with hits=8,9,10 are the majority of host. 
So try ```-t 2``` if you prefer a good recall rate with reasonable precision.
Or try ```-t 7``` if you prefer a very high precision rate.


__END__
