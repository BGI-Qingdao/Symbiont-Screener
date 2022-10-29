# Symbiont-Screener

Symbiont-Screener is a reference-free approach to identify high-confidence host’s long reads from symbionts and contaminants and overcome the low sequencing accuracy according to a trio-based screening model.

<img src="https://user-images.githubusercontent.com/38022049/198184294-45610387-79dc-4860-a3b8-4b82315f0b42.png" width="800">

## Depedences

* make
* gcc/g++ (version>4.8.5 to support c++11)
* python3 in your environment with packages :
    * numpy
    * pandas
    * sklearn
    * plotly

## Download and compile

```
git clone https://github.com/BGI-Qingdao/Symbiont-Screener  ./Symbiont-Screener
cd  ./Symbiont-Screener/sources
make
```
## Quick start with test dataset

To run Symbiont-Screener, the minimum input dataset contain one paternal NGS read file, one maternal NGS read file, and one son TGS read file.

The simulated PacBio CLR dataset can be download from https://zenodo.org/record/7257538#.Y1vuaHZByUk :

```
maternal_mix_simngs.v2.fasta.gz
paternal_mix_simngs.v2.fasta.gz
mix_simpb.fasta.gz
```

### example1: run sysc in strobemer mode with cluster. (recommend)

```
./Symbiont-Screener/easy-to-use_pipelines/sysc_strobmercluster_mode.sh --maternal maternal_mix_simngs.v2.fasta.gz \
                                                                       --paternal paternal_mix_simngs.v2.fasta.gz \
                                                                       --offspring  mix_simpb.fasta.gz
```
### example2: run sysc in strobemer only mode. (also recommend)

```
./Symbiont-Screener/easy-to-use_pipelines/sysc_strobmer_mode.sh        --maternal maternal_mix_simngs.v2.fasta.gz \
                                                                       --paternal paternal_mix_simngs.v2.fasta.gz \
                                                                       --offspring  mix_simpb.fasta.gz
```
### example3: run sysc in kmer mode with cluster.

```
./Symbiont-Screener/easy-to-use_pipelines/sysc_kmercluster_mode.sh     --maternal maternal_mix_simngs.v2.fasta.gz \
                                                                       --paternal paternal_mix_simngs.v2.fasta.gz \
                                                                       --offspring  mix_simpb.fasta.gz
```
### example4: run sysc in kmer only mode.

```
./Symbiont-Screener/easy-to-use_pipelines/sysc_kmer_mode.sh            --maternal maternal_mix_simngs.v2.fasta.gz \
                                                                       --paternal paternal_mix_simngs.v2.fasta.gz \
                                                                       --offspring  mix_simpb.fasta.gz
```

## Details of sysc

### The software design of sysc
<img src="https://user-images.githubusercontent.com/38022049/198184371-d9bc9d44-2bde-45ba-ac87-58e5bc03896f.png" width="600">

### sysc : the main entry of Symbiont Screener.

#### sysc usage

type ```./Symbioint-Screener/sysc -h``` and get:

```
Usage : sysc <action> [options]

Actions:
  +---------+-----------------------+-----------------------+
  |stage    | strobemer mode (s40)  | kmer mode (k21)       |
  +---------+-----------------------+-----------------------+
  |step01   | build_s40             | build_k21             |
  +---------+-----------------------+-----------------------+
  |step02.1 | density_s40           | density_k21           |
  |step02.2 | trio_result_s40       | trio_result_k21       |
  +---------+-----------------------+-----------------------+
  |step03.1 | cluster_s40           | cluster_k21           |
  |step03.2 | consensus_cluster_s40 | consensus_cluster_k21 |
  +---------+-----------------------+-----------------------+
```

Each action has it's own usage, for example ```./Symbioint-Screener/sysc build_s40 -h``` get:

```
Usage    :
  ./sysc build_s40 [OPTION]

Build randstrobe(2,10,30,30) based on paternal and maternal NGS reads by jellyfish.

Options  :
  basic input:
        --paternal    paternal NGS reads file in FASTA/FASTQ format.
                      file in gzip format can be accepted, but filename must end by ".gz".
        --maternal    maternal NGS reads file in FASTA/FASTQ format.
                      file in gzip format can be accepted, but filename must end by ".gz".
  resources:
        --thread      thread number.
                      [ optional, default 8 threads. ]
        --size        initial hash table by jellyfish.
                      [ optional, default 1GB. ]
  detail configuration:
        --auto_bounds (0/1) automatically calcuate lower and upper bounds based on kmer analysis.
                      [ optional, default 1; ]
        --m-lower     maternal kmer frequency table will ignore kmers with count < m-lower.
                      [ optional, default 0. ]
        --m-upper     maternal kmer frequency table will ignore kmers with count > m-upper.
                      [ optional, default 0. ]
        --p-lower     paternal kmer frequency table will ignore kmers with count < p-lower.
                      [ optional, default 0. ]
        --p-upper     paternal kmer frequency table will ignore kmers with count > p-upper.
                      [ optional, default 0. ]
```

Please see details of all other actions by ```sysc <action> -h```

### The workflows of pipelines:

Each pipeline in eays-to-use_pipelines folder uses sysc commmand but has different workflow.


```
The four available workflows of sysc :
  +---------------------------------------------------------------------------------------+-------------------------------+
  |                                       Workflows                                       |     Example pipeline          |
  +-------+------------------------------------------------------------------------+------+-------------------------------+
  |       | -> build_s40 -> density_s40 -> trio_result_s40 ----------------------> |      |  sysc_strobmer_mode.sh        | 
  |       |                                      |                                 |      |                               |
  |       |                                cluster_s40 -> consensus_cluster_s40 -> |      |  sysc_strobmercluster_mode.sh |
  | START |                                                                        | END  |                               |
  |       | -> build_k21 -> density_k21 -> trio_result_k21 ----------------------> |      |  sysc_kmer_mode.sh            |
  |       |                                      |                                 |      |                               |
  |       |                                cluster_k21 -> consensus_cluster_k21 -> |      |  sysc_kmercluster_mode.sh     |
  +-------+------------------------------------------------------------------------+------+-------------------------------+

```

## How to finetune paremeters ?

To to avoid redundant work, the pipelines in easy-to-use_pipelines folder only provide the required paremeters. 

However, sysc as the main entry, provide all detail parameters that you need.

Therefor, we recommend you create your own pipeline based on them add modify the parameters as your wish.

## What is strobemer?

See https://github.com/ksahlin/strobemers and https://github.com/BGI-Qingdao/strobemer_cpptest

## Contact us

Raising an issue is always recommended.

Email to xumengyang@genomics.com, guolidong@genomics.com, and qiyanwei@genomics.com

## Enjoy
