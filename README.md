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

* example1: run sysc in strobemer mode with cluster. (recommend)

```
./Symbiont-Screener/easy-to-use_pipelines/sysc_strobmercluster_mode.sh --maternal maternal_mix_simngs.v2.fasta.gz \
                                                                       --paternal paternal_mix_simngs.v2.fasta.gz \
                                                                       --offspring  mix_simpb.fasta.gz
```
* example2: run sysc in strobemer only mode. (also recommend)

```
./Symbiont-Screener/easy-to-use_pipelines/sysc_strobmer_mode.sh        --maternal maternal_mix_simngs.v2.fasta.gz \
                                                                       --paternal paternal_mix_simngs.v2.fasta.gz \
                                                                       --offspring  mix_simpb.fasta.gz
```
* example3: run sysc in kmer mode with cluster.

```
./Symbiont-Screener/easy-to-use_pipelines/sysc_kmercluster_mode.sh     --maternal maternal_mix_simngs.v2.fasta.gz \
                                                                       --paternal paternal_mix_simngs.v2.fasta.gz \
                                                                       --offspring  mix_simpb.fasta.gz
```
* example4: run sysc in kmer only mode.

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

#### The workflows in example pipelines:

Each pipeline in eays-to-use_pipelines folder uses sysc commmand but has different workflow.


```
The four available workflows of sysc :
  +--------------------------------------------------------------------------------+-------------------------------+
  |                                   Workflows                                    |     Example pipeline          |
  +---+------------------------------------------------------------------------+---+-------------------------------+
  |   | -> build_s40 -> density_s40 -> trio_result_s40 ----------------------> |   |  sysc_strobmer_mode.sh        | 
  | S |                     |                                                  |   |                               |
  | T |                     +--------> cluster_s40 -> consensus_cluster_s40 -> | E |  sysc_strobmercluster_mode.sh |
  | A |                                                                        | N |                               |
  | R | -> build_k21 -> density_k21 -> trio_result_k21 ----------------------> | D |  sysc_kmer_mode.sh            |
  | T |                     |                                                  |   |                               |
  |   |                     +--------> cluster_k21 -> consensus_cluster_k21 -> |   |  sysc_kmercluster_mode.sh     |
  +---+------------------------------------------------------------------------+---+-------------------------------+

```
To to avoid redundant work, the pipelines in easy-to-use_pipelines folder only provide the required paremeters. However, sysc as the main entry, provide all detail parameters that you need. Therefor, we recommend you create your own pipeline based on them add modify the parameters as your wish.


## Common Q & A 

### How to re-run after program exit by error？

First of all, please read the log and err files to find out why program exit. Invalid or incomplete parameters are the main reasons.
The sysc command create tag file to log finished steps, therefor, just re-run the same command will automaticlly skip fininshed steps and continue.

### How to re-run with different parameters？

After you modify your pipelines, please move the rename folder so that program will create new folder and restart without the influence of previous data.

For example, if you only change the parameter of ```cluster_k21```, please ```rename/move/delete``` both ```step03.1.k21``` and ```step03.2.k21``` so that program will re-run both ```cluster_k21``` and ```consensus_cluster_k21```. In this case, please keep ```step01.k21``` and ```step02.1.k21``` unchanged so that program will re-use those data.

### How to finetune parameters ?

* The ngs depth parameters

Here are the most important parameters that you shall check.

Those five parameters exist in both ```build_s40``` and ```build_k21``` actions.
```
  --auto_bounds
  --m-lower    
  --m-upper    
  --p-lower    
  --p-upper    
```
By default ```--auto_bounds=1``` so that program will infer the bounds automaticlly. 

An interactive html named ```kmer_frequency.html``` will be created in ```step02.1.k21``` or ```step02.1.s40``` folder like below :

![image](https://user-images.githubusercontent.com/8720584/198863141-3c87ff07-e5a4-47ff-9c4d-367ff1a55264.png)

If the host reads are not the majority of all your reads, you should set ```--auto_bounds=0``` and custom your depth threshold base on this html image.

* The minimum number of best-hit for consensus clustering.

In with cluster mode, an interactive html image named ```cluster.html``` will be created in ```step03.1.s40``` or ```step03.1.k21``` folder like below :

![image](https://user-images.githubusercontent.com/8720584/198863720-c4622262-85d0-4cf9-95d8-48ba69c616d0.png)

Please set the ```--min_hit``` paramter by the guide of this image.

* Other parameters

Please check our detailed usage for other paremeter like the threads, memory control ,etc... 

## What is strobemer?

See https://github.com/ksahlin/strobemers and https://github.com/BGI-Qingdao/strobemer_cpptest

## Contact us

Raising an issue is always recommended.

Email to xumengyang@genomics.com, guolidong@genomics.com, and qiyanwei@genomics.com

## Enjoy
