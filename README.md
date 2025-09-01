# Symbiont-Screener

## Description

Symbiont-Screener is a reference-free approach to identifying high-confidence host's long reads from symbionts and contaminants and overcoming the low sequencing accuracy according to a trio-based screening model.

<img src="https://user-images.githubusercontent.com/38022049/198184294-45610387-79dc-4860-a3b8-4b82315f0b42.png" width="800">

## Citing Symbiont-Screener

Symbiont-Screener: A reference-free tool to separate host sequences from symbionts for error-prone long reads

If you find MetaTrass is useful in your research, please cite：

* Online published: [Xu M, Guo L, Qi Y, Shi C, Liu X, Chen J, Han J, Deng L, Liu X and Fan G (2023) Symbiont-screener: A reference-free tool to separate host sequences from symbionts for error-prone long reads. Front. Mar. Sci. 10:1087447. doi: 10.3389/fmars.2023.1087447](https://www.frontiersin.org/articles/10.3389/fmars.2023.1087447/full)

* Preprint: [Xu M, Guo L, Shi C, Liu X, Chen J, Liu X and Fan G (2020) Symbiont-Screener: a reference-free filter to automatically separate host sequences and contaminants for long reads or co-barcoded reads by unsupervised clustering. bioRxiv 2020.10.26.354621. doi: 10.1101/2020.10.26.354621](https://www.biorxiv.org/content/10.1101/2020.10.26.354621v1)

## Dependencies
System requirement:
* make
* gcc/g++ (version>4.8.5 to support c++11)
* python3 in your environment with packages :
    * numpy
    * pandas
    * sklearn
    * plotly

Third-party software:
* [strobmers](https://github.com/ksahlin/strobemers)

## Download and compile

### from conda

```
conda install -c intel scikit-learn
conda install -c bioconda symbiontscreener
```
notice: easy-to-use_pipelines are not included in the conda version, please use sysc directly

### from source codes
 
```
git clone https://github.com/BGI-Qingdao/Symbiont-Screener  ./Symbiont-Screener
cd  ./Symbiont-Screener/sources
make
```

## Quick start with the test dataset

To run Symbiont-Screener, the minimum input dataset contains one paternal NGS read file, one maternal NGS read file, and one child‘s TGS read file.

The simulated PacBio CLR reads of the child with parental NGS reads can be downloaded from https://zenodo.org/record/7257538#.Y1vuaHZByUk :

```
maternal_mix_simngs.v2.fasta.gz
paternal_mix_simngs.v2.fasta.gz
mix_simpb.fasta.gz
```

* Usage example I: 
   * run Symbiont-Screener in **strobemer mode with clustering**. (**Recommended**:warning:)

```
   ./Symbiont-Screener/easy-to-use_pipelines/sysc_strobmercluster_mode.sh \
         --maternal maternal_mix_simngs.v2.fasta.gz \
         --paternal paternal_mix_simngs.v2.fasta.gz \
         --offspring  mix_simpb.fasta.gz
```
* Usage example II: 
   * run Symbiont-Screener in **strobemer mode without further clustering**.

```
   ./Symbiont-Screener/easy-to-use_pipelines/sysc_strobmer_mode.sh \
         --maternal maternal_mix_simngs.v2.fasta.gz \
         --paternal paternal_mix_simngs.v2.fasta.gz \
         --offspring  mix_simpb.fasta.gz
```
* Usage example III : 
   * run Symbiont-Screener in **k-mer mode with clustering**.

```
   ./Symbiont-Screener/easy-to-use_pipelines/sysc_kmercluster_mode.sh \ 
         --maternal maternal_mix_simngs.v2.fasta.gz \
         --paternal paternal_mix_simngs.v2.fasta.gz \
         --offspring  mix_simpb.fasta.gz
```
* Usage example IV: 
   * run Symbiont-Screener in **k-mer mode without further clustering**.

```
   ./Symbiont-Screener/easy-to-use_pipelines/sysc_kmer_mode.sh \ 
         --maternal maternal_mix_simngs.v2.fasta.gz \
         --paternal paternal_mix_simngs.v2.fasta.gz \
         --offspring  mix_simpb.fasta.gz
```

## Details of Symbiont-Screener

### The software design of Symbiont-Screener
<img src="https://user-images.githubusercontent.com/38022049/198184371-d9bc9d44-2bde-45ba-ac87-58e5bc03896f.png" width="600">

### sysc : the main entry of Symbiont-Screener.

* sysc usage:

   * type ```./Symbioint-Screener/sysc -h``` and get:

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

* module usage:
   * Each action has it's own usage, for example ```./Symbioint-Screener/sysc build_s40 -h``` get:

```
Usage    :
  ./sysc build_s40 [OPTION]

Build randstrobe(2,10,30,30) based on paternal and maternal NGS reads by jellyfish.

Options  :
  basic input:
        --paternal    paternal NGS read file in FASTA/FASTQ format.
                      file in gzip format can be accepted, but filename must end by ".gz".
        --maternal    maternal NGS read file in FASTA/FASTQ format.
                      file in gzip format can be accepted, but filename must end by ".gz".
  resources:
        --thread      thread number.
                      [ optional, default 8 threads. ]
        --size        initial hash table by jellyfish.
                      [ optional, default 1GB. ]
  detail configuration:
        --auto_bounds (0/1) automatically calculate lower and upper bounds based on k-mer analysis.
                      [ optional, default 1; ]
        --m-lower     assigned lower bound for maternal k-mer frequencies, which will ignore k-mers with count < m-lower.
                      [ optional, default 0. ]
        --m-upper     assigned upper bound for maternal k-mer frequencies, which will ignore k-mers with count > m-upper.
                      [ optional, default 0. ]
        --p-lower     assigned lower bound for paternal k-mer frequencies, which will ignore k-mers with count < p-lower.
                      [ optional, default 0. ]
        --p-upper     assigned upper bound for paternal k-mer frequencies, which will ignore k-mers with count > p-upper.
                      [ optional, default 0. ]
```

* Please see details of all other actions by ```sysc <action> -h```

   * The workflows in example pipelines:

     Each pipeline in easy-to-use_pipelines folder uses sysc command but has different workflows.


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
To avoid redundant work, the pipelines in the easy-to-use_pipelines folder only show basic parameters. However, sysc as the main entry, can accept more detailed parameters that you might need for each specific case. Therefore, we recommend that you create your own pipeline based on the pipelines in the easy-to-use_pipelines folder and modify these parameters.


## Common Q & A 

### How can I re-run the program after it exits by accident？

First of all, please read the log and err files to find out why the program exits. Invalid or incomplete parameters are the main reasons.
The sysc command creates tag files to log finished steps. Therefore, direct re-running the same command in the same path will automatically skip these finished steps and continue.

### How can I re-run the program but with different parameters？

After you modify your pipeline, please move or rename previous result files so that the program will create a new working folder in the current path and restart without any influence on the previous results.

For example, if you only change the parameter of ```cluster_k21```, please ```rename/move/delete``` both ```step03.1.k21``` and ```step03.2.k21``` so that program will re-run both ```cluster_k21``` and ```consensus_cluster_k21```. In this case, please keep ```step01.k21``` and ```step02.1.k21``` unchanged so that program will re-use those results.

### How can I finetune parameters?

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
By default, ```--auto_bounds=1``` so that program will infer the bounds automatically. 

An interactive html named ```kmer_frequency.html``` will be created in ```step02.1.k21``` or ```step02.1.s40``` folder as below :

![image](https://user-images.githubusercontent.com/8720584/198863141-3c87ff07-e5a4-47ff-9c4d-367ff1a55264.png)

If the host reads are not the majority of the mixed reads, you should set ```--auto_bounds=0``` and customize your depth threshold based on this html image.

* The minimum number of best-hit for consensus clustering.

In clustering mode, an interactive html image named ```cluster.html``` will be created in ```step03.1.s40``` or ```step03.1.k21``` folder like below :

![image](https://user-images.githubusercontent.com/8720584/198863720-c4622262-85d0-4cf9-95d8-48ba69c616d0.png)

Please set the ```--min_hit``` parameter by the guide of this image.

* Other parameters

Please check our detailed usage for other parameters like the threads, memory control, etc... 

## What is strobemer?

See https://github.com/ksahlin/strobemers and https://github.com/BGI-Qingdao/strobemer_cpptest

## Contact us

Raising an issue is always recommended.

Email to xumengyang@genomics.cn, guolidong@genomics.cn, and qiyanwei1@genomics.cn

## Enjoy
:)

