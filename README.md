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

## How to download and compile

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

### test01 : run sysc in strobemer mode with cluster

```
./Symbiont-Screener/easy-to-use_pipelines/sysc_strobmer_mode.sh        --maternal maternal_mix_simngs.v2.fasta.gz \
                                                                       --paternal paternal_mix_simngs.v2.fasta.gz \
                                                                       --offspring  mix_simpb.fasta.gz
```
### test02 : run sysc in strobemer mode without cluster

```
./Symbiont-Screener/easy-to-use_pipelines/sysc_strobmercluster_mode.sh --maternal maternal_mix_simngs.v2.fasta.gz \
                                                                       --paternal paternal_mix_simngs.v2.fasta.gz \
                                                                       --offspring  mix_simpb.fasta.gz
```
### test03 : run sysc in kmer mode without cluster

```
./Symbiont-Screener/easy-to-use_pipelines/sysc_kmercluster_mode.sh     --maternal maternal_mix_simngs.v2.fasta.gz \
                                                                       --paternal paternal_mix_simngs.v2.fasta.gz \
                                                                       --offspring  mix_simpb.fasta.gz
```
### test04 : run sysc in kmer mode without cluster

```
./Symbiont-Screener/easy-to-use_pipelines/sysc_kmer_mode.sh            --maternal maternal_mix_simngs.v2.fasta.gz \
                                                                       --paternal paternal_mix_simngs.v2.fasta.gz \
                                                                       --offspring  mix_simpb.fasta.gz
```

## Working flow of sysc pipelines

All pipelines in eays-to-use_pipelines folder use sysc commmand with different  workflow.

### The sysc command
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

### Details of workflows:

```
The four available workflows of sysc :
  +--------------------------------------------------------------------------------------------+-------------------------------+
  |                                   Available workflows                                      |     Corresponding pipeline    |
  +-------+-----------------------------------------------------------------------------+------+-------------------------------+
  |       | --> build_s40 --> density_s40 --> cluster_s40 --> consensus_cluster_s40 --> |      |  sysc_strobmercluster_mode.sh | 
  |       |                      |                                                      |      |                               |
  |       |                      +--------------------------> trio_result_s40 --------> |      |  sysc_strobmer_mode.sh        |
  | START |                                                                             | END  |                               |
  |       | --> build_k21 --> density_k21 --> cluster_k21 --> consensus_cluster_k21 --> |      |  sysc_kmercluster_mode.sh     |
  |       |                      |                                                      |      |                               |
  |       |                      +--------------------------> trio_result_k21 --------> |      |  sysc_kmercluster_mode.sh     |
  +-------+-----------------------------------------------------------------------------+------+-------------------------------+

```
<img src="https://user-images.githubusercontent.com/38022049/198184371-d9bc9d44-2bde-45ba-ac87-58e5bc03896f.png" width="600">

## Finetune paremeters

### min and max depth setting

### Trio density threshold

### minimam consensus result

## What is strobemer?

See https://github.com/ksahlin/strobemers and https://github.com/BGI-Qingdao/strobemer_cpptest

