# Symbiont-Screener

Separate host's sequences from contamination for TGS long reads by trio-binning markers and unsupervised clustering.

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

## Usage

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


The four available workflows of sysc :
  +-------+-----------------------------------------------------------------------------+------+
  |       | --> build_s40 --> density_s40 --> cluster_s40 --> consensus_cluster_s40 --> |      | 
  |       |                      |                                                      |      |
  |       |                      +--------------------------> trio_result_s40 --------> |      |
  | START |                                                                             | END  |
  |       | --> build_k21 --> density_k21 --> cluster_k21 --> consensus_cluster_k21 --> |      |
  |       |                      |                                                      |      |
  |       |                      +--------------------------> trio_result_k21 --------> |      |
  +-------+-----------------------------------------------------------------------------+------+

```

## How to run the example data

### prepare : generate the toy data

```
#TODO

```

### test01 : run sysc in strobemer mode with cluster

```
#TODO
```
### test02 : run sysc in strobemer mode without cluster

```
#TODO
```
### test03 : run sysc in kmer mode without cluster

```
#TODO
```
### test04 : run sysc in kmer mode without cluster

```
#TODO
```

## What is strobemer?

See https://github.com/ksahlin/strobemers and https://github.com/BGI-Qingdao/strobemer_cpptest

