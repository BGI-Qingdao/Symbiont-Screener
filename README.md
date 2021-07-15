# Symbiont-Screener

Separate host's sequences from contamination for TGS long reads by trio-binning markers and unsupervised clustering.

## Depedences

* make
* gcc/g++ (version>4.8.5 to support c++11)
* (optional) anaconda3 for numpy and sklearn packages if you use cluster function.

## How to download and install

```
git clone https://github.com/BGI-Qingdao/Symbiont-Screener  ./Symbiont-Screener
cd  ./Symbiont-Screener/main
make
```

## Usage

```
Usage : sysc <action> [options]

Actions:
  +---------+-----------------------+-----------------------+
  |stage    | s40 action name       | k21 action name       |
  +---------+-----------------------+-----------------------+
  |step01   | build_s40             | build_k21             |
  +---------+-----------------------+-----------------------+
  |step02.1 | density_s40           | density_k21           |
  |step02.2 | trio_result_s40       | trio_result_k21       |
  +---------+-----------------------+-----------------------+
  |step03.1 | cluster_s40           | cluster_k21           |
  |step03.2 | consensus_cluster_s40 | consensus_cluster_k21 |
  +---------+-----------------------+-----------------------+
  please run those actions step by step
```

What is strobemer? See https://github.com/ksahlin/strobemers and https://github.com/BGI-Qingdao/strobemer_cpptest

