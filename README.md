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

## Using strobemer mode
If you would like to use strobemers (randstrobe(2,10,10,30)) as markers (recommended) for error-tolerant matching, then please follow [pipelines/s40_pipeline/](https://github.com/BGI-Qingdao/Symbiont-Screener/tree/master/pipelines/s40_pipeline) 

What is strobemer? See https://github.com/ksahlin/strobemers and https://github.com/BGI-Qingdao/strobemer_cpptest

## Using k-mer mode
If you would like to use common k-mers (k=21) as markers , then please follow [pipelines/k21_pipeline/](https://github.com/BGI-Qingdao/Symbiont-Screener/tree/master/pipelines/k21_pipeline)
