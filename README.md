# Symbiont-Screener
Separate host's sequences from contamination for TGS long reads by trio-binning markers and unsupervised clustering.

## Using strobemers
If you would like to use strobemers (randstrobe(2,10,10,30)) as markers (recommended) for error-tolerant matching, then see [pipelines/s40_pipeline/](https://github.com/BGI-Qingdao/Symbiont-Screener/tree/master/pipelines/s40_pipeline) 

what is strobemer ? see https://github.com/ksahlin/strobemers and https://github.com/BGI-Qingdao/strobemer_cpptest

## Using k-mers
If you would like to use common k-mers (k=21) as markers ( implemented for comparison with strobemers, no longer maintain ), then see [pipelines/k21_pipeline/](https://github.com/BGI-Qingdao/Symbiont-Screener/tree/master/pipelines/k21_pipeline)
