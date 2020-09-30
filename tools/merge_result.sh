#!/bin/bash

cut -f 1 trio_density.data.txt >name.txt
paste name.txt cluster_reuslt.txt >final.result.txt 
