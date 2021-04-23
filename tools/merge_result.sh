#!/bin/bash

if [[ $# != 2 ]] ; then
    echo "Usage : merge_results.sh trio_density.data.txt cluster_reuslt.txt "
    exit 1
fi
D=$1
C=$2
echo 'name' >name.txt
cut -f 2 $D >>name.txt
paste name.txt $C >final.result.txt 
