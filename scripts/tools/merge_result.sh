#!/bin/bash

if [[ $# -lt 2 ]] ; then
    echo "Usage : merge_results.sh trio_density.data.txt cluster_reuslt.txt shortest(default 5000)"
    exit 1
fi
Density_file=$1
Cluster_file=$2
shortest=5000
if [[ $# -gt 2 ]] ; then
    shortest=$3
fi
echo 'name' >name.txt
awk -v short=$shortest '{if($4>short)print $2;}' $Density_file >>name.txt
paste name.txt $Cluster_file >final.result.txt
