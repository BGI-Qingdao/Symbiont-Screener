#!/bin/bash

if [[ $# != 4 ]] ; then
    echo "Usage :   tune_t1_t2.sh new_threshold1 new_threshold2 trio_only.results.txt prefix"
    echo "Output:   prefix.trio_only.results.txt and prefix.trio_only.filtered_readname.txt"
    exit 1
fi

THRESHOLD1=$1
THRESHOLD2=$2
TO_txt=$3
prefix=$4

awk -v T1=$THRESHOLD1 -v T2=$THRESHOLD2  '{
    if(($3>T1*1000|| $4>T1*1000) && $4>T5*1000) 
        priori=1 ;
    else
        priori=0 ;
    printf("%s\t%d\t%f\t%f\t%f\n",$1,priori,$3/1000,$4/1000,$5/1000);
}'   $TO_txt >$prefix".trio_only.results.txt"

awk '{if($2==1)print $1;}' $prefix".trio_only.results.txt" >$prefix"trio_only.filtered_readname.txt"

echo "Done. result in "$prefix".trio_only.results.txt and "$prefix"trio_only.filtered_readname.txt"
