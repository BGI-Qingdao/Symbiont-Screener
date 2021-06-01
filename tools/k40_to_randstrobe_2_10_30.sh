#!/bin/bash

SPATH=`realpath $0`
SPATH=`dirname $SPATH`
K2S="$SPATH""/../main/kmer2strobemer"

if [[ ! -e  paternal.unique.filter.mer || ! -e maternal.unique.filter.mer  || ! -e common.mer ]] ; then 
    echo "input file not exist. please link paternal.unique.filter.mer  maternal.unique.filter.mer common.mer in current folder"
    exit 1
fi

if [[ ! -e '00.kmer2strobemer_done' ]]  ; then
    $K2S --nkmer 2 --ksize 10 --wsize 30 <paternal.unique.filter.mer >paternal.strobemer || exit 1
    $K2S --nkmer 2 --ksize 10 --wsize 30 <maternal.unique.filter.mer >maternal.strobemer || exit 1
    $K2S --nkmer 2 --ksize 10 --wsize 30 <common.mer >common.strobemer || exit 1
    awk '{if(FILENAME==ARGV[1]) {m[$1]=1;} else {if($1 in m) m[$1]=0;} }END{for( x in m ) { if(m[x]==1)  print x>"paternal_only.strobemer" ; else print x>"pc.mer";}}' paternal.strobemer maternal.strobemer common.strobemer &
    awk '{if(FILENAME==ARGV[1]) {m[$1]=1;} else {if($1 in m) m[$1]=0;} }END{for( x in m ) { if(m[x]==1)  print x>"maternal_only.strobemer" ; else print x>"mc.mer";}}' maternal.strobemer paternal.strobemer common.strobemer &
    wait

    cat common.strobemer pc.mer mc.mer >common.all.strobemer || exit 1
    awk '{m[$1]=1;}END{for(x in m) print x;}' common.all.strobemer > common_uniq.strobemer ||exit 1
    date >>'00.kmer2strobemer_done'
else
    echo "skip kmer2strobemer due to 00.kmer2strobemer_done exist!"
fi
