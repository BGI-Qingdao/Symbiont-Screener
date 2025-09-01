#!/bin/bash

###############################################################################
# usage function
###############################################################################
function usage(){
echo """
Usage    :
    ./sysc consensus_cluster_k21 [OPTION]

Options  :
        --min_hit           minimum hit in best and second best cluster
        --help              print this usage message.
"""
}

###############################################################################
# basic variables 
###############################################################################

MIN_HIT=""
CLUSTER_DATA=""

###############################################################################
# parse arguments
###############################################################################

if [[ $# == 0 ]] ; then 
    usage
    exit 0
fi

echo "CMD :$0 $*"
while [[ $# > 0 ]] 
do
    case $1 in
        "-h")
            usage
            exit 0
            ;;
        "--help")
            usage
            exit 0
            ;;
        "--min_hit")
            MIN_HIT=$2
            shift
            ;;
        "--cluster_result")
            CLUSTER_DATA=` realpath $2`
            shift
            ;;
        *)
            echo "invalid params : \"$1\" . exit ... "
            exit
        ;;
    esac
    shift
done

# print arguments
echo "consensus_cluster_k21.sh log : "
echo "    min_hit               : $MIN_HIT"

# sanity check
if [[  $MIN_HIT -lt 1 ||  $MIN_HIT -gt 10 ]] ; then
    echo "ERROR : arguments invalid ... exit!!! "
    exit 1
fi

date
echo "__START__"

mkdir -p 'step03.2.k21' && cd 'step03.2.k21'
if [[ ! -e '31.step_0_done' ]]  ; then
    # in final.result.txt
    # $1    name
    # $2    priori 
    # $3    host
    # $4    best-hit-counts 
    # $5    second best-hit-counts
    awk -v T=$MIN_HIT  '{if($4+$5>=T){print $1;} }' $CLUSTER_DATA >readname.min-hit"$MIN_HIT".txt || exit 1
    date >>'31.step_0_done'
else
    echo "skip get trio.matrix  due to 31.step_0_done exist"
fi

echo "__END__"

