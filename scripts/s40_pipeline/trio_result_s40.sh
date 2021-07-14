#!/bin/bash

###############################################################################
# usage function
###############################################################################
function usage(){
echo """
Usage    :
    ./trio_result_s40_sh [OPTION]

Options  :
        --threshold1        minimum density of parental kmer density(default 0.003)
        --threshold2        minimum density of shared kmer density(default 0.01)
"""
}

###############################################################################
# basic variables 
###############################################################################
THRESHOLD1=0.003
THRESHOLD2=0.1
prefix="output"
input=""
###############################################################################
# parse arguments
###############################################################################
if [[ $# == 0 || $1 == '-h' || $1 == "--help" ]] ; then 
    usage
    exit 0
fi
echo "CMD :$0 $*"
while [[ $# > 0 ]] 
do
    case $1 in
        "--threshold1")
            THRESHOLD1=$2
            shift
            ;;
        "--threshold2")
            THRESHOLD2=$2
            shift
            ;;
        "--input")
            input=`realpath $2`
            shift
            ;;

        *)
            echo "invalid params : \"$1\" . exit ... "
            exit
        ;;
    esac
    shift
done

if [[ $input == "" || ! -e $input ]] ; then
    echo "please do this action by sysc main pipeline ! exit ..."
    exit
fi

# print arguments
echo "trio_result_s40.sh log : "
echo "    threshold1            : $THRESHOLD1"
echo "    threshold2            : $THRESHOLD2"

date
echo "__START__"

mkdir -p 'step02.2.s40' && cd 'step02.2.s40'

if [[ ! -e '11.step_2_done' ]]  ; then
    # in trio_density.data.txt
    # $2    read name
    # $4    read length
    # $8    density*1000 of pat-only
    # $9    density*1000 of mat-only
    # $10   density*1000 of shared
    # trio-only result here
    awk -v T1=$THRESHOLD1 -v T2=$THRESHOLD2  '{
        if(($8>T1*1000|| $9>T1*1000) && $10>T2*1000) 
            priori=1 ;
        else
            priori=0 ;
        printf("%s\t%d\t%f\t%f\t%f\n",$2,priori,$8/1000,$9/1000,$10/1000);
    }'   $input >trio_only.result.txt || exit 1
    awk '{if($2==1)print $1;}' $input >trio_only.filtered_readname.txt
    date >>'11.step_2_done'
else
    echo "skip get trio_only.result.txt due to 11.step_2_done exist"
fi
cd -
echo "__END__"
