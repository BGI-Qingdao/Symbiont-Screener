#!/bin/bash

###############################################################################
# usage function
###############################################################################
function usage(){
echo """
Usage    :
    ./sysc density_s40 [OPTION]

Options  :
        --offspring         Offspring sequence file.
                            gzip format file is supported but should end by '.gz' 
        --offspring_format  fasta/fastq 
                            [ optional, default fasta]
        --thread            thread num.
                            [ optional, default 8 threads. ]
"""
}

###############################################################################
# basic variables 
###############################################################################
CPU=8
OFFSPRING=""
OFFSPRING_FORMAT="fasta"
PATERNAL_MER=""
MATERNAL_MER=""
SHARED_MER=""
DENSITY_3LIB=""
DRAW_TD=""
###############################################################################
# parse arguments
###############################################################################
if [[ $# == 0 || $1 == "-h" || $1 == '--help' ]] ; then 
    usage
    exit 0
fi
echo "CMD :$0 $*"
while [[ $# > 0 ]] 
do
    case $1 in
        "--thread")
            CPU=$2
            shift
            ;;
        "--d3l")
            DENSITY_3LIB=$2
            shift
            ;;
        "--offspring")
            tmp_file=`realpath $2`
            OFFSPRING="$tmp_file"" ""$OFFSPRING"
            shift
            ;;
        "--paternal_mer")
            PATERNAL_MER=`realpath $2`
            shift
            ;;
        "--maternal_mer")
            MATERNAL_MER=`realpath $2`
            shift
            ;;
        "--shared_mer")
            SHARED_MER=`realpath $2`
            shift
            ;;
        "--offspring_format")
            OFFSPRING_FORMAT=$2
            shift
            ;;
        "--dtd")
            DRAW_TD=`realpath $2`
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
echo "density_s40 log : "
echo "    thread                : $CPU "
echo "    offspring format      : $OFFSPRING_FORMAT"
echo "    offspring    input    : $OFFSPRING"
echo "    pat-only-mer input    : $PATERNAL_MER"
echo "    mat-only-mer input    : $MATERNAL_MER"
echo "    shared-mer   input    : $SHARED_MER"

mkdir -p 'step02.1.s40' && cd 'step02.1.s40'
# sanity check
if [[  $CPU -lt 1 || -z $PATERNAL_MER || -z $MATERNAL_MER || -z $SHARED_MER || -z "$OFFSPRING" ]] ; then
    echo "ERROR : arguments invalid ... exit!!! "
    exit 1
fi

for x in $PATERNAL $MATERNAL $SHARED_MER $OFFSPRING
do
   if [[ ! -e $x ]] ; then 
       echo "ERROR : input file \"$x\" is not exist ! exit ..."
       echo "ERROR : please use the sysc main pipeline."
       exit 1
   fi
done

if [[ ! -e $DENSITY_3LIB ]] ; then 
    echo "ERROR : $DENSITY_3LIB is not exist! please run make first . exit ..."
    echo "ERROR : please use the sysc main pipeline."
    exit 1
fi
if [[ $OFFSPRING_FORMAT != 'fastq' && $OFFSPRING_FORMAT != "fasta" ]] ; then 
    echo "invalid offspring_format $OFFSPRING_FORMAT ! exit ..."
    exit 1
fi
date
echo "__START__"

if [[ ! -e '11.step_1_done' ]]  ; then
    READ_ARG=""
    for x in $OFFSPRING
    do
        READ_ARG="--read $x ""$READ_ARG"
    done
    $DENSITY_3LIB --hap $PATERNAL_MER \
                  --hap $MATERNAL_MER \
                  --hap $SHARED_MER $READ_ARG \
                  --nkmer 2 \
                  --ksize 10 \
                  --wsize 30 \
                  --format $OFFSPRING_FORMAT \
                  --thread $CPU >trio_density.data.txt 2>trio_density.log || exit 1
    $DRAW_TD || exit 1
    date >>'11.step_1_done'
else
    echo "skip density_3lib due to 11.step_1_done exist"
fi

# TODO : draw
cd -

