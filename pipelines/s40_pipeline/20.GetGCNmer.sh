#!/bin/bash

###############################################################################
# usage function
###############################################################################
function usage(){
echo """
Usage    :
    ./21.GetGCNmer.sh [OPTION]

Options  :
        --thread            thread num.
                            [ optional, default 8 threads. ]

        --offspring         Offspring sequence file.
                            gzip format file is supported but should end by '.gz' 
        --offspring_format  fasta/fastq (default fasta)
        --shortest          shortest for cluster (default 5000)
        --nmer              nmer for gc_nmer(default 2)
        --help              print this usage message.
Examples :

    ./11.GetTrioMatrix.sh   --paternal_mer paternal.mer  --maternal_mer maternal.mer --shared_mer common.mer

"""

}

###############################################################################
# basic variables 
###############################################################################
CPU=8
NMER=3
OFFSPRING_FORMAT="fasta"
OFFSPRING=""
L_SHORTEST=5000
SPATH=`dirname $0`
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
        "--thread")
            CPU=$2
            shift
            ;;
        "--offspring")
            OFFSPRING="$2"" ""$OFFSPRING"
            shift
            ;;
        "--nmer")
            NMER=$2
            shift
            ;;
        "--shortest")
            L_SHORTEST=$2
            shift
            ;;
        "--offspring_format")
            OFFSPRING_FORMAT=$2
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
echo "ContamFilter.sh log : "
echo "    thread                : $CPU "
echo "    nmer                  : $NMER"
echo "    offspring format      : $OFFSPRING_FORMAT"
echo "    offspring    input    : $OFFSPRING"
echo "    shortest              : $L_SHORTEST"
echo "BuildTrioLib.sh in  : $SPATH"

GC_NMER=$SPATH"/main/gc_nmer"
# sanity check
if [[  $CPU -lt 1 || $NMER <2  || -z "$OFFSPRING" ]] ; then
    echo "ERROR : arguments invalid ... exit!!! "
    exit 1
fi

if [[ $NMER > 4 ]] ; then 
    echo "WARN  : your use nmer=$NMER which will print a huge output file !!! "
fi

for x in $OFFSPRING
do
   if [[ ! -e $x ]] ; then 
       echo "ERROR : input file \"$x\" is not exist ! exit ..."
       exit 1
   fi
done

if [[ ! -e $GC_NMER ]] ; then 
    echo "ERROR : $GC_NMER is not exist! please run make first . exit ..."
    exit 1
fi
if [[ $OFFSPRING_FORMAT != 'fastq' && $OFFSPRING_FORMAT != "fasta" ]] ; then 
    echo "invalid offspring_format $OFFSPRING_FORMAT ! exit ..."
    exit 1
fi

date
echo "__START__"
GC_NMER_EXE=$GC_NMER

if [[ ! -e '20.step_1_done' ]]  ; then
    READ_ARG=""
    for x in $OFFSPRING
    do
        READ_ARG="--read $x ""$READ_ARG"
    done
    $GC_NMER_EXE  $READ_ARG \
	              --kmer   $NMER \
                  --format $OFFSPRING_FORMAT \
                  --thread $CPU >gc_nmer.data.txt 2>gc_nmer.log || exit 1
    date >>'20.step_1_done'
else
    echo "skip GC_NMER  due to 20.step_2_done exist"
fi

if [[ ! -e '20.step_2_done' ]]  ; then
    awk -v short=$L_SHORTEST '{if($2>short)print $0 ;}' gc_nmer.data.txt >gc_nmer.cut.data.txt
    cut -f 3- gc_nmer.cut.data.txt > gc_nmer.matrix ||exit 1
    date >>'20.step_2_done'
else
    echo "skip extract gc_nmer.matrix due to 20.step_2_done exist"
fi

echo "__END__"
