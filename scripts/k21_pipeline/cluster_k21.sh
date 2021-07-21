#!/bin/bash

###############################################################################
# usage function
###############################################################################
function usage(){
echo """
Usage    :
    ./sysc cluster_k21  [OPTION]

Options  :
        --offspring         Offspring sequence file.
                            gzip format file is supported but should end by '.gz' 
        --offspring_format  fasta/fastq (default fasta)

        --threshold1        minimum density of parental kmer density(default 0.003)
        --threshold2        minimum density of shared kmer density(default 0.01)
        --shortest          shortest for cluster (default 5000)

        --thread            thread num.
                            [ optional, default 8 threads. ]
        --help              print this usage message.
"""
}

###############################################################################
# basic variables 
###############################################################################

OFFSPRING_FORMAT="fasta"
OFFSPRING=""
TRIODATA=""
GC_NMER=""
BGM_MAIN=""

THRESHOLD1=0.003
THRESHOLD2=0.1
L_SHORTEST=5000

CPU=8
NMER=3

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
        "--threshold1")
            THRESHOLD1=$2
            shift
            ;;
        "--threshold2")
            THRESHOLD2=$2
            shift
            ;;
        "--trio_data")
            TRIODATA=`realpath $2`
            shift
            ;;
        "--dcl")
            DRAW_CL=`realpath $2`
            shift
            ;;
        "--thread")
            CPU=$2
            shift
            ;;
        "--offspring")
            offspring_file=`realpath $2`
            OFFSPRING="$offspring_file"" ""$OFFSPRING"
            shift
            ;;
        "--gc_nmer")
            GC_NMER=$2
            shift
            ;;
        "--bgmm")
            BGM_MAIN=$2
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
echo "cluster_k21.sh log : "
echo "    offspring format      : $OFFSPRING_FORMAT"
echo "    offspring input       : $OFFSPRING"
echo "    threshold1            : $THRESHOLD1"
echo "    threshold2            : $THRESHOLD2"
echo "    shortest              : $L_SHORTEST"
echo "    thread                : $CPU "

# sanity check
if [[  $CPU -lt 1 || $NMER <2  || -z "$OFFSPRING"  || -z "$TRIODATA" ]] ; then
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

mkdir -p 'step03.1.k21' && cd 'step03.1.k21'
if [[ ! -e '30.step_0_done' ]]  ; then
    # in trio_density.data.txt
    # $2    read name
    # $4    read length
    # $8    density*1000 of pat-only
    # $9    density*1000 of mat-only
    # $10   density*1000 of shared
    # trio-only result here
    awk -v T1=$THRESHOLD1 -v T2=$THRESHOLD2 -v short=$L_SHORTEST '{
       if($4>short){
            if(($8>T1*1000|| $9>T1*1000) && $10>T2*1000)
                priori=1 ;
            else
                priori=0 ;
            printf("%d\t%f\t%f\t%f\n",priori,$8/1000,$9/1000,$10/1000);
        }
    }' $TRIODATA >trio.4r.matrix || exit 1
    echo "read_name" >name.txt
    awk -v short=$L_SHORTEST '{if($4>short)print $2;}' $TRIODATA >>name.txt
    date >>'30.step_0_done'
else
    echo "skip get trio.matrix  due to 30.step_0_done exist"
fi

if [[ ! -e '30.step_1_done' ]]  ; then
    READ_ARG=""
    for x in $OFFSPRING
    do
        READ_ARG="--read $x ""$READ_ARG"
    done
    $GC_NMER $READ_ARG \
                  --kmer   $NMER \
                  --format $OFFSPRING_FORMAT \
                  --thread $CPU >gc_nmer.data.txt 2>gc_nmer.log || exit 1
    date >>'30.step_1_done'
else
    echo "skip GC_NMER  due to 30.step_2_done exist"
fi

if [[ ! -e '30.step_2_done' ]]  ; then
    awk -v short=$L_SHORTEST '{if($2>short)print $0 ;}' gc_nmer.data.txt >gc_nmer.cut.data.txt
    cut -f 3- gc_nmer.cut.data.txt > gc_nmer.matrix ||exit 1
    date >>'30.step_2_done'
else
    echo "skip extract gc_nmer.matrix due to 20.step_2_done exist"
fi

if [[ ! -e '30.step_3_done' ]]  ; then
    $BGM_MAIN -t trio.4r.matrix  \
        -m  gc_nmer.matrix \
        -r  42 \
        -l  10 >cluster_reuslt.txt 2>cluster.log || exit 1

    paste name.txt cluster_reuslt.txt >final.result.txt

    $DRAW_CL || exit 1
    echo "Result in final.result.txt"
    echo "  Format are \"name priori host best-hit-counts second best-hit-counts\""
    echo "NOTICE : reads length not great than $L_SHORTEST had beed excluded from final.result.txt"
    date >>'30.step_3_done'
else
    echo "skip cluster by bgm due to 30.step_3_done exist!"
fi

echo "__END__"
