#!/bin/bash

function usage(){
echo """
Usage   :
    ./Symbiont-Screener.sh [options]

Options :

        --paternal          paternal NGS reads file in FASTQ format.

        --maternal          maternal NGS reads file in FASTQ format.

        --offspring         Offspring sequence file.
                            gzip format file is supported but should end by '.gz' 

        --offspring_format  fasta/fastq (default fasta)

        --threshold1        minimum density of parental kmer density (default 0.001)

        --threshold2        minimum density of shared kmer density (default 0.1)

        --kmer               kmer-size (default 21. ]

        --nmer              nmer for gc_nmer(default 2)

        --sequence_platform tgs/stlfr (default tgs)

        --loop              loop number of BGMM (default 30) 

        --thread            thread num.
                            [ optional, default 8 threads. ]

        --memory            x (GB) of memory to used by meryl.
                            [ optional, default 50GB. ]

        --python3           PATH to python3 file from anaconda3 ( default python3 )

        --help              print this usage message.
"""
}

PATERNAL=
MATERNAL=
OFFSPRING=

PYTHON3='python3'
THRESHOLD1=0.001
THRESHOLD2=0.1
OFFSPRING_FORMAT='fasta'
PLAT='tgs'
KMER=21
NMER=2
CPU='8'
MEMORY=50
MEMORY=100
LOOP=30
#LOW_HIT=0.2

SPATH=`dirname $0`
STEP0=$SPATH/00.BuildTrioLib.sh
STEP11=$SPATH/11.GetTrioMatrix_tgs.sh
STEP12=$SPATH/12.GetTrioMatrix_stlfr.sh
STEP2=$SPATH/20.GetGCNmer.sh
BGM_MAIN=$SPATH/bgm/main_logic.py

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
        "--memory")
            MEMORY=$2
            shift
            ;;
        "--thread")
            CPU=$2
            shift
            ;;
        "--mer")
            KMER=$2
            shift
            ;;
        "--paternal")
            PATERNAL=$2" "$PATERNAL
            shift
            ;;
        "--maternal")
            MATERNAL=$2" "$MATERNAL
            shift
            ;;
        "--offspring")
            OFFSPRING="$2"" ""$OFFSPRING"
            shift
            ;;
        "--offspring_format")
            OFFSPRING_FORMAT=$2
            shift
            ;;
        "--nmer")
            NMER=$2
            shift
            ;;
        "--sequence_platform")
            PLAT=$2
            shift
            ;;
        "--loop")
            LOOP=$2
            shift
            ;;
        "--python3")
            PYTHON3=$2
            shift
            ;;
        *)
            echo "invalid params : \"$1\" . exit ... "
            exit
        ;;
    esac
    shift
done

###############################################################################
# main logic
###############################################################################
$STEP0 --paternal "$PATERNAL" \
    --maternal "$MATERNAL" \
    --mer $KMER \
    --thread $CPU \
    --memory $MEMORY  || exit 1 

if [[ $PLAT == 'tgs' ]] ; then 
    $STEP11 --paternal_mer paternal.mer \
        --maternal_mer maternal.mer \
        --shared_mer common.mer  \
        --offspring "$OFFSPRING" \
        --offspring_format $OFFSPRING_FORMAT  \
        --threshold1 $THRESHOLD1 \
        --threshold2 $THRESHOLD2 \
        --thread $CPU || exit 1
elif [[ $PLAT == 'stlfr' ]] ; then
    $STEP12 --paternal_mer paternal.mer \
        --maternal_mer maternal.mer \
        --shared_mer common.mer  \
        --offspring "$OFFSPRING" \
        --offspring_format $OFFSPRING_FORMAT  \
        --threshold1 $THRESHOLD1 \
        --threshold2 $THRESHOLD2 \
        --thread $CPU || exit 1
fi

$STEP2 --nmer $NMER \
       --sequence_platform $PLAT \
       --offspring "$OFFSPRING" \
       --offspring_format $OFFSPRING_FORMAT  \
       --thread $CPU || exit 1

if [[ ! -e '30.step_1_done' ]]  ; then
    $PYTHON3  $BGM_MAIN -t trio.4r.matrix  \
                        -m  gc_nmer.matrix \
                        -l $LOOP >cluster_reuslt.txt 2>cluster.log || exit 1
    date >>'30.step_1_done'
else
    echo "skip cluster by bgm due to 30.step_1_done exist!"
fi

if [[ ! -e '30.step_2_done' ]]  ; then
    cut -f 1 trio_density.data.txt >name.txt
    paste name.txt cluster_reuslt.txt >final.result.txt 

    echo "Result in final.result.txt"
    echo "  Format are \"name priori host hit-counts\""
    date >>'30.step_2_done'
else
    echo "skip get final result due to 30.step_2_done exist!"
fi
echo "__ALL DONE__"
