#!/bin/bash

###############################################################################
# usage function
###############################################################################
function usage(){
echo """
Usage    :
    ./00.BuildTrioLib.sh [OPTION]

Options  :
        --thread            thread num.
                            [ optional, default 8 threads. ]

        --memory            x (GB) of memory to used by meryl.
                            [ optional, default 50GB. ]

        --mer               mer-size
                            [ optional, default 21. ]

        --paternal          paternal NGS reads file in FASTQ format.
                            ( note : gzip format is NOT supported. )
                            [ optional, needed when --use_existing_libs off ]

        --maternal          maternal NGS reads file in FASTQ format.
                            ( note : gzip format is NOT supported. )
                            [ optional, needed when --use_existing_libs off ]

        --help              print this usage message.

Examples :

    ./00.BuildTrioLib.sh   --paternal father.fastq --maternal mater.fastq 

"""

}

###############################################################################
# basic variables 
###############################################################################
MER=21
CPU=8
MEMORY=50
PATERNAL=""
MATERNAL=""
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
        "--memory")
            MEMORY=$2
            shift
            ;;
        "--thread")
            CPU=$2
            shift
            ;;
        "--mer")
            MER=$2
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
        *)
            echo "invalid params : \"$1\" . exit ... "
            exit
        ;;
    esac
    shift
done

# print arguments
echo "ContamFilter.sh log : "
echo "    memory          : $MEMORY GB"
echo "    thread          : $CPU "
echo "    mer             : $MER "
echo "    paternal input  : $PATERNAL"
echo "    maternal input  : $MATERNAL"
echo "BuildTrioLib.sh in  : $SPATH"

MERYL=$SPATH"/merylbin/meryl"
# sanity check
if [[ $MEMORY -lt 1  || $CPU -lt 1  ]] ; then
    echo "ERROR : arguments invalid ... exit!!! "
    exit 1
fi

for x in $PATERNAL $MATERNAL
do
   if [[ ! -e $x ]] ; then 
       echo "ERROR : input file \"$x\" is not exist ! exit ..."
       exit 1
   fi
done
date
echo "__START__"

###############################################################################
# extract paternal.mer & maternal.mer & common.mer
###############################################################################
if [[ ! -e '00.step_1_done' ]]  ; then
    $MERYL count k=$MER output maternal.meryl threads=$CPU memory=$MEMORY $MATERNAL || exit 1
    $MERYL histogram maternal.meryl >maternal.hist || exit 1
    awk -f $SPATH"/find_bounds.awk" maternal.hist > maternal.bounds.txt || exit 1
    date >>'00.step_1_done'
else
    echo "skip meryl count maternal due to 00.step_1_done exist"
fi

if [[ ! -e '00.step_2_done' ]]  ; then
    $MERYL count k=$MER output paternal.meryl threads=$CPU memory=$MEMORY $PATERNAL || exit 1
    $MERYL histogram paternal.meryl >paternal.hist  || exit 1
    awk -f $SPATH"/find_bounds.awk" paternal.hist > paternal.bounds.txt || exit 1
    date >>'00.step_2_done'
else
    echo "skip meryl count paternal due to 00.step_2_done exist"
fi

MLOWER=`grep LOWER_INDEX maternal.bounds.txt| awk -F '=' '{print $2}'`
MUPPER=`grep UPPER_INDEX maternal.bounds.txt| awk -F '=' '{print $2}'`
PLOWER=`grep LOWER_INDEX paternal.bounds.txt| awk -F '=' '{print $2}'`
PUPPER=`grep UPPER_INDEX paternal.bounds.txt| awk -F '=' '{print $2}'`
echo "  the real used kmer-count bounds of maternal is [ $MLOWER , $MUPPER ] "
echo "  the real used kmer-count bounds of paternal is [ $PLOWER , $PUPPER ] "

MLOWER=$(($MLOWER-1))
MUPPER=$(($MUPPER+1))
if [[ ! -e '00.step_3_done' ]]  ; then
    $MERYL difference output mat_only.meryl threads=$CPU memory=$MEMORY maternal.meryl paternal.meryl
    $MERYL greater-than $MLOWER  mat_only.meryl output 'mat_only.gt'$MLOWER'.meryl' || exit  1
    $MERYL less-than $MUPPER  'mat_only.gt'$MLOWER'.meryl' output 'mat_only.gt'$MLOWER'.lt'$MUPPER'.meryl' || exit 1
    ln -s  'mat_only.gt'$MLOWER'.lt'$MUPPER'.meryl' mat_only.filtered.meryl
    $MERYL print mat_only.filtered.meryl | awk '{print $1}' >maternal.mer || exit  1
    date >>'00.step_3_done'
else
    echo "skip get maternal.mer due to 00.step_3_done exist"
fi


PLOWER=$(($PLOWER-1))
PUPPER=$(($PUPPER+1))
if [[ ! -e '00.step_4_done' ]]  ; then
    $MERYL difference output pat_only.meryl threads=$CPU memory=$MEMORY paternal.meryl maternal.meryl || exit 1
    $MERYL greater-than $PLOWER  pat_only.meryl output 'pat_only.gt'$PLOWER'.meryl' || exit  1
    $MERYL less-than $PUPPER  'pat_only.gt'$PLOWER'.meryl' output 'pat_only.gt'$PLOWER'.lt'$PUPPER'.meryl' || exit 1

    ln -s  'pat_only.gt'$PLOWER'.lt'$PUPPER'.meryl' pat_only.filtered.meryl
    $MERYL print pat_only.filtered.meryl | awk '{print $1}' >paternal.mer || exit 1
    date >>'00.step_4_done'
else
    echo "skip get paternal.mer due to 00.step_4_done exist"
fi

if [[ ! -e '00.step_5_done' ]]  ; then
    $MERYL intersect-sum output common.meryl threads=$CPU memory=$MEMORY paternal.meryl maternal.meryl  || exit 1
    CLOW=$((($PLOWER+$MLOWER)/2-1))
    $MERYL greater-than $CLOW output 'common.gt'$CLOW'.meryl' common.meryl || exit 1
    ln -s 'common.gt'$CLOW'.meryl' 'common.filtered.meryl'
    $MERYL print 'common.filtered.meryl'  | awk '{print $1}' >common.mer || exit 1
    date >>'00.step_5_done'
else
    echo "skip get common.mer due to 00.step_5_done exist"
fi
echo "__END__"
