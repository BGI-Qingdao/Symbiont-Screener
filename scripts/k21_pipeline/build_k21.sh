#!/bin/bash

###############################################################################
# usage function
###############################################################################
function usage(){
echo """
Usage    :
    ./sysc build_k21  [OPTION]

Options  :
        --thread            thread num.
                            [ optional, default 8 threads. ]

        --memory            x (GB) of memory to used by meryl.
                            [ optional, default 50GB. ]

        --mer               mer-size
                            [ optional, default 21. ]

        --paternal          paternal NGS reads file in FASTQ format.

        --maternal          maternal NGS reads file in FASTQ format.

        --m-lower           maternal kmer frequency table will ignore kmers with count < m-lower.
                            [ optional, default 0. ]
        --m-upper           maternal kmer frequency table will ignore kmers with count > m-upper.
                            [ optional, default 0. ]
        --p-lower           paternal kmer frequency table will ignore kmers with count < p-lower.
                            [ optional, default 0. ]
        --p-upper           paternal kmer frequency table will ignore kmers with count > p-upper.
                            [ optional, default 0. ]
        --auto_bounds       (0/1) automatically calcuate lower and upper bounds based on kmer analysis.
                            [ optional, default 1; ]
        --help              print this usage message.

Examples :

    ./sysc build_k21   --paternal father.fastq --maternal mater.fastq 
"""
}

###############################################################################
# basic variables 
###############################################################################
MER=21
CPU=8
MEMORY=50
PLOWER=0
PUPPER=0
MLOWER=0
MUPPER=0
PATERNAL=""
MATERNAL=""
AUTO_BOUNDS=1

MERYL=''
DRAW_KF=''
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
        "--meryl")
            MERYL=$2
            shift
            ;;
        "--dkf")
            DRAW_KF=$2
            shift
            ;;
        "--memory")
            MEMORY=$2
            shift
            ;;
        "--thread")
            CPU=$2
            shift
            ;;
        "--m-lower")
            MLOWER=$2
            shift
            ;;
        "--m-upper")
            MUPPER=$2
            shift
            ;;
        "--p-lower")
            PLOWER=$2
            shift
            ;;
        "--p-upper")
            PUPPER=$2
            shift
            ;;
        "--mer")
            MER=$2
            shift
            ;;
        "--auto_bounds")
            AUTO_BOUNDS=$2
            shift
            ;;
        "--paternal")
            tmp=`realpath $2`
            PATERNAL=$tmp" ""$PATERNAL"
            shift
            ;;
        "--maternal")
            tmp=`realpath $2`
            MATERNAL=$tmp" ""$MATERNAL"
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
echo "build_k21.sh log : "
echo "    memory          : $MEMORY GB"
echo "    thread          : $CPU "
echo "    mer             : $MER "
echo "    paternal input  : $PATERNAL"
echo "    maternal input  : $MATERNAL"
echo "    lower(maternal) : $MLOWER"
echo "    upper(maternal) : $MUPPER"
echo "    lower(paternal) : $PLOWER"
echo "    upper(paternal) : $PUPPER"
echo "    auto_bounds     : $AUTO_BOUNDS"

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
mkdir -p 'step01.k21' && cd 'step01.k21'
###############################################################################
# extract paternal.mer & maternal.mer & common.mer
###############################################################################
if [[ ! -e '00.step_1_done' ]]  ; then
    $MERYL  threads=$CPU memory=$MEMORY count k=$MER output maternal.meryl $MATERNAL || exit 1
    $MERYL  threads=$CPU memory=$MEMORY histogram maternal.meryl >maternal.histo || exit 1
    #awk -f $SPATH"/find_bounds.awk" maternal.histo > maternal.bounds.txt || exit 1
    awk 'BEGIN{MIN=0;MIN_INDEX=0;MAX=0;MAX_INDEX=0;STATE=0;}{i=0+$1;c=0+$2;if(S==0 ) {if(MIN==0 || c<MIN) {MIN=c ;MIN_INDEX=i;}else {S=1;} }else{ if(MAX==0 || c>MAX) { MAX=c; MAX_INDEX=i;}}}END{up_bounds=0+3*MAX_INDEX-2*MIN_INDEX;LOWER_INDEX=MIN_INDEX+1;UPPER_INDEX=up_bounds-1;printf("MIN_INDEX=%d\nMAX_INDEX=%d\nLOWER_INDEX=%d\nUPPER_INDEX=%d\n",MIN_INDEX,MAX_INDEX,LOWER_INDEX,UPPER_INDEX);}' maternal.histo > maternal.bounds.txt || exit 1
    date >>'00.step_1_done'
else
    echo "skip meryl count maternal due to 00.step_1_done exist"
fi

if [[ ! -e '00.step_2_done' ]]  ; then
    $MERYL threads=$CPU memory=$MEMORY count k=$MER output paternal.meryl $PATERNAL || exit 1
    $MERYL threads=$CPU memory=$MEMORY histogram paternal.meryl >paternal.histo  || exit 1
    #awk -f $SPATH"/find_bounds.awk" paternal.histo > paternal.bounds.txt || exit 1
    awk 'BEGIN{MIN=0;MIN_INDEX=0;MAX=0;MAX_INDEX=0;STATE=0;}{i=0+$1;c=0+$2;if(S==0 ) {if(MIN==0 || c<MIN) {MIN=c ;MIN_INDEX=i;}else {S=1;} }else{ if(MAX==0 || c>MAX) { MAX=c; MAX_INDEX=i;}}}END{up_bounds=0+3*MAX_INDEX-2*MIN_INDEX;LOWER_INDEX=MIN_INDEX+1;UPPER_INDEX=up_bounds-1;printf("MIN_INDEX=%d\nMAX_INDEX=%d\nLOWER_INDEX=%d\nUPPER_INDEX=%d\n",MIN_INDEX,MAX_INDEX,LOWER_INDEX,UPPER_INDEX);}' paternal.histo > paternal.bounds.txt || exit 1
    $DRAW_KF || exit 1
    date >>'00.step_2_done'
else
    echo "skip meryl count paternal due to 00.step_2_done exist"
fi

if [[ ! -e '00.step_3_done' ]]  ; then
    if [[ $AUTO_BOUNDS == 1 ]] ; then 
        MLOWER=`grep LOWER_INDEX maternal.bounds.txt| awk -F '=' '{print $2}'`
        MUPPER=`grep UPPER_INDEX maternal.bounds.txt| awk -F '=' '{print $2}'`
        MLOWER=$(($MLOWER-1))
        MUPPER=$(($MUPPER+1))
    fi
    echo "  the real used kmer-count bounds of maternal is [ $MLOWER , $MUPPER ] "
    $MERYL threads=$CPU memory=$MEMORY difference output mat_only.meryl maternal.meryl paternal.meryl
    $MERYL threads=$CPU memory=$MEMORY greater-than $MLOWER  mat_only.meryl output 'mat_only.gt'$MLOWER'.meryl' || exit  1
    $MERYL threads=$CPU memory=$MEMORY less-than $MUPPER  'mat_only.gt'$MLOWER'.meryl' output 'mat_only.gt'$MLOWER'.lt'$MUPPER'.meryl' || exit 1
    ln -s  'mat_only.gt'$MLOWER'.lt'$MUPPER'.meryl' mat_only.filtered.meryl
    $MERYL threads=$CPU memory=$MEMORY print mat_only.filtered.meryl | awk '{print $1}' >maternal.mer || exit  1
    date >>'00.step_3_done'
else
    echo "skip get maternal.mer due to 00.step_3_done exist"
fi

if [[ ! -e '00.step_4_done' ]]  ; then
    if [[ $AUTO_BOUNDS == 1 ]] ; then 
        PLOWER=`grep LOWER_INDEX paternal.bounds.txt| awk -F '=' '{print $2}'`
        PUPPER=`grep UPPER_INDEX paternal.bounds.txt| awk -F '=' '{print $2}'`
        PLOWER=$(($PLOWER-1))
        PUPPER=$(($PUPPER+1))
    fi
    echo "  the real used kmer-count bounds of paternal is [ $PLOWER , $PUPPER ] "
    $MERYL threads=$CPU memory=$MEMORY difference output pat_only.meryl paternal.meryl maternal.meryl || exit 1
    $MERYL threads=$CPU memory=$MEMORY greater-than $PLOWER  pat_only.meryl output 'pat_only.gt'$PLOWER'.meryl' || exit  1
    $MERYL threads=$CPU memory=$MEMORY less-than $PUPPER  'pat_only.gt'$PLOWER'.meryl' output 'pat_only.gt'$PLOWER'.lt'$PUPPER'.meryl' || exit 1

    ln -s  'pat_only.gt'$PLOWER'.lt'$PUPPER'.meryl' pat_only.filtered.meryl
    $MERYL threads=$CPU memory=$MEMORY print pat_only.filtered.meryl | awk '{print $1}' >paternal.mer || exit 1
    date >>'00.step_4_done'
else
    echo "skip get paternal.mer due to 00.step_4_done exist"
fi

if [[ ! -e '00.step_5_done' ]]  ; then
    $MERYL threads=$CPU memory=$MEMORY intersect-sum output common.meryl paternal.meryl maternal.meryl  || exit 1
    CLOW=$((($PLOWER+$MLOWER)/2-1))
    $MERYL threads=$CPU memory=$MEMORY greater-than $CLOW output 'common.gt'$CLOW'.meryl' common.meryl || exit 1
    ln -s 'common.gt'$CLOW'.meryl' 'common.filtered.meryl'
    $MERYL threads=$CPU memory=$MEMORY print 'common.filtered.meryl'  | awk '{print $1}' >common.mer || exit 1
    date >>'00.step_5_done'
else
    echo "skip get common.mer due to 00.step_5_done exist"
fi
cd -

echo "__END__"
