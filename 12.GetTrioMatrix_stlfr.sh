#!/bin/bash

###############################################################################
# usage function
###############################################################################
function usage(){
echo """
Usage    :
    ./11.GetTrioMatrix_stlfr.sh [OPTION]

Options  :
        --thread            thread num.
                            [ optional, default 8 threads. ]

        --paternal_mer      Pat-only kmer library
        --maternal_mer      Mat-only kmer library
        --shared_mer        Shared kmer library

        --offspring         Offspring sequence file.
                            gzip format file is supported but should end by '.gz' 
        --offspring_format  fasta/fastq (default fasta)

        --threshold1        minimum density of parental kmer density
        --threshold2        minimum density of shared kmer density

        --help              print this usage message.

Examples :

    ./11.GetTrioMatrix.sh   --paternal_mer paternal.mer  --maternal_mer maternal.mer --shared_mer common.mer --offspring  split_reads.1.fq --offspring  split_reads.2.fq --offspring_format fastq

"""

}

###############################################################################
# basic variables 
###############################################################################
CPU=8
PATERNAL_MER=""
MATERNAL_MER=""
SHARED_MER=""
THRESHOLD1=0.001
THRESHOLD2=0.1
OFFSPRING=""
OFFSPRING_FORMAT="fasta"

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
        "--paternal_mer")
            PATERNAL_MER=$2
            shift
            ;;
        "--maternal_mer")
            MATERNAL_MER=$2
            shift
            ;;
        "--shared_mer")
            SHARED_MER=$2
            shift
            ;;
        "--offspring_format")
            OFFSPRING_FORMAT=$2
            shift
            ;;
        "--threshold1")
            THRESHOLD1=$2
            shift
            ;;
        "--threshold2")
            THRESHOLD2=$2
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
echo "    offspring format      : $OFFSPRING_FORMAT"
echo "    offspring    input    : $OFFSPRING"
echo "    pat-only-mer input    : $PATERNAL"
echo "    mat-only-mer input    : $MATERNAL"
echo "    shared-mer   input    : $SHARED_MER"
echo "12.GetTrioMatrix_stlfr.sh in  : $SPATH"

DENSITY_3LIB=$SPATH"/main/density_3lib_stlfr"
# sanity check
if [[  $CPU -lt 1 || -z $PATERNAL_MER || -z $MATERNAL_MER || -z $SHARED_MER || -z "$OFFSPRING" ]] ; then
    echo "ERROR : arguments invalid ... exit!!! "
    exit 1
fi

for x in $PATERNAL $MATERNAL $SHARED_MER $OFFSPRING
do
   if [[ ! -e $x ]] ; then 
       echo "ERROR : input file \"$x\" is not exist ! exit ..."
       exit 1
   fi
done

if [[ ! -e $DENSITY_3LIB ]] ; then 
    echo "ERROR : $DENSITY_3LIB is not exist! please run make first . exit ..."
    exit 1
fi
if [[ $OFFSPRING_FORMAT != 'fastq' && $OFFSPRING_FORMAT != "fasta" ]] ; then 
    echo "invalid offspring_format $OFFSPRING_FORMAT ! exit ..."
    exit 1
fi
date
echo "__START__"

###############################################################################
# extract paternal.mer & maternal.mer & common.mer
###############################################################################
if [[ ! -e '12.step_1_done' ]]  ; then
    READ_ARG=""
    for x in $OFFSPRING
    do
        READ_ARG="--read $x ""$READ_ARG"
    done
    $DENSITY_3LIB --hap0 $PATERNAL_MER \
                  --hap1 $MATERNAL_MER \
                  --hap2 $SHARED_MER \
                  $READ_ARG \
                  --format $OFFSPRING_FORMAT \
                  --thread $CPU >trio_density_stlfr.data.txt 2>trio_density_stlfr.log || exit 1
    date >>'12.step_1_done'
else
    echo "skip density_3lib due to 12.step_1_done exist"
fi

if [[ ! -e '12.step_2_done' ]]  ; then
    # $2    total kmer
    # $3    count of pat-only
    # $4    count of mat-only
    # $5    count of shared
    awk '{
            if(($3/$2+$4/$2)>=T1 && $5/$2>T2) 
                priori=1 ;
            else
                priori=0 ;
            printf("%d\t%f\t%f\t%f\n",priori,$3/$2,$4/$2,$5/$2);
         }'  T1=$THRESHOLD1 T2=$THRESHOLD2 trio_density_stlfr.data.txt >trio.4r.stlfr.matrix || exit 1
    date >>'12.step_2_done'
else
    echo "skip meryl count paternal due to 12.step_2_done exist"
fi

echo "__END__"
