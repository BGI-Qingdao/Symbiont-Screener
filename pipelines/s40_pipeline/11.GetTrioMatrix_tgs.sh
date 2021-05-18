#!/bin/bash

###############################################################################
# usage function
###############################################################################
function usage(){
echo """
Usage    :
    ./11.GetTrioMatrix.sh [OPTION]

Options  :
        --thread            thread num.
                            [ optional, default 8 threads. ]

        --paternal_mer      Pat-only kmer library
        --maternal_mer      Mat-only kmer library
        --shared_mer        Shared kmer library

        --offspring         Offspring sequence file.
                            gzip format file is supported but should end by '.gz' 
        --offspring_format  fasta/fastq (default fasta)
        --shortest          shortest for cluster (default 5000)
        --threshold1        minimum density of parental kmer density
        --threshold2        minimum density of shared kmer density

        --help              print this usage message.

Examples :

    ./11.GetTrioMatrix.sh   --paternal_mer paternal.mer \\
                            --maternal_mer maternal.mer \\
                            --shared_mer common.mer \\
                            --offspring test.fa
"""

}

###############################################################################
# basic variables 
###############################################################################
CPU=8
PATERNAL_MER=""
MATERNAL_MER=""
SHARED_MER=""
THRESHOLD1=0.002
THRESHOLD2=0.1
OFFSPRING=""
OFFSPRING_FORMAT="fasta"
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
        "--shortest")
            L_SHORTEST=$2
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
echo "    pat-only-mer input    : $PATERNAL_MER"
echo "    mat-only-mer input    : $MATERNAL_MER"
echo "    shared-mer   input    : $SHARED_MER"
echo "    threshold1            : $THRESHOLD1"
echo "    threshold2            : $THRESHOLD2"
echo "    shortest              : $L_SHORTEST"
echo "11.GetTrioMatrix_tgs.sh in:  $SPATH"

DENSITY_3LIB=$SPATH"/main/density_3lib_strobemer"
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
    date >>'11.step_1_done'
else
    echo "skip density_3lib due to 11.step_1_done exist"
fi

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
    }'   trio_density.data.txt >trio_only.result.txt || exit 1
    awk '{if($2==1)print $1;}' trio_only.result.txt >trio_only.filtered_readname.txt
    # for cluster :
    awk -v T1=$THRESHOLD1 -v T2=$THRESHOLD2 -v short=$L_SHORTEST '{
            if($4>short){
                if(($8>T1*1000|| $9>T1*1000) && $10>T2*1000) 
                    priori=1 ;
                else
                    priori=0 ;
                printf("%d\t%f\t%f\t%f\n",priori,$8/1000,$9/1000,$10/1000);
            }
    }'   trio_density.data.txt >trio.4r.matrix || exit 1
    date >>'11.step_2_done'
else
    echo "skip extract trio.4r.matrix  paternal due to 11.step_2_done exist"
fi

echo "__END__"
