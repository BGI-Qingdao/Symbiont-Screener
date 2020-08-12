#!/bin/bash

###############################################################################
# usage function
###############################################################################
function usage(){
echo """
Usage    :
    ./ContamFilter.sh [OPTION]

Filter contamination from filial tgs reads by parental kmer sets.

Options  :
        --filial            filial TGS reads file in FASTA/Q format.
                            file in gzip format can be accepted, but filename must end by ".gz".

        --format            fasta/fastq . set the format of --filial.
                            [ optional, default fasta. ]

        --thread            thread num.
                            [ optional, default 8 threads. ]

        --memory            x (GB) of memory to used by meryl.
                            [ optional, default 50GB. ]

        --mer               mer-size
                            [ optional, default 21. ]

        --use_existing_lib  on/off
                            [ optional, default off]

        --paternal          paternal NGS reads file in FASTQ format.
                            ( note : gzip format is NOT supported. )
                            [ optional, needed when --use_existing_libs off ]

        --maternal          maternal NGS reads file in FASTQ format.
                            ( note : gzip format is NOT supported. )
                            [ optional, needed when --use_existing_libs off ]

        --paternal_mer      existing paternal specific kmer lib file.
                            [ optional, needed when --use_existing_libs on ]

        --maternal_mer      existing maternal specific kmer lib file.
                            [ optional, needed when --use_existing_libs on ]

        --common_mer        existing parental common kmer lib file.
                            [ optional, needed when --use_existing_libs on ]

        --density_common    the minimum density of common_kmer
                            [ optional, default 0.01]

        --density_specific  the minimum density of parental specific kmer
                            [ optional, default 0.00005]

        --count_specific    the minimum count of parental specific kmer
                            [ optional, default 2]

        --help              print this usage message.

Examples :

    ./ContamFilter.sh --filial son.fasta  --paternal father.fastq --maternal mater.fastq 

    # if the filial read file follow fastq format :
    ./ContamFilter.sh --paternal father.fastq --maternal mater.fastq --filial son.fastq --format fastq

    # if there are more than one filial read files :
    ./ContamFilter.sh --paternal father.fastq --maternal mater.fastq --filial son.L01.fasta --filial son.L02.fasta

    # use existing libs
    ./ContamFilter.sh --filial son.fasta  --use_existing_lib on  --paternal_mer p.kmer\\
                                                                 --maternal_mer m.kmer\\
                                                                 --common_mer   c.kmer
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
FILIAL=""
PATERNAL_MER=''
MATERNAL_MER=''
COMMON_MER=''
DENSITY_COMMON=0.01
DENSITY_SPECIFIC=0.00005
COUNT_SPECIFIC=2
FORMAT='fasta'
USE_EXISTING_LIB='off'
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
            PATERNAL=$2
            shift
            ;;
        "--maternal")
            MATERNAL=$2
            shift
            ;;
        "--filial")
            FILIAL=$2" "$FILIAL
            shift 
            ;;
        "--format")
            FORMAT=$2
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
        "--common_mer")
            COMMON_MER=$2
            shift
            ;;
        "--density_common")
            DENSITY_COMMON=$2
            shift
            ;;
        "--density_specific")
            DENSITY_SPECIFIC=$2
            shift
            ;;
        "--count_specific")
            COUNT_SPECIFIC=$2
            shift
            ;;
        "--use_existing_lib")
            USE_EXISTING_LIB=$2
            if [[ $USE_EXISTING_LIB != "on" && $USE_EXISTING_LIB != "off" ]] ; then 
                echo "--use_existing_lib must be on or off but this is $USE_EXISTING_LIB ! exit ..."
                exit 1
            fi
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
echo "    filial input    : $FILIAL"
echo "    filial format   : $FORMAT"
echo "    memory          : $MEMORY GB"
echo "    thread          : $CPU "
echo "    mer             : $MER "
echo "    use_existing_lib: $USE_EXISTING_LIB"
echo "    paternal input  : $PATERNAL"
echo "    maternal input  : $MATERNAL"
echo "    patternal_mer   : $PATERNAL_MER"
echo "    matternal_mer   : $MATERNAL_MER"
echo "    common_mer      : $COMMON_MER"
echo "ContamFilter.sh in  : $SPATH"

MERYL=$SPATH"/bin/meryl"
CLASSIFY=$SPATH"/classify_3lib"
ANALYSIS=$SPATH"/analysis_kmercount.sh"

# sanity check
if [[ $MEMORY -lt 1  || $CPU -lt 1 || -z $FILIAL || $MER -lt 11 ]]
    echo "ERROR : arguments invalid ... exit!!! "
    exit 1
fi

if [[ $FORMAT != 'fasta' && $FORMAT != 'fastq' ]] ; then 
    echo "ERROR : format invalid ... exit!!!"
    exit 1
fi

if [[ $USE_EXISTING_LIB == 'on' ]] ; then 
    if [[ -z $PATERNAL_MER || -z $MATERNAL_MER || -z $COMMON_MER ]] ; then 
        echo "ERROR : --paternal_mer && --maternal_mer && --commmon_mer are needed when --use_existing_lib on"
        echo "exit ..."
        exit 1
    fi
else 
    if [[ -z $PATERNAL || -z $MATERNAL  ]] ; then 
        echo "ERROR : --paternal && --maternal are needed when --use_existing_lib off"
        echo "exit ..."
        exit 1
    fi
fi

if [[ ! -e $CLASSIFY ]] ; then 
    echo "ERROR : please run \"make\" command in $SPATH before using this script! exit..."
    exit 1
fi

for x in $FILIAL
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
if [[ $USE_EXISTING_LIB == "off" ]] ; then 
    if [[ ! -e 'step_1_done' ]]  ; then
        $MERYL count k=$MER output maternal.meryl threads=$CPU memory=$MEMORY $MATERNAL || exit 1
        $MERYL histogram maternal.meryl >maternal.hist || exit 1
        awk -f $SPATH"/find_bounds.awk" maternal.hist > maternal.bounds.txt || exit 1
        date >>'step_1_done'
    else
        echo "skip meryl count maternal due to step_1_done exist"
    fi

    if [[ ! -e 'step_2_done' ]]  ; then
        $MERYL count k=$MER output paternal.meryl threads=$CPU memory=$MEMORY $PATERNAL || exit 1
        $MERYL histogram paternal.meryl >paternal.hist  || exit 1
        awk -f $SPATH"/find_bounds.awk" paternal.hist > paternal.bounds.txt || exit 1
        date >>'step_2_done'
    else
        echo "skip meryl count paternal due to step_2_done exist"
    fi

    MLOWER=`grep LOWER_INDEX maternal.bounds.txt| awk -F '=' '{print $2}'`
    MUPPER=`grep UPPER_INDEX maternal.bounds.txt| awk -F '=' '{print $2}'`
    PLOWER=`grep LOWER_INDEX paternal.bounds.txt| awk -F '=' '{print $2}'`
    PUPPER=`grep UPPER_INDEX paternal.bounds.txt| awk -F '=' '{print $2}'`
    echo "  the real used kmer-count bounds of maternal is [ $MLOWER , $MUPPER ] "
    echo "  the real used kmer-count bounds of paternal is [ $PLOWER , $PUPPER ] "

    if [[ ! -e 'step_3_done' ]]  ; then
        $MERYL difference output mat_only.meryl threads=$CPU memory=$MEMORY maternal.meryl paternal.meryl
        MLOWER=$(($MLOWER-1))
        MUPPER=$(($MUPPER+1))
        $MERYL greater-than $MLOWER  mat_only.meryl output 'mat_only.gt'$MLOWER'.meryl' || exit  1
        $MERYL less-than $MUPPER  'mat_only.gt'$MLOWER'.meryl' output 'mat_only.gt'$MLOWER'.lt'$MUPPER'.meryl' || exit 1
        ln -s  mat_only.gt'$MLOWER'.lt'$MUPPER'.meryl mat_only.filtered.meryl
        $MERYL print mat_only.filtered.meryl | awk '{print $1}' >maternal.mer
        date >>'step_3_done'
    else
        echo "skip get maternal.mer due to step_3_done exist"
    fi


    if [[ ! -e 'step_4_done' ]]  ; then
        $MERYL difference output pat_only.meryl threads=$CPU memory=$MEMORY paternal.meryl maternal.meryl || exit 1
        PLOWER=$(($PLOWER-1))
        PUPPER=$(($PUPPER+1))
        $MERYL greater-than $PLOWER  pat_only.meryl output 'pat_only.gt'$PLOWER'.meryl' || exit  1
        $MERYL less-than $PUPPER  'pat_only.gt'$PLOWER'.meryl' output 'pat_only.gt'$PLOWER'.lt'$PUPPER'.meryl' || exit 1

        ln -s  pat_only.gt'$PLOWER'.lt'$PUPPER'.meryl pat_only.filtered.meryl
        $MERYL print pat_only.filtered.meryl | awk '{print $1}' >paternal.mer
        date >>'step_4_done'
    else
        echo "skip get paternal.mer due to step_4_done exist"
    fi

    if [[ ! -e 'step_5_done' ]]  ; then
        $MERYL intersect-sum output common.meryl threads=$CPU memory=$MEMORY paternal.meryl maternal.meryl  || exit 1
        CLOW=$((($PLOWER+$MLOWER)/2-1))
        $MERYL greater-than $CLOW output 'common.gt'$CLOW'.meryl'
        ln -s 'common.gt'$CLOW'.meryl' 'common.filtered.meryl'
        $MERYL print 'common.filtered.meryl'  | awk '{print $1}' >common.mer
        date >>'step_5_done'
    else
        echo "skip get common.mer due to step_5_done exist"
    fi
fi

###############################################################################
# extract reads
###############################################################################
for x in $FILIAL
do 
    READ="$READ"" --read ""$x"
done

if [[ USE_EXISTING_LIB == "off" ]] ; then
    PATERNAL_MER='paternal.mer'
    MATERNAL_MER='maternal.mer'
    COMMON_MER='common.mer'
fi

if [[ ! -e "step_6_done" ]] ; then
    $CLASSIFY --hap $PATERNAL_MER --hap $MATERNAL_MER --hap $COMMON_MER \
        --thread $CPU  $READ --format $FORMAT >phasing.out  2>phasing.log || exit 1
    date >> "step_6_done"
else
    echo "skip get phasing.out due to step_6_done exist"
fi

if [[ ! -e "step_7_done" ]] ; then
    grep  haplotype phasing.out | awk DENSITY_COMMON=$DENSITY_COMMON DENSITY_SPECIFIC=$DENSITY_SPECIFIC COUNT_SPECIFIC=$COUNT_SPECIFIC  '{if($10>=DENSITY_COMMON && ( $8+$9>DENSITY_SPECIFIC) && ($11+$12>COUNT_SPECIFIC))  print $2}' >filtered.readname.txt
    date >> "step_7_done"
else
    echo "skip get filtered.readname.txt due to step_7_done exist"
fi

if [[ ! -e "step_8_done" ]] ; then
    # extract the reads
    if [[ $FORMAT == 'fasta' ]] ; then
        for x in $FILIAL
        do
            name=`basename $x`
            if [[ ${name: -3} == ".gz" ]] ; then
                name=${name%%.gz}
                gzip -dc $x | awk  -F '>| '  ' {if( FILENAME == ARGV[1] ) { s[$1]=1} else { if( NF>1){ if ($2 in s ){ print $0 ; c=1;} else {c=0} } else { if(c==1) { print $0 ; }  } } }' filtered.readname.txt  - >$name".filtered.fasta"
            else 
                awk  -F '>| '  ' {if( FILENAME == ARGV[1] ) { s[$1]=1} else { if( NF>1){ if ($2 in s ){ print $0 ; c=1;} else {c=0} } else { if(c==1) { print $0 ; }  } } }' filtered.readname.txt  $x >$name".filered.fasta"
            fi
        done
    else
        for x in $FILIAL
        do
            name=`basename $x`
            if [[ ${name: -3} == ".gz" ]] ; then
                name=${name%%.gz}
                gzip -dc $x | awk  -F '>|@| '  ' {if( FILENAME == ARGV[1] ) { s[$1]=1} else { if(FNR %4==1 && NF>1){ if ($2 in s ){ print $0 ; c=1;} else {c=0} } else { if(c==1) { print $0 ; }  } } }' filtered.readname.txt  - >$name".filtered.fasta"
            else 
                awk  -F '>|@| '  ' {if( FILENAME == ARGV[1] ) { s[$1]=1} else { if(FNR %4==1 && NF>1){ if ($2 in s ){ print $0 ; c=1;} else {c=0} } else { if(c==1) { print $0 ; }  } } }' filtered.readname.txt  $x >$name".filered.fasta"
            fi
        done
    fi
    date >> "step_8_done"
else
    echo "skip get reads due to step_8_done exist."
fi
echo "phase reads done"
date
echo "__END__"
