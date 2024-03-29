#!/bin/bash
ABSPATH=`realpath $0`
PIPEDIR=`dirname $ABSPATH`
#################################################################
# required settings 
#################################################################
PATNGS=''
MATNGS=''
SONTGS=''
MAINDIR=`realpath $PIPEDIR/../`
###############################################################################
# parse arguments
###############################################################################
function usage(){
echo """
usage of example pipeline - sysc_strobmer_mode.sh:
    sysc_strobmer_mode.sh <--maternal mat.ngs> \\
                          <--paternal pat.ngs> \\
                          <--offspring son.tgs.fasta> \\
                          [--sspath path-to-SybointScreener]
"""
}

echo "CMD :$0 $*"
if [[ $# == 0 || $1 == '-h' || $1 == "--help" ]] ; then
    usage
    exit 0
fi
while [[ $# > 0 ]]
do
    case $1 in
        "--maternal")
            MATNGS=$2
            shift
            ;;
        "--paternal")
            PATNGS=$2
            shift
            ;;
        "--offspring")
            SONTGS=$2
            shift
            ;;
        "--sspath")
            MAINDIR=$2
            shift
            ;;
        *)
            echo "invalid params : \"$1\" . exit ... "
            exit
        ;;
    esac
    shift
done

SEQKIT=$MAINDIR/thirdparty/seqkitbin/seqkit
SYSC=$MAINDIR/sysc
# sanity check
if [[ ! -e $SEQKIT || ! -e $SYSC ]] ; then
    echo "ERROR : please provide valid sspath."
    exit 1
fi
# sanity check
if [[ $MATNGS == '' || $PATNGS == "" || $SONTGS == "" ]] ; then
    echo "ERROR : please provide all three input files."
    exit 1
fi
# sanity check
if [[ ! -e $MATNGS  || ! -e $PATNGS  || ! -e $SONTGS ]] ; then
    echo "ERROR : please provide valid input files."
    exit 1
fi

#################################################################
# the strobmer s40 mode with BGMM clustering, feel free to edit
#################################################################

$SYSC build_s40   --thread 8 \
                  --paternal $PATNGS \
                  --maternal $MATNGS \
                  --auto_bounds 1 \
                  1>build_s40.log  2>build_s40.err

$SYSC density_s40 --thread 8 \
                  --offspring $SONTGS \
                  1>density_s40.log 2>density_s40.err 

$SYSC trio_result_s40 1>trio_result_s40.log 2>trio_result_s40.err

#################################################################
# extract reads basd on sysc result
#################################################################
$SEQKIT grep -f step02.2.s40/trio_only.filtered_readname.txt \
                $SONTGS >host.fa

$SEQKIT grep -v -f step02.2.s40/trio_only.filtered_readname.txt \
                $SONTGS > symbiont.fa
