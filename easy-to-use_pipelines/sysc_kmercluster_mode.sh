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
    sysc_kmercluster_mode.sh     <--maternal mat.ngs> \\
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
# the strobmer k21 mode with BGMM clustering, feel free to edit
#################################################################

$SYSC build_k21   --thread 8 \
                  --paternal $PATNGS \
                  --maternal $MATNGS \
                  --auto_bounds 1 \
                  1>build_k21.log  2>build_k21.err

$SYSC density_k21 --thread 8 \
                  --offspring $SONTGS \
                  1>density_k21.log 2>density_k21.err 

$SYSC trio_result_k21 1>trio_result_k21.log 2>trio_result_k21.err

$SYSC cluster_k21 --thread 8 \
                  --offspring $SONTGS \
                  1>cluster_k21.log  2>cluster_k21.err

$SYSC consensus_cluster_k21 --min_hit 10 \
                  1>consensus_cluster_k21.log \
                  2>consensus_cluster_k21.err

#################################################################
# extract reads basd on sysc result
#################################################################
$SEQKIT grep -f step03.2.k21/readname.min-hit10.txt \
                $SONTGS >host.fa

$SEQKIT grep -f step03.2.k21/readname.other.min-hit10.txt \
                $SONTGS > symbiont.fa
