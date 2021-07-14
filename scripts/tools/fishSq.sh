#!/bin/bash

SEQKIT=seqkit

if [[ ! -e $SEQKIT ]] ; then
    echo "please assign SEQKIT path first by edit this script."
    echo "try conda install seqkit"
    exit 1
fi

function usage(){
    echo """
    Usage :
    fish_target <-f final.result.txt> <-t threshold> <-i input.fa> <-o output.fa>

    Brief :
    filter the target sequences ( \"hit-count\"+\"second-hit-count\">\"threshold\" )
    """
}

###############################################################################
# parse arguments
###############################################################################
Fresult=""
Ifa=""
Ofa=""
T=1
echo "CMD :$0 $*"
if [[ $# != 8 ]] ; then
    usage
    exit 1
fi
while [[ $# > 0 ]]
do
    case $1 in
        "-f")
            Fresult=$2
            shift
            ;;
        "-i")
            Ifa=$2
            shift
            ;;
        "-o")
            Ofa=$2
            shift
            ;;
        "-t")
            T=$2
            shift
            ;;
        *)
            echo "invalid params : \"$1\" . exit ... "
            exit 0
            ;;
    esac
    shift
done

if [[ ! -e  $Fresult || ! -e $Ifa || $T -lt 1 ]] ; then
    echo "invalid parameters !!!  exit ..."
    usage
    exit 1
fi
if [[ -e $Ofa ]] ; then
    echo "Warn : will overwrite $Ofa now ..."
fi
###############################################################################
# logic
###############################################################################
rName="target.gt"$T".name.txt"
awk -v t=$T '{if((int($4)+int($5))>t) print $1 }' $Fresult > $rName
$SEQKIT grep -f $rName <$Ifa >$Ofa
echo "Done ..."
