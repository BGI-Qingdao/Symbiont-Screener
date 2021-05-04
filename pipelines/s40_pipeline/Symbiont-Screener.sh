#!/bin/bash

function usage(){
echo """
Usage   :
    ./Symbiont-Screener.sh [options]

Options :

   Basic parameters:

        --paternal          contaminated paternal NGS read file in FASTQ format.

        --maternal          contaminated maternal NGS read file in FASTQ format.

        --offspring         contaminated offspring's TGS read file.
                            gzip format file is supported but should end by '.gz'

        --offspring_format  fasta/fastq (default fasta)

        --thread            thread num.
                            [ optional, default 8 threads. ]

  For marker generation:
  
        --size              initial hash size for jellyfish 
                            [ optional, default 10GB. ]
			    
        --low_depth          estimated lower depth for k-mer histogram (default 0)

        --high_depth        estimated higher depth for k-mer histogram (default 0)
                            this pipeline will automatically choose lower and higher depth threasholds when both --low_depth and --high_depth are not set.
                            if the user estimates that sequencing coverage or depth of the host is around x , then please set low_depth=x/4 and high_depth=x*[3 or 5]

  For trio-binning-based detection:

        --threshold1        minimum of parental-specific kmer density (default 0.001)
                            for ONT reads(error rate~=15%), we recommand 0.001.
                            for PacBio reads(error rate<5%), we recommand 0.002-0.005.

        --threshold2        minimum of shared kmer density (default 0.1)

  For BGMM-clustering-based detection:

        --cluster           (1/0) use clustering or not. default(0)

        --shortest          length threshold for clustering ( default 5000 )
                            only reads with lenghs > the shortest will be used for clustering.
                            shorter reads ( <=5k ) often create noisy points and hamper accurate clustering.

        --loop              number of BGMM clustering for consensus (default 10)

        --python3           python3 path ( default python3 )

        --seed              random seed ( default 42 )

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
KMER=40
NMER=3
CPU='8'
SIZE='10G'
L_SHORTEST=5000
LOOP=10
RSEED=42
CLUSTER=0
L_DEPTH=0
H_DEPTH=0
SPATH=`dirname $0`
STEP0=$SPATH/00.BuildTrioLib/build_trio_kmer.sh
K2S=$SPATH/main/kmer2strobemer
STEP1=$SPATH/11.GetTrioMatrix_tgs.sh
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
        "--size")
            SIZE=$2
            shift
            ;;
        "--seed")
            RSEED=$2
            shift
            ;;
        "--thread")
            CPU=$2
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
        "--high_depth")
            H_DEPTH=$2
            shift
            ;;
	"--shortest")
            L_SHORTEST=$2
            shift
            ;;
        "--low_depth")
            L_DEPTH=$2
            shift
            ;;
        "--offspring_format")
            OFFSPRING_FORMAT=$2
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
        "--threshold1")
            THRESHOLD1=$2
            shift
            ;;
        "--threshold2")
            THRESHOLD2=$2
            shift
            ;;
	 "--cluster")
            CLUSTER=$2
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
if [[ $L_DEPTH == 0 && $H_DEPTH == 0 ]] ; then
       $STEP0 --paternal "$PATERNAL" \
               --maternal "$MATERNAL" \
               --mer $KMER \
               --thread $CPU \
               --auto_bounds \
               --size $SIZE  || exit 1
else
       if [[  $L_DEPTH -lt 1 || $H_DEPTH -lt 1 || $L_DEPTH -gt $H_DEPTH ]] ; then
               echo "L_DEPTH or H_DEPTH error ... exit !"
               exit 1
       fi
       $STEP0 --paternal "$PATERNAL" \
               --maternal "$MATERNAL" \
               --mer $KMER \
               --thread $CPU \
               --m-lower $L_DEPTH \
               --p-lower $L_DEPTH \
               --m-upper $H_DEPTH \
               --p-upper $H_DEPTH \
               --size $SIZE  || exit 1
fi

if [[ ! -e '00.kmer2strobemer_done' ]]  ; then
    $K2S --nkmer 2 --ksize 10 --wsize 30 <paternal.unique.filter.mer >paternal.strobemer || exit 1
    $K2S --nkmer 2 --ksize 10 --wsize 30 <maternal.unique.filter.mer >maternal.strobemer || exit 1
    $K2S --nkmer 2 --ksize 10 --wsize 30 <common.mer >common.strobemer || exit 1
    awk '{if(FILENAME==ARGV[1]) {m[$1]=1;} else {if($1 in m) m[$1]=0;} }END{for( x in m ) { if(m[x]==1)  print x>"paternal_only.strobemer" ; else print x>"pc.mer";}}' paternal.strobemer maternal.strobemer common.strobemer &
    awk '{if(FILENAME==ARGV[1]) {m[$1]=1;} else {if($1 in m) m[$1]=0;} }END{for( x in m ) { if(m[x]==1)  print x>"maternal_only.strobemer" ; else print x>"mc.mer";}}' maternal.strobemer paternal.strobemer common.strobemer &
    wait

    cat common.strobemer pc.mer mc.mer >common.all.strobemer || exit 1
    awk '{m[$1]=1;}END{for(x in m) print x;}' common.all.strobemer > common_uniq.strobemer ||exit 1
    date >>'00.kmer2strobemer_done'
else
    echo "skip kmer2strobemer due to 00.kmer2strobemer_done exist!"
fi

$STEP1 --paternal_mer paternal_only.strobemer \
    --maternal_mer maternal_only.strobemer \
    --shared_mer common_uniq.strobemer  \
    --offspring "$OFFSPRING" \
    --offspring_format $OFFSPRING_FORMAT  \
    --threshold1 $THRESHOLD1 \
    --threshold2 $THRESHOLD2 \
    --shortest $L_SHORTEST \
    --thread $CPU || exit 1

if [[ $CLUSTER == 1 ]] ; then 
    $STEP2 --nmer $NMER \
           --offspring "$OFFSPRING" \
           --offspring_format $OFFSPRING_FORMAT  \
           --shortest $L_SHORTEST \
           --thread $CPU || exit 1

    if [[ ! -e '30.step_1_done' ]]  ; then
        $PYTHON3  $BGM_MAIN -t trio.4r.matrix  \
                        -m  gc_nmer.matrix \
                        -r  $RSEED \
                        -l $LOOP >cluster_reuslt.txt 2>cluster.log || exit 1
        date >>'30.step_1_done'
    else
        echo "skip cluster by bgm due to 30.step_1_done exist!"
    fi

    if [[ ! -e '30.step_2_done' ]]  ; then
        echo "read_name" >name.txt
        awk -v short=$L_SHORTEST '{if($4>short)print $2;}' trio_density.data.txt >>name.txt
        paste name.txt cluster_reuslt.txt >final.result.txt 
 
        echo "Result in final.result.txt"
        echo "  Format are \"name priori host best-hit-counts second best-hit-counts\""
        echo "NOTICE : reads length not great than $L_SHORTEST had beed excluded from final.result.txt"
        date >>'30.step_2_done'
    else
        echo "skip get final result due to 30.step_2_done exist!"
    fi
    echo "__ALL DONE__"
else
    echo "__DONE__ without cluster"
fi
