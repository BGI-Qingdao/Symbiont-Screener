#!/bin/bash

if [[ $1 == "-h" || $1 == "--help" ]] ; then 
    echo "Usage : $0 final.result.txt reads.fa"
    exit 
fi

awk -f '>| ' '{
    if( FILENAME == ARGV[1] ) {
        if( $3 == 1 ) t[$1]=1;
    }else {
        if( NF > 1 ) {
            if( $2 in t ) pass = 1;
            else pass = 0;
        }
        if( pass == 1 ) print $0;
    }
}' $1 $2
