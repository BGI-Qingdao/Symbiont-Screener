#!/bin/bash

if [[ $1 == "-h" || $1 == "--help" ]] ; then 
    echo "Usage : $0 final.result.txt reads.fq"
    exit 
fi

awk -f '#|/' '{
    if( FILENAME == ARGV[1] ) {
        if( $3 == 1 ) s[$1]=1;
    } else { 
        if( FNR%4==1 && NF>1 ){
            total = total + 1 ;
            if ( $2 in s ){ 
                print $0 ;
                used = used + 1;
                c=1;
            } else {
                c=0
            }
        } else {
            if(c==1) { 
                print $0 ; 
            }
        }
    }
}' $1 $2

