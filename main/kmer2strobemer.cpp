#include "strobemer/strobemer.h"
#include <getopt.h>
#include <iostream>
#include <string>
#include <cassert>

void printUsage() {
    std::cerr<<"Usage :  kmer2strobemer <-n/--nkmer n> <-k/--ksize k> <-w/--wsize w>  [-m|--minstrobe(default randstrobe)]  0<input.kmers 1>output.stronemers 2>log.txt \n";
}

int main(int argc ,char ** argv) {
    int n , k ,w;
    strobemer_type type = strobemer_type::randstrobe;
    static struct option long_options[] = {
        {"nkmer",required_argument, NULL, 'n'},
        {"ksize",required_argument, NULL, 'k'},
        {"wsize",required_argument,NULL, 'w'},
        {"minstrobe",required_argument,NULL, 'm'},
        {"help",  no_argument,       NULL, 'h'},
        {0, 0, 0, 0}
    };
    static char optstring[] = "n:k:w:mh";
    while(1) {
        int c = getopt_long(argc, argv, optstring, long_options, NULL);
        if (c<0) break;
        switch (c){
            case 'n':
                n=atoi(optarg);
                break;
            case 'k':
                k=atoi(optarg);
                break;
            case 'w':
                w=atoi(optarg);
                break;
            case 'm':
                type = strobemer_type::minstrobe;
                break;
            case 'h':
                printUsage();
                return 0;
            default :
                printUsage();
                return -1;
        }
    }
    if( n<2 || k<2 || w<=k || type  == strobemer_type::unknow ) {
        printUsage();
        return -1 ;
    }
    strobemer::init(n,k,w,type);
    std::string line;
    char rc_line[1000]; // support max-kmer-size 1000 bp
    long long line_num = 1 ;
    while(!std::getline(std::cin,line).eof()){
        assert(line.size() == strobemer::strobmer_span());
        reverse_complete(line.c_str(),line.size(),rc_line);
        strobemer buff;
        strobemer::chop_strobemer(line.c_str(),line.size(),&buff);
        if( buff.valid) {
            std::cout<<buff.to_string()<<'\n';
        } else {
            std::cerr<<"N found in line "<<line_num<<" kmer="<<line<<'\n';
        }
        strobemer::chop_strobemer(rc_line,line.size(),&buff);
        if( buff.valid) {
            std::cout<<buff.to_string()<<'\n';
        } else {
            std::cerr<<"N found in line "<<line_num<<" kmer="<<line<<'\n';
        }
    }
    return 0;
}
