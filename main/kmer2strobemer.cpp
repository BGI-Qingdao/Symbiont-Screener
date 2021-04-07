#include "strobemer/strobemer.h"
#include <getopt.h>
#include <iostream>
#include <string>
#include <cassert>

void printUsage() {
    std::cerr<<"Usage :  kmer2strobemer <-n/--nkmer n> <-k/--ksize k> <-w/--wsize w>  <-m|-r/--minstrobe|--randstrobe>  0<input.kmers 1>output.stronemers 2>log.txt \n";
}

int main(int argc ,char ** argv) {
    int n , k ,w;
    strobemer_type type = strobemer_type::unknow;
    static struct option long_options[] = {
        {"nkmer",required_argument, NULL, 'n'},
        {"ksize",required_argument, NULL, 'k'},
        {"wsize",required_argument,NULL, 'w'},
        {"minstrobe",required_argument,NULL, 'm'},
        {"randstrobe",required_argument,NULL, 'r'},
        {"help",  no_argument,       NULL, 'h'},
        {0, 0, 0, 0}
    };
    static char optstring[] = "n:k:w:mrh";
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
            case 'r':
                type = strobemer_type::randstrobe;
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
    long long line_num = 1 ;
    while(!std::getline(std::cin,line).eof()){
        assert(line.size() == strobemer::strobmer_span());
        strobemer buff;
        strobemer::chop_strobemer(line.c_str(),line.size(),&buff);
        if( buff.valid) {
            std::cout<<buff.to_string()<<'\n';
        } else {
            std::cerr<<"N found in line "<<line_num<<" kmer="<<line<<'\n';
        }
    }
    return 0;
}
