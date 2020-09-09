#include <iostream>
#include <map>
#include <unordered_set>
#include <fstream>
#include <sstream>
#include <cassert>
#include <ctime>
#include <thread>
#include <stack>
#include <mutex>
#include <vector>
#include <chrono>
#include <getopt.h>
#include "gzstream/gzstream.h"
#include <cmath>
#include <set>

int g_nmer=4;

inline bool isA(char c) { return c == 'A' || c == 'a' ; }
inline bool isG(char c) { return c == 'G' || c == 'g' ; }
inline bool isC(char c) { return c == 'C' || c == 'c' ; }
inline bool isT(char c) { return c == 'T' || c == 't' ; }

void logtime(){
    time_t now = time(0);
    char* dt = ctime(&now);
    std::cerr<<dt<<std::endl;
}
//
// kmer relate functions
//
std::map<char,char> g_oppo ;
void InitMap(){
    g_oppo['a']='T';
    g_oppo['A']='T';
    g_oppo['g']='C';
    g_oppo['G']='C';
    g_oppo['c']='G';
    g_oppo['C']='G';
    g_oppo['t']='A';
    g_oppo['T']='A';
    g_oppo['n']='N';
    g_oppo['N']='N';
}

std::string reverse_complement(const std::string & kmer)
{
    std::string ret = kmer;
    for(int i = 0 ; i< (int)kmer.size(); i++)
        ret[kmer.size()-i-1] = g_oppo[kmer[i]];
    return ret ;
}

std::set<std::string> GetAllMer(int n ) 
{
    std::set<std::string> ret ;
    if( n == 1 ) {
        ret.insert("A");
        ret.insert("G");
        ret.insert("C");
        ret.insert("T");
        return ret ;
    } else {
        auto prev_set = GetAllMer(n-1);
        for( const auto & x : prev_set ) {
            ret.insert(x+"A");
            ret.insert(x+"G");
            ret.insert(x+"C");
            ret.insert(x+"T");
        }
        return ret ;
    }
}

struct NmerFreq
{
    std::map<std::string,int> cache;

    void DealSeq(const std::string & a) 
    {
        std::string seq = a;

        for(int i = 0 ; i < (int)seq.size() ; i++)
        {
            seq[i] = std::toupper(seq[i]);
        }

        for( int i = 0 ; i <(int)seq.size()-n+1 ;i++ )
        {
            std::string nmer=seq.substr(i,n);
            auto itr = cache.find(nmer) ;
            if( itr == cache.end() ) cache[nmer]=1;
            else itr->second ++ ;
        }

        seq = reverse_complement(seq);

        for( int i = 0 ; i <(int)seq.size()-n+1 ;i++ )
        {
            std::string nmer=seq.substr(i,n);
            auto itr = cache.find(nmer) ;
            if( itr == cache.end() ) cache[nmer]=1;
            else itr->second ++ ;
        }
    }
    static int n ;
    static std::set<std::string> allmer;

    static void GenAll(int n) {
        NmerFreq::n = n ;
        allmer = GetAllMer(n);
    }

    void Print() const{
        int t=0;
        for ( const auto & p: cache )
            t+=p.second;
        for( const auto & tmp : allmer )
            if ( cache.find(tmp) != cache.end() )
                std::cout<<'\t'<< float(cache.at(tmp))/float(t);
            else
                std::cout<<"\t0";
    }

};

int NmerFreq::n ;
std::set<std::string> NmerFreq::allmer ;

float GC_content(const std::string & read) {
    int gc =0 ; int total=0 ;
    for( char c : read ){
        total ++;
        if( isG(c) || isC(c) )
            gc ++ ;
    }
    return float(gc)/float(total) ;
}


//
// Output cache
//
struct OutputCache {
    struct ReadInfo {
        long long  read_id;
        std::string name;
        int read_length ;
        float gc_content;
        NmerFreq nmer_freq;
    };

    std::mutex mm;
    std::map<int , ReadInfo> output_cache;

    void SaveOutput(const ReadInfo & data){
        mm.lock();
        output_cache[data.read_id] = data;
        mm.unlock();
    }

    void PrintOutput() const {
        for( const auto & pair : output_cache) {
            const auto & data = pair.second;
            std::cout<<data.name;
            std::cout<<'\t'<<data.gc_content;
            data.nmer_freq.Print();
            std::cout<<'\n';
        }
    }
};

//
//reads relate functions
//

struct ReadBasic{
    long long id ;
    std::string head ;
    std::string seq ;
};
struct MultiThread {
    int t_nums ;
    bool end;
    bool busy;
    OutputCache data;

    std::vector< std::stack<ReadBasic> >  caches;
    std::mutex * locks;
    std::thread ** threads; 

    MultiThread(int t_num){
        t_nums = t_num ;
        busy = false;
        end=false;

        locks = new std::mutex[t_num];
        threads = new std::thread*[t_num];
        for(int i = 0 ; i< t_num ; i++){
            caches.push_back(std::stack< ReadBasic >());
            threads[i] = new std::thread([this,i](){int index=i ;Worker(index); });
        }
    }
    void wait(){
        for(int i = 0 ; i <t_nums; i++){
            threads[i]->join();
            delete threads[i];
        }
    }
    ~MultiThread(){
        delete [] threads;
        delete [] locks;
    }

    void Worker(int index){
        ReadBasic job;
        while(true){
            locks[index].lock();
            if( caches[index].empty() ){
                busy = false ;
                locks[index].unlock();
                if(end) return ;
                std::this_thread::sleep_for(std::chrono::microseconds(10));
                continue;
            }
            if( ! caches[index].empty() ){
                std::swap(job,caches[index].top());
                if( caches[index].size() > 10000 ) busy = true ;
                else if ( caches[index].size() < 3000 ) busy = false ;
                caches[index].pop();
                locks[index].unlock();
                process_reads(job.head,job.seq,job.id);
            } else 
                locks[index].unlock();
        }
    }

    void process_reads(const std::string & head ,
            const std::string & seq , long long  readid) {
        OutputCache::ReadInfo tmp ;
        tmp.read_id = readid ;
        tmp.name = head.substr(1) ;
        tmp.read_length = seq.size();
        tmp.gc_content = GC_content(seq);
        tmp.nmer_freq.DealSeq(seq);
        data.SaveOutput(tmp);
    }
    void submit(std::string & head ,std::string & seq, long long id){
        //while( busy ) { std::this_thread::sleep_for(std::chrono::seconds(1));}
        int tid = id % t_nums;
        ReadBasic tmp ;
        std::swap(tmp.head,head);
        std::swap(tmp.seq,seq);
        tmp.id = id ;
        locks[tid].lock();
        caches[tid].push(std::move(tmp));
        locks[tid].unlock();
    }
};

std::istream * openfile( const std::string & file ){
    std::istream *in ;
    bool gz_file = false;
    if( file.size() > 3 ) {
        int end=file.size() ;
        if( file[end-3] == '.' && file[end-2] == 'g' && file[end-1]=='z' ) {
            gz_file = true ;
        }
    }
    if ( gz_file )
        in = new igzstream(file.c_str());
    else 
        in = new std::ifstream(file);
    return in ;
}

void processFastq(const std::string & file,int t_num){
    std::string head;
    std::string seq;
    std::string tmp;
    MultiThread mt(t_num);
    std::istream *in = openfile(file);
    long long id = 0 ;
    while(!std::getline(*in,head).eof()){
        if( head[0] == '>' ) {
            std::cerr<<"fasta detected . ERROR . please use \"--format fasta\". exit ... "<<std::endl;
            exit(1);
        }
        std::getline(*in,seq);
        std::getline(*in,tmp);
        std::getline(*in,tmp);
        mt.submit(head, seq, id );
        id ++ ;
    }
    mt.end=true;
    delete in ;
    mt.wait();
    mt.data.PrintOutput();
}

void processFasta(const std::string & file,int t_num){
    std::string head;
    std::string seq;
    std::string tmp;
    MultiThread mt(t_num);
    std::istream *in = openfile(file);
    long long id = 0 ;
    while(!std::getline(*in,tmp).eof()){
        if( tmp.empty() ) continue ;
        if( tmp[0] == '@' || tmp[0] == '+' ) {
            std::cerr<<"fasta detected . ERROR . please use \"--format fastq\". exit ... "<<std::endl;
            exit(1);
        }
        if( tmp[0] == '>' ) {
            if( id > 0 ) {
                  mt.submit(head,seq,id );
            }
            std::swap(head,tmp);
            seq = "";
            id ++ ;
        }
        else{
            seq += tmp ;
        }
    }
    mt.submit(head, seq, id );
    mt.end=true;
    delete in ;
    mt.wait();
    mt.data.PrintOutput();
}

void printUsage(){
    std::cerr<<"Uasge :\n\tgc_nmer  --read read.fa [--kmer k(default 4)] [--read read_2.fa] [--thread t_num (8 default) ] [--format fasta/fastq (default fasta)]"<<std::endl;
    std::cerr<<"notice : --read accept file in gzip format , but file must end by \".gz\""<<std::endl;
    std::cerr<<"warn   : --read default only accept fasta read."<<std::endl;
    std::cerr<<"         add --format fastq if --read refer to fastq file."<<std::endl;
}

//
// Main function
//
int main(int argc ,char ** argv ) {
    InitMap();
    static struct option long_options[] = {
        {"read",  required_argument, NULL, 'r'},
        {"format",required_argument, NULL, 'f'},
        {"thread",required_argument, NULL, 't'},
        {"kmer",required_argument,NULL, 'k'},
        {"help",  no_argument,       NULL, 'h'},
        {0, 0, 0, 0}
    };
    static char optstring[] = "r:t:f:k:h";
    std::vector<std::string> read;
    std::string format="fasta";
    int t_num=8;
    while(1){
        int c = getopt_long(argc, argv, optstring, long_options, NULL);
        if (c<0) break;
        switch (c){
            case 'r':
                read.push_back(std::string(optarg));
                break;
            case 'f':
                format = std::string(optarg);
                break;
            case 't':
                t_num = atoi(optarg);
                break;
            case 'k':
                g_nmer = atoi(optarg);
                break;
            case 'h':
            default :
                printUsage();
                return -1;
        }
    }
    if( g_nmer<1 || read.empty()|| t_num< 1) {
        printUsage();
        return -1;
    }
    std::cerr<<" K = "<<g_nmer<<" . the greater K you use, the more memory will be used!"<<std::endl;
    if( format != "fasta" && format != "fastq" ){
        std::cerr<<" ERROR : invalid format : ["<<format<<"] . exit ...\n";
        return -1;
    }
    logtime();
    std::cerr<<"__START__"<<std::endl;
    NmerFreq::GenAll(g_nmer);
    for(const auto r : read ){
        std::cerr<<"__process read: "<<r<<std::endl;
        logtime();
        if( format == "fasta")
            processFasta(r,t_num);
        else 
            processFastq(r,t_num);
        std::cerr<<"__process read done__"<<std::endl;
        logtime();
    }
    std::cerr<<"__END__"<<std::endl;
    logtime();
}
