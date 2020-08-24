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

std::string reverse_complement(const std::string & kmer){
    std::string ret = kmer;
    for(int i = 0 ; i< (int)kmer.size(); i++)
        ret[kmer.size()-i-1] = g_oppo[kmer[i]];
    return ret ;
}

//
// load & cache maternal unique kmer & paternal unique kmer
//
#define HAPLOTYPES 3
std::unordered_set<std::string> g_kmers[HAPLOTYPES];
int g_K=0;
int g_binsize = 1000;
float g_sd_fac = 2.0f ;
int total_kmers[HAPLOTYPES];
void load_kmers(const std::string & file,int index){
    std::ifstream ifs(file);
    std::string line;
    int total_kmer = 0 ;
    g_kmers[index].reserve(10000000);
    if(index==0){
        std::getline(ifs,line);
        g_K = line.size();
        g_kmers[index].insert(line);
        g_kmers[index].insert(reverse_complement(line));
        total_kmer++;
    }
    while(!std::getline(ifs,line).eof()){
        g_kmers[index].insert(line);
        g_kmers[index].insert(reverse_complement(line));
        total_kmer++;
    }
    total_kmers[index] = total_kmer ;
    std::cerr<<"Recorded "<<total_kmer<<" haplotype "<<index<<" specific "<<g_K<<"-mers\n"; 
}

struct ReadCalc
{
    float gc_fac;

    float aa_ratio; // AA TT
    float at_ratio; // AT
    float ag_ratio; // AG CT
    float ac_ratio; // AC GT

    float gg_ratio; // GG CC
    float ga_ratio; // GA TC
    float gc_ratio; // GC

    float ta_ratio; // TA
    float tg_ratio; // TG CA

    float cg_ratio; // CG
};

inline bool isA(char c) { return c == 'A' || c == 'a' ; }
inline bool isG(char c) { return c == 'G' || c == 'g' ; }
inline bool isC(char c) { return c == 'C' || c == 'c' ; }
inline bool isT(char c) { return c == 'T' || c == 't' ; }

ReadCalc calc_read(const std::string & read ) {
    ReadCalc ret ;
    ret.gc_fac = 0;
    ret.aa_ratio=0;
    ret.at_ratio=0;
    ret.ag_ratio=0;
    ret.ac_ratio=0;
    ret.gg_ratio=0;
    ret.ga_ratio=0;
    ret.gc_ratio=0;
    ret.ta_ratio=0;
    ret.tg_ratio=0;
    ret.cg_ratio=0;
    int total = 0;
    char prev = 0;
    for( char c : read ){
        total ++;
        if( isG(c) || isC(c) )
            ret.gc_fac ++ ;
        if( prev == 0 )  prev = c ;
        else {
            if( isA(prev) && isA(c) ) ret.aa_ratio ++ ;
            if( isA(prev) && isT(c) ) ret.at_ratio ++ ;
            if( isA(prev) && isG(c) ) ret.ag_ratio ++ ;
            if( isA(prev) && isC(c) ) ret.ac_ratio ++ ;

            if( isG(prev) && isA(c) ) ret.ga_ratio ++ ;
            if( isG(prev) && isT(c) ) ret.gc_ratio ++ ;
            if( isG(prev) && isG(c) ) ret.gg_ratio ++ ;
            if( isG(prev) && isC(c) ) ret.ac_ratio ++ ;

            if( isT(prev) && isA(c) ) ret.ta_ratio ++ ;
            if( isT(prev) && isT(c) ) ret.aa_ratio ++ ;
            if( isT(prev) && isG(c) ) ret.tg_ratio ++ ;
            if( isT(prev) && isC(c) ) ret.ga_ratio ++ ;

            if( isC(prev) && isA(c) ) ret.tg_ratio ++ ;
            if( isC(prev) && isT(c) ) ret.ag_ratio ++ ;
            if( isC(prev) && isG(c) ) ret.cg_ratio ++ ;
            if( isC(prev) && isC(c) ) ret.gg_ratio ++ ;

            prev = c ;
        }
    }
    ret.gc_fac  /=total;
    ret.aa_ratio/=total-1;
    ret.at_ratio/=total-1;
    ret.ag_ratio/=total-1;
    ret.ac_ratio/=total-1;
    ret.gg_ratio/=total-1;
    ret.ga_ratio/=total-1;
    ret.gc_ratio/=total-1;
    ret.ta_ratio/=total-1;
    ret.tg_ratio/=total-1;
    ret.cg_ratio/=total-1;

    return ret ;
}

struct MeanCalc
{
    public:
        void Push( float i ){ data.push_back(i);}

        float Mean(float sd_range = 2.0f) {
            if( data.size() == 0 ) return 0 ;
            float total = 0;
            for( float i : data ) total+=i;
            float mean = total / data.size() ;
            float sd = 0 ;
            for( float i : data ) sd += ( i -mean ) * ( i-mean) ;
            sd /= data.size();
            sd = sqrt(sd);
            int count=0;
            total = 0;
            for( float i : data ) {
                if( fabs(i-mean) < sd_range*sd ) {
                    total += i ;
                    count ++ ;
                }
            }
            if( count == 0 ) return 0 ;
            return total / count ;
        }
    private:
        std::vector<float> data;
};
//
// Output cache
//
struct OutputCache {
    struct ReadInfo {
        long long  read_id;
        std::string name;
        int hapCounts[HAPLOTYPES];
        float density[HAPLOTYPES];
        int read_length ;
        ReadCalc rc_data;
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
            double readHapCount = 0 ;
            double secondBest = 0 ;
            std::string readHap = "amibigous";
            double probablity[HAPLOTYPES];
            //double density[HAPLOTYPES];
            for( int i = 0 ; i<HAPLOTYPES ; i++ ) {
                probablity[i] = double(data.hapCounts[i]) / double(g_kmers[i].size());
                //density[i] = double(data.hapCounts[i]) / double( data.read_length -g_K  + 1 );
            }
            for( int i = 0 ; i<HAPLOTYPES ; i++ ) {
                if( probablity[i] > 0 and probablity[i] < readHapCount and probablity[i] > secondBest)
                    secondBest = probablity[i];

                if( probablity[i] > 0 and probablity[i] > readHapCount){
                    readHap = "haplotype"+std::to_string(i);
                    secondBest = readHapCount;
                    readHapCount = probablity[i];
                }
            }
            if( secondBest == 0 and readHapCount != 0 )
                std::cout<<"Read\t"<<data.name<<"\t"<<readHap<<"\t"<<data.read_length;
            else if( readHapCount == 0 and secondBest == 0 )
                std::cout<<"Read\t"<<data.name<<"\t"<<"ambiguous\t"<<data.read_length;
            else if( readHapCount == 0 and secondBest != 0 ) {
                printf("Not possible!\n");
                exit(1);
            }
            else if ( readHapCount / secondBest > 1 )
                std::cout<<"Read\t"<<data.name<<"\t"<<readHap<<"\t"<<data.read_length;
            else
                std::cout<<"Read\t"<<data.name<<"\t"<<"ambiguous\t"<<data.read_length;

            for(int i = 0 ; i<HAPLOTYPES ; i++ ) {
                std::cout<<"\t"<<probablity[i];
            }
            for(int i = 0 ; i<HAPLOTYPES ; i++ ) {
                std::cout<<"\t"<<data.density[i];
            }
            for(int i = 0 ; i<HAPLOTYPES ; i++ ) {
                std::cout<<"\t"<<data.hapCounts[i];
            }

            std::cout<<'\t'<<data.rc_data.gc_fac;
            std::cout<<'\t'<<data.rc_data.aa_ratio;
            std::cout<<'\t'<<data.rc_data.at_ratio;
            std::cout<<'\t'<<data.rc_data.ag_ratio;
            std::cout<<'\t'<<data.rc_data.ac_ratio;
            std::cout<<'\t'<<data.rc_data.gg_ratio;
            std::cout<<'\t'<<data.rc_data.ga_ratio;
            std::cout<<'\t'<<data.rc_data.gc_ratio;
            std::cout<<'\t'<<data.rc_data.ta_ratio;
            std::cout<<'\t'<<data.rc_data.tg_ratio;
            std::cout<<'\t'<<data.rc_data.cg_ratio;

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
        for(int i = 0 ; i< HAPLOTYPES; i++ ) tmp.hapCounts[i]=0;
        int prev_bin_index = 0 ;
        int bin_count[HAPLOTYPES] ;
        for( int j = 0 ; j< HAPLOTYPES; j++ ){
            bin_count[j] = 0 ;
        }
        int curr_bin_index = 0 ;
        MeanCalc mc[HAPLOTYPES];
        for(int i = 0 ; i <(int)seq.size()-g_K+1;i++){
            std::string kmer = seq.substr(i,g_K);
            curr_bin_index = i / g_binsize ;
            if( curr_bin_index != prev_bin_index ) {
                for( int j = 0 ; j< HAPLOTYPES; j++ ){
                    mc[j].Push(bin_count[j]);
                    bin_count[j] = 0 ;
                }
                prev_bin_index = curr_bin_index ;
            }
            for( int j = 0 ; j< HAPLOTYPES; j++ )
            {
                if( g_kmers[j].find(kmer) != g_kmers[j].end() ){
                    tmp.hapCounts[j] ++ ;
                    bin_count[j] ++ ;
                }
            }
        }
        for( int j = 0 ; j< HAPLOTYPES; j++ ){
            if( bin_count[j] > 0 || float ((seq.size()-g_K+1)%g_binsize)/float(g_binsize) > 0.5 ){
                mc[j].Push(bin_count[j]);
            }
            tmp.density[j]=mc[j].Mean(g_sd_fac);
        }
        tmp.rc_data = calc_read(seq);
        data.SaveOutput(tmp);
    }
    void submit(std::string & head ,std::string & seq, long long id){
        while( busy ) { std::this_thread::sleep_for(std::chrono::seconds(1));}
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
    std::cerr<<"Uasge :\n\tclassify_read --hap hap0.kmer --hap hap1.kmer --read read.fa [--read read_2.fa] [--thread t_num (8 default) ] [--format fasta/fastq (default fasta)] [--bin_size binsize (1000 default)] [--sd_fac sd-factor (default 2.0) ]"<<std::endl;
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
        {"hap",   required_argument, NULL, 'p'},
        {"read",  required_argument, NULL, 'r'},
        {"format",required_argument, NULL, 'f'},
        {"thread",required_argument, NULL, 't'},
        {"bin_size",required_argument,NULL, 'b'},
        {"sd_fac",required_argument,NULL, 's'},
        {"help",  no_argument,       NULL, 'h'},
        {0, 0, 0, 0}
    };
    static char optstring[] = "p:r:t:f:s:b:h";
    std::vector<std::string> haps;
    std::vector<std::string> read;
    std::string format="fasta";
    int t_num=8;
    while(1){
        int c = getopt_long(argc, argv, optstring, long_options, NULL);
        if (c<0) break;
        switch (c){
            case 's':
                g_sd_fac = atof(optarg);
                break;
            case 'b':
                g_binsize = atoi(optarg);
                break;
            case 'p':
                haps.push_back(std::string(optarg));
                break;
            case 'r':
                read.push_back(std::string(optarg));
                break;
            case 'f':
                format = std::string(optarg);
                break;
            case 't':
                t_num = atoi(optarg);
                break;
            case 'h':
            default :
                printUsage();
                return -1;
        }
    }
    if( haps.size() != HAPLOTYPES || read.empty()|| t_num< 1) {
        printUsage();
        return -1;
    }
    if( format != "fasta" && format != "fastq" ){
        std::cerr<<" ERROR : invalid format : ["<<format<<"] . exit ...\n";
        return -1;
    }
    logtime();
    std::cerr<<"__START__"<<std::endl;
    for( int i = 0 ; i<HAPLOTYPES ; i ++ ) {
        std::cerr<<"__load hap"<<i<<" kmers__"<<std::endl;
        logtime();
        load_kmers(haps[i],i);
    }
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
