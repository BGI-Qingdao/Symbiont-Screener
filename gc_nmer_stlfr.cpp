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

int g_nmer=2;

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

    void Init(){

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
    void Add(const NmerFreq & other){
        for( const auto & p : other.cache){
            if( cache.find(p.first ) == cache.end() )
                cache[p.first] = p.second ;
            else 
                cache[p.first] += p.second ;
        }
    }
};

int NmerFreq::n ;
std::set<std::string> NmerFreq::allmer ;

struct BarcodeGC_content {
    int total ;
    int gc ;
    void Init() { total = 0 ; gc = 0 ; }
    void GC_count(const std::string & read ) {
        for( char c : read ){
            total ++;
            if( isG(c) || isC(c) )
                gc ++ ;
        }
    }

    void Add(const BarcodeGC_content & other ){
        total += other.total ;
        gc += other.gc ;
    }

    float GC_content() const { 
        if (total == 0 ) return 0 ;
        return float(gc)/float(total);
    }
};

//
// Output cache
//
struct BarcodeInfo {
    std::string barcode_name;
    BarcodeGC_content gc;
    NmerFreq nmer_freq;
    void Add(const BarcodeInfo& other){
        assert(barcode_name==other.barcode_name);
        gc.Add(other.gc);
        nmer_freq.Add(other.nmer_freq);
    }
};
//
// barcode haplotype relate functions
//
struct BarcodeCache {
    std::map<std::string, BarcodeInfo> barcode_cache;
    void InitBarcode(const std::string & barcode){
        if( barcode_cache.find(barcode) == barcode_cache.end() ){
            barcode_cache[barcode].barcode_name = barcode ;
            barcode_cache[barcode].gc.Init();
            barcode_cache[barcode].nmer_freq.Init();
        }
    }
    void UpdateSeq(const std::string & barcode , const std::string & seq){
        barcode_cache[barcode].gc.GC_count(seq);
        barcode_cache[barcode].nmer_freq.DealSeq(seq);
    }
    void Add(const BarcodeCache & other){
        for(const auto & pair : other.barcode_cache) {
            InitBarcode(pair.first);
            barcode_cache[pair.first].Add(pair.second);
        }
    }
};

void printBarcodeInfos(const BarcodeCache& fdata){
    for( const auto & pair : fdata.barcode_cache){
        std::cout<<pair.first;
        std::cout<<'\t'<<pair.second.gc.GC_content();
        pair.second.nmer_freq.Print();
        std::cout<<'\n';
    }
}

//
//reads relate functions
//

// @return : barcode
// A stLFR read's head looks like :
//   @V300017823L1C001R051096800#203_1533_1069/1
//                barcode str :  203_1533_1069
std::string parseName(const std::string & head){
    int s=-1, e=-1;
    for( int i = 0 ; i< (int)head.size(); i++ ){
        if( head[i] == '#' ) s=i;
        if( head[i] == '/' ) e=i;
    }
    return head.substr(s+1,e-s-1);
}

#define max_buffer 1024
struct Buffer{
    std::array<std::string,max_buffer> heads;
    std::array<std::string,max_buffer> seqs ;
    void Init() { size = 0 ; }
    int size ;
};

struct MultiThread {
    int t_nums ;
    bool end;
    bool busy;
    void Worker(int index){
        //int miss = 0;
        //int hit = 0 ;
        Buffer buffer;
        //std::pair<std::string,std::string> job;
        while(true){
            locks[index].lock();
            if( caches[index].empty() ){
                //miss ++ ;
                busy = false ;
                locks[index].unlock();
                if(end) { 
                    //std::cerr<<"thread="<<index<<" miss="<<miss<<" hit="<<hit<<std::endl;
                    return ;
                }
                std::this_thread::sleep_for(std::chrono::microseconds(10));
                continue;
            }
            if( ! caches[index].empty() ){
                //hit ++ ;
                if( caches[index].size() > 300 ) busy = true ;
                else if ( caches[index].size() < 50 ) busy = false ;
                std::swap(buffer,caches[index].top());
                caches[index].pop();
                locks[index].unlock();
                for( int i = 0 ; i < buffer.size ; i ++ )
                    process_reads(buffer.heads.at(i),buffer.seqs.at(i),index);
            } else 
                locks[index].unlock();
        }
    }
    MultiThread(int t_num){
        t_nums = t_num ;
        barcode_caches = new BarcodeCache[t_num];
        locks = new std::mutex[t_num];
        threads = new std::thread*[t_num];
        busy = false;
        end=false;
        for(int i = 0 ; i< t_num ; i++){
            caches.push_back(std::stack<Buffer>());
            //caches.push_back(std::stack< std::pair<std::string,std::string> >());
            threads[i] = new std::thread([this,i](){int index=i ;Worker(index); });
        }
    }
    ~MultiThread(){
        delete [] threads;
        delete [] locks;
        delete [] barcode_caches;
    }
    void process_reads(const std::string & head ,
                         const std::string & read , int index) {
        int vote[3] ; vote[0]=0;vote[1]=0;vote[2]=0;
        std::string barcode = parseName(head);
        barcode_caches[index].InitBarcode(barcode);
        barcode_caches[index].UpdateSeq(barcode,read);
    }
    //void submit(const std::string & head ,const std::string & seq){
    void submit(const Buffer & buffer ){
        static long index = 0;
        while( busy ) { std::this_thread::sleep_for(std::chrono::seconds(1));}
        index ++ ;
        int id = index % t_nums;
        locks[id].lock();
        caches[id].push(buffer);//std::make_pair(head,seq));
        locks[id].unlock();
    }
    void wait(){
        for(int i = 0 ; i <t_nums; i++){
            threads[i]->join();
            delete threads[i];
        }
    }
    void collectBarcodes(){
        for(int i = 0 ; i<t_nums ;i++)
            final_data.Add(barcode_caches[i]);
    }
    //std::vector<std::stack< std::pair<std::string,std::string> >>  caches;
    std::vector<std::stack<Buffer> > caches;
    std::mutex * locks;
    std::thread ** threads; 
    BarcodeCache * barcode_caches;
    BarcodeCache final_data;
};

void processFastq(const std::string & file,int t_num,BarcodeCache& data){
    std::string head;
    std::string seq;
    std::string tmp;
    MultiThread mt(t_num);
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
    Buffer buffer ;
    buffer.Init();
    while(!std::getline(*in,head).eof()){
        std::getline(*in,seq);
        //mt.submit(head,seq);
        buffer.heads.at(buffer.size) = head ;
        buffer.seqs.at(buffer.size) = seq;
        buffer.size ++ ;
        if ( buffer.size >= max_buffer ){
            mt.submit(buffer);
            buffer.Init();
        }
        std::getline(*in,tmp);
        std::getline(*in,tmp);
    }
    if( buffer.size > 0 ){
        mt.submit(buffer);
        buffer.Init();
    }
    mt.end=true;
    mt.wait();
    mt.collectBarcodes();
    data.Add(mt.final_data);
}

void printUsage(){
    std::cerr<<"\n\
Uasge :\n\
    classify --hap0 hap0 --hap1 hap1 --read read1.fq [options]\n\
\n\
Options:\n\
        -h/--help                       print this uasge and exit.\n\
        -p/--hap0                       unshared kmer set of hap0.\n\
        -m/--hap1                       unshared kmer set of hap1.\n\
        -c/--hap2                       unshared kmer set of hap2.\n\
        -r/--read                       filial reads in fastq format. gzip file must be ended by \".gz\".\n\
        -t/--thread   (8 default)       thread number to used.\n\
        -w/--weight0  (1.0 default)     weight of hap0.\n\
        -u/--weight1  (1.0 default)     weight of hap1.\n\
        -x/--weight2  (1.0 default)     weight of hap2.\n\
        -f/--adaptor_f                  forward adaptor sequence.\n\
                                        default \"CTGTCTCTTATACACATCTTAGGAAGACAAGCACTGACGACATGA\"\n\
        -q/--adaptor_r                  reverse adaptor sequence.\n\
                                        default \"TCTGCTGAGTCGAGAACGTCTCTGTGAGCCAAGGAGTTGCTCTGG\"\n\
\n\
Examples:\n\
    ./classify --hap0 p.kmers --hap1 m.kmers --read input.fastq.gz\n\
\n\
    ./classify --hap0 p.kmers --hap1 m.kmers --read input.L01.fastq.gz --read input.L02.fastq.gz\n\
\n\
    ./classify --hap0 p.kmers --hap1 m.kmers --read input.L01.fastq.gz --read input.L02.fastq.gz -t 24 --weight1 1.04 -f CTGTCTCTTATACACATCTTAGGAAGACAA -q TCTGCTGAGTCGAGAACGTCTCTG\n\
\n\
Output format:\n\
barcode\thaplotype(0/1/-1)\tkmer_count_hap0\tkmer_count_hap1\n\
\n\
Usage done.\n\
";
}


void TestAll(){
    assert(parseName("VSDSDS#XXX_xxx_s/1")=="XXX_xxx_s");
    assert(parseName("VSDSDS#Libxxx_XXX_xxx_s/1")=="Libxxx_XXX_xxx_s");
    Kmer::InitFilter(5);
    auto str1=BaseStr::str2BaseStr("AGCTC");
    int  t1[] = { '\000','\003','\001','\002','\001'};
    for( int i = 0 ; i <5 ; i++ ) assert(str1[i] == t1[i]);
    auto str2 = BaseStr::str2BaseStr("GAGCT");
    int  t2[] = {'\003','\000','\003','\001','\002'};
    for( int i = 0 ; i <5 ; i++ ) assert(str2[i] == t2[i]);
                         // 00 1101 1001
    assert(Kmer::str2Kmer(BaseStr::str2BaseStr("AGCTC")).low == 0xD9);
    assert(Kmer::str2Kmer(BaseStr::str2BaseStr("AGCTC")).high == 0);
    assert(Kmer::str2Kmer(BaseStr::str2BaseStr("GAGCT")).low == 0xD9);
    assert(Kmer::str2Kmer(BaseStr::str2BaseStr("GAGCT")).high == 0);
    auto kmers = Kmer::chopRead2Kmer(BaseStr::str2BaseStr("GAGCTA"));
    assert(kmers.size() == 2 );
    // GAGCT->AGCTC 00 1101 1001 -> 0xD9
    // AGCTA        00 1101 1000 -> 0xD8
    assert(kmers[0].high == 0 );
    assert(kmers[0].low == 0xD9 );
    assert(kmers[1].high == 0 );
    assert(kmers[1].low == 0xD8 );
    std::string k1 = BaseStr::BaseStr2Str(Kmer::ToBaseStr(kmers[0]));
    std::string k2 = BaseStr::BaseStr2Str(Kmer::ToBaseStr(kmers[1]));
    assert(k1 == "AGCTC");
    assert(k2 == "AGCTA");
}

//
// Main function
//

int main(int argc ,char ** argv ){
    TestAll();
    static struct option long_options[] = {
        {"read", required_argument,  NULL, 'r'},
        {"thread",required_argument, NULL, 't'},
        {"help",  no_argument,       NULL, 'h'},
        {0, 0, 0, 0}
    };
    static char optstring[] = "r:t:h";
    std::vector<std::string> read;
    int t_num=8;
    while(1){
        int c = getopt_long(argc, argv, optstring, long_options, NULL);
        if (c<0) break;
        switch (c){
            case 'r':
                read.push_back(std::string(optarg));
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
    if(  read.empty()|| t_num< 1) {
        printUsage();
        return -1;
    }
    std::cerr<<"__START__"<<std::endl;
    logtime();
    BarcodeCache data;
    for(const auto r : read ){
        std::cerr<<"__process read: "<<r<<std::endl;
        processFastq(r,t_num,data);
        logtime();
        std::cerr<<"__process read done__"<<std::endl;
    }
    std::cerr<<"__print result__"<<std::endl;
    printBarcodeInfos(data);
    logtime();
    std::cerr<<"__END__"<<std::endl;
}
