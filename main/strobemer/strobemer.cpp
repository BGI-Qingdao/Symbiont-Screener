#include "strobemer.h"
#include <cassert>
#include <cstring>
#include <iostream>

/**********************************************************
 *
 * hash kmer into uint64
 *
 * *******************************************************/
// copy from minimap2:sketch.c :
static inline uint64_t hash64(uint64_t key, uint64_t mask)
{
    key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
    key = key ^ key >> 28;
    key = (key + (key << 31)) & mask;
    return key;
}//hash64

/**********************************************************
 *
 * char tables
 *
 * *******************************************************/
char complementary_base[256] = {
    'T', 'G', 'C', 'A',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'T', 'N', 'G',  'N', 'N', 'N', 'C',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'A', 'A', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'T', 'N', 'G',  'N', 'N', 'N', 'C',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'A', 'A', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N'
}; //complementary_base

static unsigned char seq_nt4_table[256] = {
    0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
}; //seq_nt4_table

/**********************************************************
 *
 * to get reverse-complementary sequences
 *
 * *******************************************************/
int inline rc_index(int index, int len){
    return len - 1 - index ;
}

void reverse_complete( const char * seq , int len, char * buff ){
    for( int i = 0 ; i < len ; i++ ){
        buff[i] = complementary_base[seq[rc_index(i,len)]];
    }
}

/**********************************************************
 *
 * to get hashed-kmers
 *
 * *******************************************************/
class binary_kmer {
    public:
        static int      ksize;
        static uint64_t kmask;
    public:
        static void InitK(int k) {
            assert(k>2 && k<33);
            ksize = k;
            kmask=(1ULL<<2*ksize) - 1;
        } // InitK

    public:
        binary_kmer() {
            kvalue=0; // forward kmer
            kspan = 0;
        } // Kmer
        void inline AddChar(char nucleotide){
            int c = seq_nt4_table[(uint8_t)nucleotide];
            if( c>3 ) { // hit invalid nucleotide
                kspan = 0;
                kvalue=0;
                return ;
            }
            kvalue = (kvalue << 2 | c) & kmask;
            kspan ++ ;
            if( kspan > ksize)  kspan = ksize ;;
        } // AddChar
        inline bool valid() const {
            return kspan>=ksize ;      /*check kmer length*/
        } // valid
        inline uint64_t khash64() const {
            assert(valid());
            return hash64(kvalue, kmask);
        } // khash64
    private:
        uint64_t kvalue;
        int      kspan;
}; // class binary_kmer

int      binary_kmer::ksize;
uint64_t binary_kmer::kmask;

struct hash_kmer {
    uint64_t hash_forward;
    bool valid;
    void init(const binary_kmer & bkmer){
        valid = bkmer.valid();
        if( ! valid ) return ;
        hash_forward = bkmer.khash64() ;
    }
    //
    // @Brief : chop sequence if length l into l-k+1 kmers.
    // @In    : seq and len
    // @Out   : buff
    // @Warn  : user must ensure the buff is large enough for l-k+1 kmers.
    //
    static void chop_kmers(const char * seq , int len, hash_kmer * buffer){
        binary_kmer bkmer;
        for( int i = 0 ; i<len ; i++ ){
            bkmer.AddChar(seq[i]);
            int start_pos = i-bkmer.ksize+1 ;
            if(start_pos >= 0){
                buffer[start_pos].init(bkmer);
            }
        }
    } //chop_kmers
}; // class hash_kmer



int strobemer::nkmer;
int strobemer::ksize;
int strobemer::wsize;
int strobemer::span ;
int strobemer::kspan;
strobemer_type strobemer::type ;

void strobemer::init(int n, int k , int w , strobemer_type t){
    assert(n>1);
    assert(w>k);
    assert(t == strobemer_type::minstrobe || t == strobemer_type::randstrobe );
    nkmer = n;
    ksize = k;
    wsize = w;
    span=k+(n-1)*w;
    kspan=w-k+1;
    type = t ;
    binary_kmer::InitK(k);
    std::cerr<<"INFO: init strobemer with"
        <<" n="<<n
        <<" k="<<k
        <<" w="<<w
        <<" t="<<((t==strobemer_type::minstrobe) ? "minstrobe" : "randstrobe")
        <<" span="<<span
        <<" kspan="<<kspan
        <<std::endl;
}

void strobemer::chop_strobemer(const char * seq,int len, strobemer * buff){
    if( type == strobemer_type::minstrobe ) {
        chop_minstrobe(seq,len,buff);
    } else if ( type == strobemer_type::randstrobe) {
        chop_randstrobe(seq,len,buff);
    } else {
        assert(0);
    }
}

void strobemer::chop_minstrobe(const char * seq,int len, strobemer * buff){
    // sanity check ...
    assert(seq!=nullptr);
    assert(len>=strobemer::span);
    // alloc buffers
    char * rc_seq = new char[len];
    hash_kmer* hkmer_buffer = new hash_kmer[len-binary_kmer::ksize+1];
    // sanity check ...
    assert(rc_seq!=nullptr);
    assert(hkmer_buffer!=nullptr);
    // fill buffers
    reverse_complete(seq,len,rc_seq);
    hash_kmer::chop_kmers(seq,len,hkmer_buffer);
    // clean buffer
    for( int i=0 ; i<=len-strobemer::span; i++ ){
        buff[i].valid = false;
    }

    for( int i=0 ; i<=len-strobemer::span; i++ ){
        // [i,w1_s) [w1_s , w1_e) [w1_e,w2_e)
        bool valid = true;
        // check the validation of span area
        for( int j = i ; j < i+strobemer::kspan; j++ ){
            if( ! hkmer_buffer[j].valid ) {
                i = j ; // iterator jump all infected area.
                valid =false ;
                break;
            }
        }
        if( ! valid ) continue ;
        buff[i].valid = true;
        //
        // -----------------------------------------------------------------------
        //                                     k            w       ....
        // 1. construct forward strobemer in [i,w1_s) [w1_s , w1_e) ....
        int p0 = i ;
        strncpy(buff[i].kmer_forward, &seq[p0],strobemer::ksize);
        // 1.1 find the ki minimizer in [wi_s , wi_e)
        int wi_s = i+strobemer::ksize;
        int wi_e = wi_s+strobemer::kspan ;
        int k_shift = strobemer::ksize;
        for( int ki=1; ki<strobemer::nkmer ; ki++){
            int p_next = -1;
            uint64_t temp_min ;
            uint64_t h_now;
            for( int j = wi_s ; j < wi_e ; j++ ){
                h_now = hkmer_buffer[j].hash_forward ; // minstrobe calculate hash value independantly
                if ( p_next == -1 || h_now < temp_min ) {
                    temp_min = h_now;
                    p_next = j ;
                }
            }
            strncpy(&buff[i].kmer_forward[k_shift], &seq[p_next],strobemer::ksize);
            wi_s=wi_s+strobemer::wsize;
            wi_e=wi_s+strobemer::kspan;
            k_shift += strobemer::ksize;
        }
    }
    // free memory
    delete [] hkmer_buffer;
    delete [] rc_seq;
}

void strobemer::chop_randstrobe(const char * seq,int len, strobemer * buff){
    // sanity check ...
    assert(seq!=nullptr);
    assert(len>=strobemer::span);
    // alloc buffers
    char * rc_seq = new char[len];
    hash_kmer* hkmer_buffer = new hash_kmer[len-binary_kmer::ksize+1];
    // sanity check ...
    assert(rc_seq!=nullptr);
    assert(hkmer_buffer!=nullptr);
    // fill buffers
    reverse_complete(seq,len,rc_seq);
    hash_kmer::chop_kmers(seq,len,hkmer_buffer);
    // clean buffer
    for( int i=0 ; i<=len-strobemer::span; i++ ){
        buff[i].valid = false;
    }

    for( int i=0 ; i<=len-strobemer::span; i++ ){
        // [i,w1_s) [w1_s , w1_e) [w1_e,w2_e)
        bool valid = true;
        // check the validation of span area
        for( int j = i ; j < i+strobemer::kspan; j++ ){
            if( ! hkmer_buffer[j].valid ) {
                i = j ; // iterator jump all infected area.
                valid =false ;
                break;
            }
        }
        if( ! valid ) continue ;
        buff[i].valid = true;
        //
        // -----------------------------------------------------------------------
        //                                     k            w       ....
        // 1. construct forward strobemer in [i,w1_s) [w1_s , w1_e) ....
        int p0 = i ;
        strncpy(buff[i].kmer_forward, &seq[p0],strobemer::ksize);
        uint64_t h_prev = hkmer_buffer[i].hash_forward ;
        // 1.1 find the ki minimizer in [wi_s , wi_e)
        int wi_s = i+strobemer::ksize;
        int wi_e = wi_s+strobemer::kspan ;
        int k_shift = strobemer::ksize;
        for( int ki=1; ki<strobemer::nkmer ; ki++){
            int p_next = -1;
            uint64_t temp_min ;
            uint64_t h_now;
            for( int j = wi_s ; j < wi_e ; j++ ){
                h_now = h_prev ^hkmer_buffer[j].hash_forward ;
                if ( p_next == -1 || h_now < temp_min ) {
                    temp_min = h_now;
                    p_next = j ;
                }
            }
            strncpy(&buff[i].kmer_forward[k_shift], &seq[p_next],strobemer::ksize);
            wi_s=wi_s+strobemer::wsize;
            wi_e=wi_s+strobemer::kspan;
            k_shift += strobemer::ksize;
            h_prev = h_now ;
        }
    }
    // free memory
    delete [] hkmer_buffer;
    delete [] rc_seq;
}

strobemer::strobemer() {
    kmer_forward = new char[nkmer*ksize];
    valid = false ;
}
strobemer::~strobemer(){
    delete [] kmer_forward;
}
strobemer::strobemer(const strobemer & other){
    kmer_forward = new char[nkmer*ksize];
    strncpy(kmer_forward,other.kmer_forward,nkmer*ksize);
}
strobemer & strobemer::operator =(const strobemer & other){
    if( &other == this) return *this;
    assert(kmer_forward!=nullptr);
    strncpy(kmer_forward,other.kmer_forward,nkmer*ksize);
    return *this;
}
void strobemer::print() const{
    for(int i = 0 ; i <nkmer*ksize; i++)
        std::cout<<kmer_forward[i];
    std::cout<<'\n';
}

std::string strobemer::to_string() const {
    return std::string(kmer_forward,0,nkmer*ksize);
}
