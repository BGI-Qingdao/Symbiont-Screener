#ifndef __STROBEMER_H_
#define __STROBEMER_H_

#include <stdint.h>
#include <string>

// type of strobemer
enum strobemer_type{
    unknow=0,
    minstrobe=1,
    randstrobe=2
};

// k1-w-w format strobemer
//
// example :
//
// c++ code in example.cpp : 
// ``` 
// #include "strobemer.h"
// 
//#include "strobemer.h"
//
//int main() {
//     char seq[101]="ATGGGCAGAGTTTGACGTAGTCAATGCTTATGAACGAACGCTCCAATATGAATCAGCTCGTGATTTTTGCTGTAAAAATCGTAGCATACTGTTTGATAAA";
//     //strobemer::init(3,13,20,strobemer_type::minstrobe); // n=3,k=13,w1=w2=20,type=minstrobe
//     strobemer::init(3,13,20,strobemer_type::randstrobe);  // n=3,k=13,w1=w2=20,type=randstrobe
//     int number = 100-strobemer::strobmer_span()+1;
//     strobemer * buff = new strobemer[number];
//     strobemer::chop_strobemer(seq,100,buff);
//     for(int i = 0 ; i< number ; i++ ){
//         if(buff[i].valid)
//             buff[i].print();
//     }
//}
// ```
//
// command :
//
// g++ -c strobemer.cpp -o strobemer.o
// g++ example.cpp strobemer.o -o example
// ./example
//
class strobemer {
    public:
        static void init(int n, int k, int w , strobemer_type t) ;
        static void chop_strobemer(const char * seq,int len, strobemer * buff);
        static int strobmer_span() { return span ; }
    private:
        static int nkmer;
        static int ksize;
        static int wsize;
        static int span ; //W=k+2w=7+2*12=31
        static int kspan; //kspan=w-k+1; define a minimizer of how many kmers.
        static strobemer_type type;
    public:
        bool valid ;
        strobemer(); 
        ~strobemer();
        strobemer(const strobemer & other);
        strobemer & operator =(const strobemer & other);
        void print() const;
        std::string to_string() const ;
    private:
        char * kmer_forward ;//[nkmer][ksize];
    private:
        static void chop_minstrobe(const char * seq,int len, strobemer * buff);
        static void chop_randstrobe(const char * seq,int len, strobemer * buff);
}; // class strobemer


#endif //__MINSTROBE_H_
