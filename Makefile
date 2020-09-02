.PHONY: all clean

all: classify_3lib classify_stlfr_3lib gc_nmer

gc_nmer: gc_nmer.cpp gzstream/gzstream.C gzstream/gzstream.h
	g++ -c -g  gzstream/gzstream.C -I./gzstream -lz -o gzstream.o
	g++ -g -std=c++11 gc_nmer.cpp gzstream.o -lz -lpthread -o gc_nmer 


classify_3lib : classify_3lib.cpp gzstream/gzstream.C gzstream/gzstream.h
	g++ -c -g  gzstream/gzstream.C -I./gzstream -lz -o gzstream.o
	g++ -g -std=c++11 classify_3lib.cpp gzstream.o -lz -lpthread -o classify_3lib


classify_stlfr_3lib : classify_stlfr_3lib.cpp gzstream/gzstream.C gzstream/gzstream.h kmer/kmer.h
	g++ -c -g  gzstream/gzstream.C -I./gzstream -lz -o gzstream.o
	g++ -g -std=c++11  classify_stlfr_3lib.cpp gzstream.o -lz -lpthread -o classify_stlfr_3lib


clean :
	rm -f classify_3lib classify_stlfr_3lib
	rm -f *.o gzstream/*.o
