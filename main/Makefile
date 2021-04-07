.PHONY: all clean

all: density_3lib density_3lib_stlfr gc_nmer gc_nmer_stlfr  kmer2strobemer

gc_nmer: gc_nmer.cpp gzstream/gzstream.C gzstream/gzstream.h
	g++ -c -g  gzstream/gzstream.C -I./gzstream -lz -o gzstream.o
	g++ -g -std=c++11 gc_nmer.cpp gzstream.o -lz -lpthread -o gc_nmer 

gc_nmer_stlfr: gc_nmer_stlfr.cpp gzstream/gzstream.C gzstream/gzstream.h
	g++ -c -g  gzstream/gzstream.C -I./gzstream -lz -o gzstream.o
	g++ -g -std=c++11 gc_nmer_stlfr.cpp gzstream.o -lz -lpthread -o gc_nmer_stlfr

density_3lib : density_3lib.cpp gzstream/gzstream.C gzstream/gzstream.h
	g++ -c -g  gzstream/gzstream.C -I./gzstream -lz -o gzstream.o
	g++ -g -std=c++11 density_3lib.cpp gzstream.o -lz -lpthread -o density_3lib


density_3lib_stlfr : density_3lib_stlfr.cpp gzstream/gzstream.C gzstream/gzstream.h kmer/kmer.h
	g++ -c -g  gzstream/gzstream.C -I./gzstream -lz -o gzstream.o
	g++ -g -std=c++11  density_3lib_stlfr.cpp gzstream.o -lz -lpthread -o density_3lib_stlfr

kmer2strobemer:
	g++ -c -g strobemer/strobemer.cpp -I./strobemer -o strobemer.o
	g++ -g kmer2strobemer.cpp strobemer.o -o kmer2strobemer

clean :
	rm -f density_3lib density_3lib_stlfr gc_nmer gc_nmer_stlfr
	rm -f *.o gzstream/*.o
