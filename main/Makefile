.PHONY: all clean

all:   gc_nmer   kmer2strobemer density_3lib_strobemer

gc_nmer: gc_nmer.cpp gzstream/gzstream.C gzstream/gzstream.h
	g++ -c -g  gzstream/gzstream.C -I./gzstream -lz -o gzstream.o
	g++ -g -std=c++11 gc_nmer.cpp gzstream.o -lz -lpthread -o gc_nmer 


kmer2strobemer: kmer2strobemer.cpp strobemer/strobemer.cpp strobemer/strobemer.h
	g++ -c -std=c++11 -g strobemer/strobemer.cpp -I./strobemer -o strobemer.o
	g++ -g -std=c++11 kmer2strobemer.cpp strobemer.o -o kmer2strobemer

density_3lib_strobemer: density_3lib_strobemer.cpp gzstream/gzstream.C gzstream/gzstream.h strobemer/strobemer.cpp strobemer/strobemer.h
	g++ -c -g  gzstream/gzstream.C -I./gzstream -lz -o gzstream.o
	g++ -c -std=c++11 -g strobemer/strobemer.cpp -I./strobemer -o strobemer.o
	g++ -g -std=c++11 density_3lib_strobemer.cpp gzstream.o strobemer.o -lz -lpthread -o density_3lib_strobemer
clean :
	rm -f   gc_nmer  kmer2strobemer density_3lib_strobemer
	rm -f *.o gzstream/*.o strobemer/*.o
