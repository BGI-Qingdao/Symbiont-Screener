.PHONY: all clean cpbin

CC  ?= gcc
CXX ?= g++

CXXFLAGS += -std=c++11 -I./gzstream -I./strobemer

LD_FLAGS += -lz -lpthread

all: cpbin density_3lib gc_nmer  kmer2strobemer density_3lib_strobemer

cpbin: kmer2strobemer density_3lib_strobemer density_3lib gc_nmer
	mkdir -p ../symbiontscreenerbin
	cp kmer2strobemer density_3lib_strobemer density_3lib gc_nmer ../symbiontscreenerbin/

gzstream.o: gzstream/gzstream.C gzstream/gzstream.h
	${CXX} -c gzstream/gzstream.C ${CXXFLAGS} ${LD_FLAGS} -o gzstream.o

strobemer.o: strobemer/strobemer.cpp strobemer/strobemer.h
	${CXX} ${CXXFLAGS} -c strobemer/strobemer.cpp ${LD_FLAGS} -o strobemer.o

gc_nmer: gc_nmer.cpp gzstream.o
	${CXX} ${CXXFLAGS}  gc_nmer.cpp gzstream.o ${LD_FLAGS} -o gc_nmer

density_3lib : density_3lib.cpp gzstream.o
	${CXX} ${CXXFLAGS} density_3lib.cpp gzstream.o ${LD_FLAGS} -o density_3lib

kmer2strobemer: kmer2strobemer.cpp strobemer.o
	${CXX} ${CXXFLAGS} kmer2strobemer.cpp strobemer.o ${LD_FLAGS} -o kmer2strobemer

density_3lib_strobemer: density_3lib_strobemer.cpp gzstream.o strobemer.o
	${CXX} ${CXXFLAGS}  density_3lib_strobemer.cpp gzstream.o strobemer.o ${LD_FLAGS} -o density_3lib_strobemer

clean :
	rm -f density_3lib density_3lib_stlfr gc_nmer gc_nmer_stlfr kmer2strobemer density_3lib_strobemer
	rm -f *.o gzstream/*.o strobemer/*.o
