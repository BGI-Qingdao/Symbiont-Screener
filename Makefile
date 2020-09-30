 
.PHONY: all clean

all:
	cd main && make && cd -

clean:
	cd main && make clean
