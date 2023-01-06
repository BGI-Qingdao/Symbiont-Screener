.PHONY:  all symbiontscreenerscripts symbiontscreenerscripts symbiontscreenerbinary

all:  symbiontscreenerbin symbiontscreenerscripts symbiontscreenerbinary 

symbiontscreenerbin:
	mkdir symbiontscreenerbin

symbiontscreenerbinary:
	cd sources && make CC=${CC} CXX=${CXX}

symbiontscreenerscripts:
	cp -r scripts/* symbiontscreenerbin/

condainstall: all
	echo  "Installed into ${PREFIX}"
	mkdir -p  ${PREFIX}/bin
	cp sysc ${PREFIX}/bin
	cp -r symbiontscreenerbin ${PREFIX}/bin/
	chmod a+x ${PREFIX}/bin/sysc
	chmod a+x ${PREFIX}/bin/symbiontscreenerbin/*
	chmod a+x ${PREFIX}/bin/symbiontscreenerbin/*/*.sh
	chmod a+x ${PREFIX}/bin/symbiontscreenerbin/*/*.py
	mkdir -p  ${PREFIX}/info
	cp LICENSE ${PREFIX}/info/

