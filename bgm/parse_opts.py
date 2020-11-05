#!/usr/bin/python3
# -*- coding: UTF-8 -*-
import sys
import getopt
class OPTs:
    def __init__(self):
        self.trio_file=""
        self.mer2_file=""
        self.loop_num =30
        self.clusters_number=5
        self.debug = False

    def Usage(self):
        print("usage: main_logic.py -t <trio_file> -m <2mer2_file> [-l loop number (default 30)] [-c cluters number (default 5)]");

    def Parse(self,argv):
        print("Info : argv = %s " % str(argv), file=sys.stderr)
        try:
            opts, args = getopt.getopt(argv,"hdt:m:l:c:",["help","debug","trio=","2mer=","loop=","cluters="])
        except getopt.GetoptError:
            self.Usage()
            sys.exit(2)
        #print(opts)
        for opt, arg in opts:
            if opt in ('-d' , "--debug"):
                self.debug=True
            elif opt in ('-h' , "--help"):
                self.Usage()
                sys.exit()
            elif opt in ("-t", "--trio"):
                self.trio_file=arg
            elif opt in ("-c", "--clusters"):
                self.clusters_number=int(arg)
            elif opt in ("-m", "--2mer"):
                self.mer2_file=arg
            elif opt in ("-l", "--loop"):
                self.loop_num=int(arg)
        if( self.clusters_number < 2 ) :
            print("Error: -c must greater than 1")
            sys.exit(3)
        if( self.trio_file == '' or self.mer2_file == '' ):
            print("Error: -t and -m are necessary !!! exit ..." ,file=sys.stderr)
            sys.exit(2)

        if ( self.loop_num>=5 and self.loop_num < 20 ):
            print("WARN : loop_num is %d but >=20 is recommand."%self.loop_num ,file=sys.stderr)

        if (  self.loop_num<5 ):
            print("Error: loop_num is %d less than 5 !!! exit ..."%self.loop_num ,file=sys.stderr)
            sys.exit(4)
