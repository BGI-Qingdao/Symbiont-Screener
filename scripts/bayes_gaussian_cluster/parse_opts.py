#!/usr/bin/python3
# -*- coding: UTF-8 -*-
import sys
import getopt
class OPTs:
    def __init__(self):
        self.trio_file=""
        self.mer2_file=""
        self.loop_num =10
        self.threads=30
        self.rseed=42
        self.clusters_number=5
        self.debug = False

    def Usage(self):
        print("usage: main_logic.py -t <trio_matrix> -m <3mer_matrix> [-l loop number (default 10)] [-c cluters number (default 5)] [-r random seed(default(42)] [-s threads(default 30)]");

    def Parse(self,argv):
        print("Info : argv = %s " % str(argv), file=sys.stderr)
        try:
            opts, args = getopt.getopt(argv,"hdt:m:l:c:r:s:",["help","debug","trio=","2mer=","loop=","cluters=","rseed=","threads="])
        except getopt.GetoptError:
            self.Usage()
            sys.exit(2)
        for opt, arg in opts:
            if opt in ('-r' , "--rseed"):
                self.rseed=int(arg)
            if opt in ('-d' , "--debug"):
                self.debug=True
            elif opt in ('-h' , "--help"):
                self.Usage()
                sys.exit()
            elif opt in ("-s", "--threads"):
                self.threads=arg
            elif opt in ("-t", "--trio"):
                self.trio_file=arg
            elif opt in ("-c", "--clusters"):
                self.clusters_number=int(arg)
            elif opt in ("-m", "--2mer"):
                self.mer2_file=arg
            elif opt in ("-l", "--loop"):
                self.loop_num=int(arg)
        if( self.clusters_number < 1 ) :
            print("Error: -c must greater than 0")
            sys.exit(3)
        if( self.trio_file == '' or self.mer2_file == '' ):
            print("Error: -t and -m are necessary !!! exit ..." ,file=sys.stderr)
            sys.exit(2)
