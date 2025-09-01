import numpy as np

class MatrixLoader:

    def __init__( self, trio_file , mer2_file ,debug = False ):
        self.trio = trio_file
        self.mer2 = mer2_file
        self.debug = debug

    def Load(self):
        self.T = np.loadtxt(self.trio,dtype=float)
        self.M = np.loadtxt(self.mer2,dtype=float)
        if not self.debug :
            self.formula_predict = self.T[:,0]
            self.X = np.hstack( (self.T[:,1:] ,self.M ) )
        else:
            self.formula_predict = self.T[:,1]
            self.X = np.hstack( (self.T[:,2:] ,self.M ) )
            self.Y =  self.T[:,0]
