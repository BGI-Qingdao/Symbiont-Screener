import numpy as np
import matplotlib.pyplot as plt
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.datasets import load_iris
from sklearn.feature_selection import SelectFromModel

X=np.loadtxt("5r.matrix",dtype=float)
W=np.loadtxt("2mer.matrix",dtype=float)
Y=X[:,0]
Y=Y==1
X=X[:,2:]
X=np.hstack((X,W))
forest = ExtraTreesClassifier(n_estimators=50)
forest.fit(X, Y)
importances = forest.feature_importances_
std = np.std([tree.feature_importances_ for tree in forest.estimators_],
axis=0)
print(forest.feature_importances_)
print(std)
fig=plt.figure()
plt.title("Feature importances")
plt.bar(range(X.shape[1]), importances,color="r", yerr=std, align="center")
plt.xticks(range(0,20), ['F-only','M-only'
,'Shared'
,'GC-content'
,
'AA',
'AC',
'AG',
'AT',
'CA',
'CC',
'CG',
'CT',
'GA',
'GC',
'GG',
'GT',
'TA',
'TC',
'TG',
'TT'  ],rotation=66)
fig.subplots_adjust(left=0.1 , right=0.95, bottom=0.18, top=0.9)#,wspace=0.5,hspace=0.5)
plt.show()
