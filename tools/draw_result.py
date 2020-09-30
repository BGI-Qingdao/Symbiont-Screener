import numpy as np
import matplotlib as mpl
#mpl.use('agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.decomposition import PCA
#np.seterr(invalid='ignore')
#plt.switch_backend("agg")
##############################
data=np.loadtxt("5r.matrix",dtype=float)
W=np.loadtxt("2mer.matrix",dtype=float)
Z=np.loadtxt("result.txt",dtype=float)

Y=data[:,0]
predict_Y=data[:,1]
X=data[:,2:]
X=np.hstack((X,W))
pca=PCA(n_components=20,whiten=True)
X2=pca.fit_transform(X)
print(pca.components_)
print(pca.explained_variance_ratio_)
np.random.seed(0)
idx=np.random.randint(Y.shape[0]-1,size=5000)
draw_point_Y=Z[idx]
draw_point_X=X2[idx,:]
predict_Y=predict_Y[idx]
p_class=np.where(predict_Y==1)
one_class=np.where(draw_point_Y==1)
class2=np.where(draw_point_Y==0)

fig = plt.figure()
plt.scatter(draw_point_X[class2,0],draw_point_X[class2,1],edgecolor='blue',c='',label='others',marker='.')
plt.scatter(draw_point_X[one_class,0],draw_point_X[one_class,1],edgecolor='red',c='',label='host',marker='.')
plt.scatter(draw_point_X[p_class,0],draw_point_X[p_class,1],c='k',marker='.' ,label='priori',edgecolor='')
#plt.legend()
handles, labels = plt.gca().get_legend_handles_labels()
order = [2,1,0]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
plt.xlabel("feature 1")
plt.ylabel("feature 2")
plt.savefig('./r_01_01.png')

fig = plt.figure()
plt.scatter(draw_point_X[one_class,0],draw_point_X[one_class,1],edgecolor='red',c='',label='host',marker='.',alpha=0.5)
plt.scatter(draw_point_X[p_class,0],draw_point_X[p_class,1],c='k',marker='.', label='priori',edgecolor='',alpha=0.5)
plt.scatter(draw_point_X[class2,0],draw_point_X[class2,1],edgecolor='blue',c='',label='others',marker='.',alpha=0.5)
#plt.legend()
handles, labels = plt.gca().get_legend_handles_labels()
order = [1,0,2]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
plt.xlabel("feature 1")
plt.ylabel("feature 2")
plt.savefig('./r_01_02.png')


fig = plt.figure()
plt.scatter(draw_point_X[class2,0],draw_point_X[class2,2],edgecolor='blue',c='',label='others',marker='.')
plt.scatter(draw_point_X[one_class,0],draw_point_X[one_class,2],edgecolor='red',c='',label='host',marker='.')
plt.scatter(draw_point_X[p_class,0],draw_point_X[p_class,2],c='k',marker='.', label='priori',edgecolor='')
plt.xlabel("feature 1")
plt.ylabel("feature 3")
plt.savefig('./r_02_01.png')

fig = plt.figure()
plt.scatter(draw_point_X[one_class,0],draw_point_X[one_class,2],edgecolor='red',c='',label='host',marker='.')
plt.scatter(draw_point_X[p_class,0],draw_point_X[p_class,2],c='k',marker='.' ,label='priori',edgecolor='')
plt.scatter(draw_point_X[class2,0],draw_point_X[class2,2],edgecolor='blue',c='',label='others',marker='.')
plt.xlabel("feature 1")
plt.ylabel("feature 3")
plt.savefig('./r_02_02.png')


fig = plt.figure()
plt.scatter(draw_point_X[class2,1],draw_point_X[class2,2],edgecolor='blue',c='',label='others',marker='.')
plt.scatter(draw_point_X[one_class,1],draw_point_X[one_class,2],edgecolor='red',c='',label='host',marker='.')
plt.scatter(draw_point_X[p_class,1],draw_point_X[p_class,2],c='k',marker='.', label='priori',edgecolor='')
plt.xlabel("feature 2")
plt.ylabel("feature 3")
plt.savefig('./r_12_01.png')

fig = plt.figure()
plt.scatter(draw_point_X[one_class,1],draw_point_X[one_class,2],edgecolor='red',c='',label='host',marker='.')
plt.scatter(draw_point_X[p_class,1],draw_point_X[p_class,2],c='k',marker='.', label='priori',edgecolor='')
plt.scatter(draw_point_X[class2,1],draw_point_X[class2,2],edgecolor='blue',c='',label='others',marker='.')
plt.xlabel("feature 2")
plt.ylabel("feature 3")
plt.savefig('./r_12_02.png')

#fig = plt.figure()
#ax = Axes3D(fig)
#plt.scatter(draw_point_X[class9,0],draw_point_X[class9,1],draw_point_X[class9,2],c='purple',marker='.')
#plt.scatter(draw_point_X[class8,0],draw_point_X[class8,1],draw_point_X[class8,2],c='purple',label='Y983&Y319',marker='.')
#plt.scatter(draw_point_X[class7,0],draw_point_X[class7,1],draw_point_X[class7,2],c='gold',marker='.')
#plt.scatter(draw_point_X[class6,0],draw_point_X[class6,1],draw_point_X[class6,2],c='gold',label='LB8&M71',marker='.')
#plt.scatter(draw_point_X[class5,0],draw_point_X[class5,1],draw_point_X[class5,2],c='green',marker='.')
#plt.scatter(draw_point_X[class4,0],draw_point_X[class4,1],draw_point_X[class4,2],c='green',label='43&44',marker='.')
#plt.scatter(draw_point_X[class3,0],draw_point_X[class3,1],draw_point_X[class3,2],c='blue',marker='.')
#plt.scatter(draw_point_X[class2,0],draw_point_X[class2,1],draw_point_X[class2,2],c='blue',label='Ecoli&R12',marker='.')
#plt.scatter(draw_point_X[one_class,0],draw_point_X[one_class,1],draw_point_X[one_class,2],c='red',label='chr19',marker='.')
#plt.show()
