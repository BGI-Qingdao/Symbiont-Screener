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
Y=data[:,0]
X=data[:,2:]
X=np.hstack((X,W))
pca=PCA(n_components=20,whiten=True)
X2=pca.fit_transform(X)
print(pca.components_)
print(pca.explained_variance_ratio_)

np.random.seed(0)
idx=np.random.randint(Y.shape[0]-1,size=5000)
draw_point_Y=Y[idx]
draw_point_X=X2[idx,:]
classes=[];
color=['black','red','darkorange','gold','green','lime','turquoise','darkcyan','royalblue','navy']
labels= ['none','chr19','E.coli'   ,'R12','43-1A','44A' ,'LB8'      ,'M714'    ,'Y319'     ,'Y983']
for i in range(0,10):
    classes.append(np.where(draw_point_Y==i))

fig = plt.figure()
for i in range(1,10):
    indexs=classes[i]
    plt.scatter(draw_point_X[indexs,0],draw_point_X[indexs,1],c=color[i],label=labels[i],marker='.')
plt.legend()
plt.xlabel("feature 1")
plt.ylabel("feature 2")
plt.savefig('./01_01.png')

fig = plt.figure()
for i in range(1,10):
    indexs=classes[10-i]
    plt.scatter(draw_point_X[indexs,0],draw_point_X[indexs,1],c=color[10-i],label=labels[10-i], marker='.')
plt.legend()
plt.xlabel("feature 1")
plt.ylabel("feature 2")
handles, labels = plt.gca().get_legend_handles_labels()
order = [8,7,6,5,4,3,2,1,0]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
plt.savefig('./01_02.png')

fig = plt.figure()
for i in range(1,10):
    indexs=classes[i]
    plt.scatter(draw_point_X[indexs,0],draw_point_X[indexs,2],c=color[i],marker='.')
#plt.legend()
plt.xlabel("feature 1")
plt.ylabel("feature 3")
plt.savefig('./02_01.png')

fig = plt.figure()
for i in range(1,10):
    indexs=classes[10-i]
    plt.scatter(draw_point_X[indexs,0],draw_point_X[indexs,2],c=color[10-i],marker='.')
#plt.legend()
plt.xlabel("feature 1")
plt.ylabel("feature 3")
plt.savefig('./02_02.png')

fig = plt.figure()
for i in range(1,10):
    indexs=classes[i]
    plt.scatter(draw_point_X[indexs,1],draw_point_X[indexs,2],c=color[i],marker='.')
#plt.legend()
plt.xlabel("feature 2")
plt.ylabel("feature 3")
plt.savefig('./12_01.png')

fig = plt.figure()
for i in range(1,10):
    indexs=classes[10-i]
    plt.scatter(draw_point_X[indexs,1],draw_point_X[indexs,2],c=color[10-i],marker='.')
#plt.legend()
plt.xlabel("feature 2")
plt.ylabel("feature 3")
plt.savefig('./12_02.png')

##fig = plt.figure()
##ax = Axes3D(fig)
##plt.scatter(draw_point_X[class9,0],draw_point_X[class9,1],draw_point_X[class9,2],c='purple',marker='.')
##plt.scatter(draw_point_X[class8,0],draw_point_X[class8,1],draw_point_X[class8,2],c='purple',label='Y983&Y319',marker='.')
##plt.scatter(draw_point_X[class7,0],draw_point_X[class7,1],draw_point_X[class7,2],c='gold',marker='.')
##plt.scatter(draw_point_X[class6,0],draw_point_X[class6,1],draw_point_X[class6,2],c='gold',label='LB8&M71',marker='.')
##plt.scatter(draw_point_X[class5,0],draw_point_X[class5,1],draw_point_X[class5,2],c='green',marker='.')
##plt.scatter(draw_point_X[class4,0],draw_point_X[class4,1],draw_point_X[class4,2],c='green',label='43&44',marker='.')
##plt.scatter(draw_point_X[class3,0],draw_point_X[class3,1],draw_point_X[class3,2],c='blue',marker='.')
##plt.scatter(draw_point_X[class2,0],draw_point_X[class2,1],draw_point_X[class2,2],c='blue',label='Ecoli&R12',marker='.')
##plt.scatter(draw_point_X[one_class,0],draw_point_X[one_class,1],draw_point_X[one_class,2],c='red',label='chr19',marker='.')
##plt.show()
