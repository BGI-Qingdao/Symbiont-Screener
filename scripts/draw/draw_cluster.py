#!/usr/bin/env python3

import numpy as np
import pandas as pd
import plotly.express as px

pdata=pd.read_csv('cluster_result.txt',sep='\t')
#formula_predict consensus_result        best_count      second_best_count
pdata['hit'] = pdata['best_count'] + pdata['second_best_count']

cache=np.zeros((2,11))
for _ , row in pdata.iterrows():
    if row['formula_predict']==1:
        cache[0,row['hit']]+=1
    else:
        cache[1,row['hit']]+=1

fdata=pd.DataFrame({'hit':
        [ 0,1,2,3,4,5,6,7,8,9,10,0,1,2,3,4,5,6,7,8,9,10],
        'type':['trio_only','trio_only','trio_only',
   'trio_only','trio_only','trio_only',
   'trio_only','trio_only','trio_only',
   'trio_only','trio_only',
   'cluster','cluster','cluster','cluster',
   'cluster','cluster','cluster','cluster',
   'cluster','cluster','cluster']})
fdata['count']=cache.reshape(-1)

fig = px.bar(fdata, x="hit", y="count",color='type')
fig.write_html('cluster.html')


