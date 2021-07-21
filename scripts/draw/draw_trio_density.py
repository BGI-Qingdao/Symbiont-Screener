#!/usr/bin/env python3
import numpy as np
import pandas as pd
import plotly.express as px

df=pd.read_csv('trio_density.data.txt',header=None,sep="\t")
df.columns=['read','readname','t1','length','p1','p2','p3','d1','d2','d3','n1','n2','n3']

df[['d1','d2','d3']]=  df[['d1','d2','d3']].applymap(np.int64)

df['t1'] = df[["d1", "d2"]].max(axis=1)
df['t2'] = df['d3']

t1_max = int(np.max(df['t1']))
t2_max = int(np.max(df['t2']))

t1_range = t1_max if t1_max<30  else 30
t2_range = t2_max if t2_max<300 else 300

fig = px.density_heatmap(df, x="t2",
                y="t1",
                marginal_x="histogram",
                marginal_y="histogram",
                nbinsy=t1_max//1+1,
                nbinsx=t2_max//10+1,
                range_x=(0,t2_range),
                range_y=(0,t1_range),
                )
fig.write_html('heatmap_t1_t2.html')
