#!/usr/bin/env python3

import pandas as pd
import plotly.express as px

pdata=pd.read_csv('paternal.histo',header=None,sep=' ')
pdata.columns=['frequency_of_kmers', 'number_of_distinct_kmers']
pdata['source']='paternal'

mdata=pd.read_csv('maternal.histo',header=None,sep=' ')
mdata.columns=['frequency_of_kmers', 'number_of_distinct_kmers']
mdata['source']='maternal'

pmdata=pd.concat([pdata,mdata],ignore_index=True)
fig = px.bar(pmdata, x="frequency_of_kmers", y=["number_of_distinct_kmers"],color="source", barmode="group" ,title="kmer frequency histogram",range_x=(0,200))

fig.write_html('kmer_frequency.html')
