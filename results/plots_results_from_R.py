# -*- coding: utf-8 -*-
"""
Created on Wed Aug 23 13:09:37 2017

@author: cdonnat
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
name='change_pnt_SBM_85'
info='SBM graph, change point at T=6 and T=13'
dist=pd.DataFrame.from_csv('/Users/cdonnat/Dropbox/Distances/results/'+name+'.csv')
dist=dist.T
dist.index=range(dist.shape[0])
C=dist.corr('kendall')
sns.heatmap(C)
plt.savefig('/Users/cdonnat/Dropbox/Distances/write_up/plot/eigenT40_N'+name+'_correlation.pdf', bbox_inches='tight')

plt.figure()
cmap=plt.get_cmap('gnuplot')
colors=[cmap(xx) for xx in np.linspace(0,1,dist.shape[1])]

colors2=['red','blue','yellow','black','orange','purple','violet','green','lightblue','forrestgreen','cyan','gold','chocolate','indianred','mediumseagreen','magenta','grey','royalblue']

sns.set_context('paper',font_scale=1.3)
sns.set_style('white')
for i in range(6):
    if dist.columns.values[i]=='Spanning Trees':
        plt.plot(dist.iloc[:,i]/np.max(dist.iloc[:,i]),label=str(dist.columns.values[i])+ ' /' +'%.3f'%np.max(dist.iloc[:,i]),c=colors2[i])
    else:
        plt.plot(dist.iloc[:,i],label=dist.columns.values[i],c=colors2[i])
plt.plot(dist['ST norm'],label='ST norm',c=colors[9])
plt.xlabel('Time')
plt.ylabel('Distance')
plt.title('Distances for  small consecutive random changes,\n'+info)
plt.legend(loc='center left',bbox_to_anchor=(1,0.5))
plt.savefig('/Users/cdonnat/Dropbox/Distances/write_up/plot/eigenT40_N'+name+'_first_plots.pdf', bbox_inches='tight')

colors2=['red','blue','yellow','black','orange','purple','violet','green','lightblue','forrestgreen','cyan','gold','chocolate','indianred','mediumseagreen','magenta','grey','royalblue']
plt.figure()
sns.set_context('paper',font_scale=1.3)
sns.set_style('white')
for i in [6,7,8,10,11]:
        plt.plot(dist.iloc[:,i]/np.max(dist.iloc[:,i]),label=str(dist.columns.values[i])+ ' /' +'%.3f'%np.max(dist.iloc[:,i]),c=colors2[i])
plt.xlabel('Time')
plt.ylabel('Distance')
plt.title('Distances for  small consecutive random changes,\n'+info)
plt.legend(loc='center left',bbox_to_anchor=(1,0.5))
plt.savefig('/Users/cdonnat/Dropbox/Distances/write_up/plot/eigenT40_N'+name+'_plots_poly.pdf', bbox_inches='tight')


plt.figure()
sns.set_context('paper',font_scale=1.3)
sns.set_style('white')
for i in range(12,18):
        plt.plot(dist.iloc[:,i]/np.max(dist.iloc[:,i]),label=str(dist.columns.values[i])+ ' /' +'%.3f'%np.max(dist.iloc[:,i]),c=colors2[i])
plt.plot(dist['ST norm'],label='ST norm',c=colors[9])
plt.xlabel('Time')
plt.ylabel('Distance')
plt.title('Distances for  small consecutive random changes,\n'+info)
plt.legend(loc='center left',bbox_to_anchor=(1,0.5))
plt.savefig('/Users/cdonnat/Dropbox/Distances/write_up/plot/eigenT40_N'+name+'_plots_eigen_func.pdf', bbox_inches='tight')



