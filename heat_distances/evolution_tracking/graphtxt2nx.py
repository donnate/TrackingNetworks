# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 18:31:43 2017

@author: cdonnat
This script takes a sequence of  graphs produced by one of the functions in main_tests.R and
outputs the distances as computed by the heat diffusion method.
"""

import numpy as np
import networkx as nx
import sys
import seaborn as sb
sys.path.append( '/Users/cdonnat/Dropbox/Distances/heat_distances/')
from shapes.shapes import *
from heat_diffusion import *
from distances_between_graphs import *
from graph_tools import *
from evolution_tracking.distances_evolution  import *

#name='smooth_rdm_30_0.4'
name='test_smooth_realistic_changes_60_2_Island'
type_g='Island'
b=0
#name='test_smooth_realistic_changes_60_2_Island_'+str(b+1)+'distances.txt'
path2graph='/Users/cdonnat/Dropbox/Distances/tests_synthetic_data/generated_graphs/'
graph_name='test_smooth_realistic_changes_60_2_Island_'+str(b+1)+'.txt'
path2figs='/Users/cdonnat/Dropbox/Distances/tests_synthetic_data/figs/'
figname=type_g+'AUCdistances.pdf'
D=pd.DataFrame.from_csv(path2graph+graph_name,header=None,index_col=None,sep=' ')
N=int(np.sqrt(D.shape[1]))
T=D.shape[0]
NetSeq={t:nx.from_numpy_matrix(D.iloc[t,:].values.reshape((N,N))) for t in range(T)}

taus=range(1,25,2)
distances,distances_chi,distances_chi_agg,coeff,chi_agg=track_evolution(NetSeq,diffusion_type="heat",tau=taus)
cmap=sb.plt.get_cmap('gnuplot')
x=np.linspace(0,1,len(taus))
colors={taus[i]: cmap(x[i]) for i in range(len(taus))}
fig, ax = plt.subplots()
sb.set(style='ticks')
sb.set_context("paper", font_scale=1.5)
sb.set_style('white')
for tau in taus:
    #distances[tau],distances_chi[tau],coeff[tau],chi_agg[tau]=track_evolution(NetSeq,diffusion_type="heat",tau=tau)
    plt.plot(distances[tau],c=colors[tau],label=tau)
plt.title('Evolution of the AUC Distances between consecutive \n graphs for different values of the scale')
plt.legend(loc='center left', bbox_to_anchor=(1.0,0.5))
plt.xlabel('Time t')
plt.ylabel(r'Distance between graphs $G_t$ and $G_{t+1}$')
plt.tight_layout()
ax.spines['left'].set_position('zero')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_position('zero')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.savefig(path2figs+figname, bbox_inches='tight')



figname=type_g+name+'Chidistances'+type_g+'.pdf'
fig, ax = plt.subplots()
cmap=sb.plt.get_cmap('gnuplot')
x=np.linspace(0,1,len(taus))
colors={taus[i]: cmap(x[i]) for i in range(len(taus))}
sb.set(style='ticks')
sb.set_context("paper", font_scale=1.5)
sb.set_style('white')
for tau in taus:
    #distances[tau],distances_chi[tau],coeff[tau],chi_agg[tau]=track_evolution(NetSeq,diffusion_type="heat",tau=tau)
    plt.plot(distances_chi[tau],c=colors[tau],label=tau)
plt.plot([1.0/len(taus)*d for d in distances_chi_agg],c='magenta',label="aggregated characteristic vector")
plt.title('Evolution of the Characteristic Distances between consecutive \n graphs for different values of the scale')
plt.legend(loc='center left', bbox_to_anchor=(1.0,0.5))
plt.xlabel('Time t')
plt.ylabel(r'Distance between graphs $G_t$ and $G_{t+1}$')
plt.tight_layout()
ax.spines['left'].set_position('zero')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_position('zero')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.savefig(path2figs+figname, bbox_inches='tight')



var={}

### Assess the variability of the changes
B=15

distances={tau:[] for tau in taus}
distances_chi={tau:[] for tau in taus}
distances_chi_agg=[]
for b in range(B):
    name='80_0.2_boot'+str(b+1)+'distances.txt'
    path2graph='/Users/cdonnat/Dropbox/Distances/tests_synthetic_data/generated_graphs/'
    graph_name='80_0.2_boot'+str(b+1)+'.txt'
    path2figs='/Users/cdonnat/Dropbox/Distances/tests_synthetic_data/figs/'
    D=pd.DataFrame.from_csv(path2graph+graph_name,header=None,index_col=None,sep=',')
    N=int(np.sqrt(D.shape[1]))
    T=D.shape[0]
    NetSeq={t:nx.from_numpy_matrix(D.iloc[t,:].values.reshape((N,N))) for t in range(T)}
    taus=range(1,25,2)
    dist,dist_chi,dist_chi_agg,coeff,chi_agg=track_evolution(NetSeq,diffusion_type="heat",tau=taus)
    for tau in taus:
        distances[tau]+=dist[tau]
        distances_chi[tau]+=dist_chi[tau]
    distances_chi_agg+=dist_chi_agg

fig, ax = plt.subplots()
sb.set_context("paper", font_scale=1.5)
sb.set_style('white')
figname='Hist_chi_agg_'+type_g+'.pdf'
plt.hist(distances_chi_agg,bins=10)   
plt.title('Variability of the characteristic distances \n between consecutive graphs  (60 samples) ')
plt.legend(loc='center left', bbox_to_anchor=(1.0,0.5))
plt.ylabel('Number of samples')
plt.xlabel(r'Distance between graphs $G_t$ and $G_{t+1}$')
plt.tight_layout()
ax.spines['left'].set_position('zero')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_position('zero')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.savefig(path2figs+figname, bbox_inches='tight')     
   



sharpe_agg=np.std(distances_chi_agg)/np.mean(distances_chi_agg)
sharpe_chi=[np.std(distances_chi[k])/np.mean(distances_chi[k]) for k in distances.keys()]
sharpe=[np.std(distances[k])/np.mean(distances[k]) for k in distances.keys()]
fig,ax=plt.subplots()
sb.set_style('white')
plt.plot(taus,sharpe,label='AUC',c='black')
plt.plot(taus,sharpe_chi,label='chi',color='red')
plt.scatter(taus[-1]+1,sharpe_agg,label='chi agg',color='blue',s=100)
plt.xlabel('Scale')
plt.ylabel('Std/Mean')
plt.tight_layout()
ax.spines['left'].set_position('zero')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_position('zero')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.legend(loc='center left', bbox_to_anchor=(1.0,0.5))
plt.savefig(path2figs+type_g+'sharpe.pdf', bbox_inches='tight')    



fig, ax = plt.subplots()
figname='Hist_all_chi_'+type_g+'.pdf'
for t in taus:
    plt.hist(distances_chi[t],bins=10,color=colors[t],label=str(t), alpha=0.5)
plt.title('Variability of the characteristic distances \n between consecutive graphs  (60 samples) ')
plt.legend(loc='center left', bbox_to_anchor=(1.0,0.5))
plt.ylabel('Number of samples')
plt.xlabel(r'Distance between graphs $G_t$ and $G_{t+1}$')
plt.tight_layout()
ax.spines['left'].set_position('zero')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_position('zero')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.savefig(path2figs+figname, bbox_inches='tight') 


taus=range(1,25,2)
fig, ax = plt.subplots()
plt.boxplot([distances_chi[t] for t in taus]+[[1.0/len(taus)*d for d in distances_chi_agg]])
figname='Boxplot_all_chi_'+type_g+'.pdf'
plt.xlabel('Scale')
plt.ylabel(r'Distance between graphs $G_t$ and $G_{t+1}$')
plt.tight_layout()
a=[t for t in taus]
a+=['agg. chi']
ax.set_xticklabels(a)
ax.spines['left'].set_position('zero')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_position('zero')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.savefig(path2figs+figname, bbox_inches='tight') 




fig, ax = plt.subplots()
plt.boxplot([distances[t] for t in taus])
figname='Boxplot_all_AUC_'+type_g+'.pdf'
plt.xlabel('Scale')
plt.ylabel(r'Distance between graphs $G_t$ and $G_{t+1}$')
plt.tight_layout()
ax.set_xticklabels(taus)
ax.spines['left'].set_position('zero')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_position('zero')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.savefig(path2figs+figname, bbox_inches='tight') 




name='test_smooth_realistic_changes_60_4_SBM_3'
dist=pd.DataFrame.from_csv(path2graph+name+'distances.txt',index_col=0,header=None,sep=' ')
dist=dist.iloc[1:,:]
dist.columns=range(dist.shape[1])
#info=type_g0+' graph, change point at T=9'
#info=type_g0+' graph, change point at T=6 and T=13'
info="SBM graph, smooth evolution"
for i in range(dist.shape[1]):
    dist[i]=pd.to_numeric(dist[i])
dist=dist.T
plt.figure()
figname=type_g+'other_distances.pdf'
colors2=['red','blue','yellow','black','orange','purple','violet','green','lightblue','forrestgreen','cyan','gold','chocolate','indianred','mediumseagreen','magenta','grey','royalblue']
sb.set_context('paper',font_scale=1.3)
sb.set_style('white')
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
plt.savefig(path2figs+figname, bbox_inches='tight')

colors2=['red','blue','yellow','black','orange','purple','violet','green','lightblue','green','cyan','gold','chocolate','indianred','mediumseagreen','magenta','grey','royalblue']
plt.figure()
figname=type_g+'other_distances2.pdf'
sb.set_context('paper',font_scale=1.3)
sb.set_style('white')
for i in [6,7,8,10,11]:
        plt.plot(dist.iloc[:,i]/np.max(dist.iloc[:,i]),label=str(dist.columns.values[i])+ ' /' +'%.3f'%np.max(dist.iloc[:,i]),c=colors2[i])
plt.xlabel('Time')
plt.ylabel('Distance')
plt.title('Distances for  small consecutive random changes,\n'+info)
plt.legend(loc='center left',bbox_to_anchor=(1,0.5))
plt.savefig(path2figs+figname, bbox_inches='tight')


colors2=['red','blue','yellow','black','orange','purple','violet','lightblue','green','cyan','gold','chocolate','indianred','mediumseagreen','magenta','grey','royalblue']

plt.figure()
figname=type_g+'other_distances3.pdf'
sb.set_context('paper',font_scale=1.3)
sb.set_style('white')
it_c=0
for column in ['Jaccard', 'Hamming', 'IM', 'HIM','alpha=0.5,order=5', 'alpha=0.9,order=3','f(l)=exp(-1.2l)']:
        plt.plot(dist[column]/np.max(dist[column]),label=str(column)+ ' /' +'%.3f'%np.max(dist[column]),color=colors2[it_c])
        it_c+=1
plt.xlabel('Time')
plt.ylabel('Distance')
plt.title('Distances for  small consecutive random changes,\n'+info)
plt.legend(loc='center left',bbox_to_anchor=(1,0.5))
plt.savefig(path2figs+figname, bbox_inches='tight')

plt.figure()
figname=type_g+'other_distances3_unscaled.pdf'
sb.set_context('paper',font_scale=1.3)
sb.set_style('white')
it_c=0
for column in ['Jaccard', 'Hamming', 'IM', 'HIM','alpha=0.5,order=5', 'alpha=0.9,order=3','f(l)=exp(-1.2l)']:
        plt.plot(dist[column],label=str(column),color=colors2[it_c])
        it_c+=1
plt.xlabel('Time')
plt.ylabel('Distance')
plt.title('Distances for  small consecutive random changes,\n'+info)
plt.legend(loc='center left',bbox_to_anchor=(1,0.5))
plt.savefig(path2figs+figname, bbox_inches='tight')