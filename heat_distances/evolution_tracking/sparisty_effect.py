# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 18:47:54 2017

@author: cdonnat
"""

ext='spar_'
import numpy as np
import networkx as nx
import sys
sys.path.append( '/Users/cdonnat/Dropbox/Distances/heat_distances/')
from shapes.shapes import *
from heat_diffusion import *
from distances_between_graphs import *
from graph_tools import *
from evolution_tracking.distances_evolution  import *

#name='smooth_rdm_30_0.4'

b=0.2
type_g='spar_rdm'
name='80_0.2_spar_'+str(b)+'distances.txt'
path2graph='/Users/cdonnat/Dropbox/Distances/tests_synthetic_data/generated_graphs/'
graph_name='80_0.2_spar_0.2.txt'
path2figs='/Users/cdonnat/Dropbox/Distances/tests_synthetic_data/figs/'
figname=type_g+'AUCdistances.pdf'
D=pd.DataFrame.from_csv(path2graph+graph_name,header=None,index_col=None,sep=',')
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



figname=type_g+'Chidistances'+type_g+'.pdf'
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
B=[0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6]

distances={tau:{} for tau in taus}
distances_chi={tau:{} for tau in taus}
distances_chi_agg={}
for b in B:
    name='80_'+str(b)+'_spar_'+str(b)+'distances.txt'
    path2graph='/Users/cdonnat/Dropbox/Distances/tests_synthetic_data/generated_graphs/'
    graph_name='80_'+str(b)+'_spar_'+str(b)+'.txt'
    path2figs='/Users/cdonnat/Dropbox/Distances/tests_synthetic_data/figs/'
    D=pd.DataFrame.from_csv(path2graph+graph_name,header=None,index_col=None,sep=',')
    N=int(np.sqrt(D.shape[1]))
    T=D.shape[0]
    NetSeq={t:nx.from_numpy_matrix(D.iloc[t,:].values.reshape((N,N))) for t in range(T)}
    taus=range(1,25,2)
    dist,dist_chi,dist_chi_agg,coeff,chi_agg=track_evolution(NetSeq,diffusion_type="heat",tau=taus)
    for tau in taus:
        distances[tau][b]=dist[tau]
        distances_chi[tau][b]=dist_chi[tau]
    distances_chi_agg[b]=dist_chi_agg
    
    
cmap=sb.plt.get_cmap('hot')
x=np.linspace(0,1,len(B))
colors={B[i]: cmap(x[i]) for i in range(len(B))}
#colors['chi aggr.']=cmap(x[len(B)])

for tau in taus:
    figname='spar_plot_AUC_'+str(tau)+'.pdf'
    fig, ax = plt.subplots()
    sb.set(style='ticks')
    sb.set_context("paper", font_scale=1.5)
    sb.set_style('white')
    #distances[tau],distances_chi[tau],coeff[tau],chi_agg[tau]=track_evolution(NetSeq,diffusion_type="heat",tau=tau)
    for b in B:    
        plt.plot(distances[tau][b],c=colors[b],label='s='+str(b))
    plt.title('Evolution of the AUC Distances between consecutive \n graphs for different sparsity levels s')
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
    
    
    
    figname='spar_plot_chi_'+str(tau)+'.pdf'
    fig, ax = plt.subplots()
    sb.set(style='ticks')
    sb.set_context("paper", font_scale=1.5)
    sb.set_style('white')
    #distances[tau],distances_chi[tau],coeff[tau],chi_agg[tau]=track_evolution(NetSeq,diffusion_type="heat",tau=tau)
    for b in B:    
        plt.plot(distances_chi[tau][b],c=colors[b],label='s='+str(b))
    plt.title('Evolution of the AUC Distances between consecutive \n graphs for different sparsity levels s')
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
    
    
figname='spar_plot_chi_agg_.pdf'
fig, ax = plt.subplots()
sb.set(style='ticks')
sb.set_context("paper", font_scale=1.5)
sb.set_style('white')
#distances[tau],distances_chi[tau],coeff[tau],chi_agg[tau]=track_evolution(NetSeq,diffusion_type="heat",tau=tau)
for b in B:    
    plt.plot(distances_chi_agg[b],c=colors[b],label='s='+str(b))
plt.title('Evolution of the  aggr. characteristic distances between consecutive \n graphs for different sparsity levels s')
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
    
    
sharpe_agg={b: np.std(distances_chi_agg[b])/np.mean(distances_chi_agg[b]) for b in B}
sharpe_chi={t:{b: np.std(distances_chi[t][b])/np.mean(distances_chi[t][b])for b in B} for t in taus} 
sharpe={t:{b:np.std(distances[t][b])/np.mean(distances[t][b]) for b in B} for t in taus}
for t in taus:
    fig,ax=plt.subplots()
    sb.set_style('white')
    plt.plot(B,[sharpe[t][b] for b in  B],label='AUC',c='black')
    plt.plot(B,[sharpe_chi[t][b] for b in  B],label='chi',color='red')
    #plt.scatter(taus[-1]+1,sharpe_agg,label='chi agg',color='blue',s=100)
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
    plt.savefig(path2figs+type_g+'sharpe_'+str(b)+'.pdf', bbox_inches='tight')

#figname='spar_plot_var_'+str(tau)+'.pdf'
#fig, ax = plt.subplots()
#sb.set(style='ticks')
#sb.set_context("paper", font_scale=1.5)
#sb.set_style('white')
##distances[tau],distances_chi[tau],coeff[tau],chi_agg[tau]=track_evolution(NetSeq,diffusion_type="heat",tau=tau)
#plt.plot([[distances[tau][b] for b in B])
#plt.title('Evolution of the AUC Distances between consecutive \n graphs for different sparsity levels s')
#plt.legend(loc='center left', bbox_to_anchor=(1.0,0.5))
#plt.xlabel('Time t')
#plt.ylabel(r'Distance between graphs $G_t$ and $G_{t+1}$')
#plt.tight_layout()
#ax.spines['left'].set_position('zero')
#ax.spines['right'].set_color('none')
#ax.spines['bottom'].set_position('zero')
#ax.spines['top'].set_color('none')
#ax.xaxis.set_ticks_position('bottom')
#ax.yaxis.set_ticks_position('left')
#plt.savefig(path2figs+figname, bbox_inches='tight') 
