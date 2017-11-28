# -*- coding: utf-8 -*-
"""
Created on Sun Apr 30 11:19:13 2017

@author: cdonnat
"""

#### Compute tests and compares sequences of evolving graph sequences:
import numpy as np
import networkx as nx
import sys
#sys.path.append( '../../heat_distances/')
from shapes.shapes import *
from heat_diffusion import *
from distances_between_graphs import *
from graph_tools import *
from characteristic_functions import *
from clustering_via_distances import *

def track_evolution(NetSeq,diffusion_type="heat",**kwargs):
    Time=len(NetSeq)
    
    Distances_chi_agg=[None]*(Time -1)
    coeff=[None]*Time
    G=NetSeq[0]
    chi_agg={}
    
    if diffusion_type=="ad_hoc":
        ### check for arguments
        if "level" in kwargs.keys():
            lev=kwargs["level"]
        else:
            lev=3
        coeff[0]=compute_total_diffusion(nx.adjacency_matrix(G).todense(),lev,agg=False)
    elif diffusion_type=="heat":
        if "tau" in kwargs.keys():
            s=kwargs["tau"]
        else:
            s=range(1,11,2)
        if type(s)==int or type(s)==float:
            s=[s]
        print "s is ",s
        coeff[0]=heat_diffusion(G,taus=s,type_graph="nx")
        #coeff[0]=coeff[0].as_matrix()
        chi_agg[0]=featurize_characteristic_function(coeff[0],t=range(0,100,3))
    else:
        print "Diffusion type not recognized. Switching to heat"
        return track_evolution(NetSeq,diffusion_type="heat")
    Distances={ss:[None]*(Time -1) for ss in s} 
    Distances_chi={ss:[None]*(Time -1) for ss in s} 
    for i in range(Time-1):
        if diffusion_type=="ad_hoc":
            coeff[i+1]=compute_total_diffusion(nx.adjacency_matrix(NetSeq[i+1]).todense(),lev,agg=False)
        else:
            coeff[i+1]=heat_diffusion(NetSeq[i+1],taus=s,type_graph="nx")
            #coeff[i+1]=coeff[i+1].as_matrix()
        chi_points=range(0,100,3)
        chi_agg[i+1]=featurize_characteristic_function(coeff[i+1],t=chi_points)
        for ss in range(len(s)):
            _,Distances[s[ss]][i]=compute_gen_auc(coeff[i][ss],coeff[i+1][ss])
            a=ss*len(chi_points)
            b=(ss+1)*len(chi_points)
            _,Distances_chi[s[ss]][i]=compute_gen_auc(chi_agg[i][:,a:b],chi_agg[i+1][:,a:b],mode_diff="L2")
        _,Distances_chi_agg[i]=compute_gen_auc(chi_agg[i],chi_agg[i+1],mode_diff="L2")
        
    return Distances,Distances_chi, Distances_chi_agg,coeff,chi_agg
    
