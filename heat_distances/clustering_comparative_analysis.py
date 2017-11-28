# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 09:44:33 2017

@author: cdonnat
Functions for comparing the different structural-equivalent recovery techniques
"""
import numpy as np
import networkx as nx 
import pandas as pd
import matplotlib.pyplot as plt
import sklearn as sk
import matplotlib.pyplot as plt
import sys
import os
from SC import *
from distances_signature import *
from shapes.shapes import *
from heat_diffusion import *
from clustering_via_distances import *
import graph_tools 
from utils import *
#from stats.statistic_checks import *
from purity import *
from roleX import *
from characteristic_functions import *
from cluster_analysis import *
import sys
#sys.path.append( '/struc2vec_alg/src')


def clustering_comparative_analysis(G, colors, heat_print,nb_clust,dirpath='/Users/cdonnat/Desktop/structural_equivalents',graph_nodes_name='nodes.txt',graph_edges_name='edges.txt',agg_max=5,struc2vec_input='test_edges.txt',struc2vec_output='testy3.txt'):
    ''' Pipeline for comparing the performance of the  methods for characterization of the topology
    (i.e RoleX, struc2vec (to be implemented), chi with different scale or aggregation range)
    
    

    INPUT:
    -----------------------------------------------------------------------------------------------
    G:          the nx graph
    colors:     the true "topological" coloring of the graph
    heat_print: corresponding heat_print
    agg_max:    max aggregation of the characteristic function
    
    
    OUTPUT:
    -----------------------------------------------------------------------------------------------
    D:          dictionary of the distances between nodes (key= name of one method)
    chi:        dictionary of the characteristic functions
    Perf:       pandas DataFrame with all the performance metrics
      
    '''
    graph_list_rep=[[i,colors[i]] for i in range(nx.number_of_nodes(G))]
    np.savetxt(dirpath+'/struc2vec_alg/graph_synthetic_experiments/'+graph_nodes_name,graph_list_rep,fmt='%i')
    nx.write_edgelist(G,dirpath+'/struc2vec_alg/graph_synthetic_experiments/'+graph_edges_name,data=False)
    
    D={}
    perf={}
    labels_pred={}
    chi={}
    roleX_phi,D["RoleX"],labels_pred,_,perf["RoleX"]=cluster_analysis(G,colors,algo="RoleX",type_analysis="PCA",k=3,nb_clust=nb_clust)
    chi["RoleX"]=roleX_phi
    struc2vec_phi,D["struc2vec"],labels_pred2,_,perf["struc2vec"]=cluster_analysis(G,colors,algo="struc2vec",type_analysis="PCA",k=3,dirpath=dirpath,struc2vec_input=struc2vec_input,struc2vec_output=struc2vec_output,nb_clust=nb_clust)
    chi["struc2vec"]=struc2vec_phi
    for tau in range(18):
        chi["chi"+str(tau)],D["chi"+str(tau)],labels_pred,_,perf["chi"+str(tau)]=cluster_analysis(G,colors,algo="heat",type_analysis="PCA",k=3,nb_clust=nb_clust,heat_print=heat_print, tau=tau)
    
    chi_agg=featurize_characteristic_function(heat_print[:5],t=range(0,200,3))
    D["chi_agg1_"+str(agg_max)]=distance_nodes(chi_agg)  
    labels_pred,hom,comp=clustering_representation(chi_agg,colors,nb_clust,type_red="pca", type_label="names",name_addendum='heat aggregated ('+str(agg_max)+' first scale values)',plot=True,savefig=False)
    perf["chi_agg1_"+str(agg_max)]=np.zeros(3+15)
    perf["chi_agg1_"+str(agg_max)][0]=hom
    perf["chi_agg1_"+str(agg_max)][1]=comp
    perf["chi_agg1_"+str(agg_max)][3:]=np.mean(purity(D["chi_agg1_5"],colors,16),0) 
    
    chi_agg=featurize_characteristic_function(heat_print[[1,2,5,6,7,8,9]],t=range(0,200,3))
    D["chi_agg1_mid"]=distance_nodes(chi_agg)  
    labels_pred,hom,comp=clustering_representation(chi_agg,colors,nb_clust,type_red="pca", type_label="names",name_addendum='heat aggregated (mid scale values)',plot=True,savefig=False)
    perf["chi_agg1_mid"]=np.zeros(3+15)
    perf["chi_agg1_mid"][0]=hom
    perf["chi_agg1_mid"][1]=comp
    perf["chi_agg1_mid"][3:]=np.mean(purity(D["chi_agg1_mid"],colors,16),0) 
    Perf=pd.DataFrame(np.zeros((len(perf),18)),index=perf.keys())
    Perf.columns=["homogeneity","completeness","score"]+["purity"+str(k+1) for k in range(15)]
    for k in perf.keys():
        Perf.loc[k,:]=perf[k]
    return D, chi, Perf
    
