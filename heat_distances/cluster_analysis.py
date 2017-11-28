# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 23:15:24 2017

@author: cdonnat

Set of functions to perform the analysis and compare the different results
"""
import numpy as np
import networkx as nx 
import pandas as pd
import matplotlib.pyplot as plt
import sklearn as sk
import matplotlib.pyplot as plt
import os
from SC import *
from distances_signature import *
from shapes.shapes import *
from heat_diffusion import *
from clustering_via_distances import *
from graph_tools import *
#from stats.statistic_checks import *
from purity import *
from characteristic_functions import *
from roleX import *
#from struc2vec_alg.src.main import *




def read_struc2vec_output(path_to_struct2vec,skip_header=True):
    ''' reads the struc2vec output to  convert it into the proper feature format
    INPUT
    -----------------------------------------------------------------------
    path_to_struct2vec :  path to the file containing all the embeddings (string)
    skip_header        :  boolean: should the first line of th struc2vec embedding file  be skipped?
    
    OUTPUT
    -----------------------------------------------------------------------
    feat_str2vec       :  the corresponding feature matrix (row: feature vector for node i)
    '''
    handle_map=open(path_to_struct2vec);
    feat_str2vec={}
    itt=0
    for line in handle_map:
        line=line.strip("\n")
        lst=line.split(" ")
        if itt==0:
            if not skip_header:
                N=int(lst[0])
                H=int(lst[1])
                feat_str2vec=np.zeros((N,H-1))
        elif itt==1:
            if skip_header:
                N=int(lst[0])
                H=int(lst[1])
                feat_str2vec=np.zeros((N,H-1))
            else:
                ind=lst[0]
            feat_str2vec[ind,:]=[float(lst[i]) for i in range(1,H)]    
            
        else:
            ind=lst[0]
            feat_str2vec[ind,:]=[float(lst[i]) for i in range(1,H)]
        itt+=1
    return feat_str2vec


def cluster_analysis(G,colors,algo,type_analysis="PCA",k=3,**kwargs):
    ''' function for computing the feature vectors, clustering nodes by their topological status and analyzing the induced clusters
        INPUT
        -----------------------------------------------------------------------
        G               : nx Graph
        colors          : colors of the nodes/ true labelling of the nodes (list)
        algo            : what algorithm should we test ("RoleX",  "heat" or "struc2vec")
        type_analysis   : what analysis should we perform (PCA (default), or clust)
        k               : number of neighbors to consider when computing the purity
        dirpath         : path to directory   
        graph_nodes_name: name of the file in which to write the list of nodes and their colors
        graph_edges_name: name of the file in which to write the edge list
        kwargs          : additional arguments (such as the number of clusters)
        
        OUTPUT
        -----------------------------------------------------------------------
        chi             : feature vector
        labels_pred     : outcome of the clustering
        neighbors       : dictionary of each node's k nearest topological neighbor
        perf            : score measuring the clusters's performance: average homogeneity, completeness, success rate and purity
    '''
    
    
    print "k=",k
    Adj_graph=nx.adjacency_matrix(G).todense()
    N=nx.number_of_nodes(G)
    
    if algo=="RoleX":
        Gi=igraph.Graph.Adjacency((Adj_graph > 0).tolist())
        test_RoleX=extract_rolx_roles(Gi,roles=kwargs["nb_clust"])
        chi=test_RoleX[0]
    elif algo=="heat":
        chi=featurize_characteristic_function_selected_mode(kwargs['heat_print'],kwargs['tau'],t=range(0,150,3))
    elif algo=="struc2vec":
        cmd='cd '+kwargs['dirpath']+'/struc2vec_alg; python src/main.py --input '+'graph_synthetic_experiments/'+kwargs["struc2vec_input"]+' --output '+'emb/'+kwargs["struc2vec_output"]
        os.system(cmd)        
        chi=read_struc2vec_output(kwargs['dirpath']+'/struc2vec_alg/emb/'+kwargs["struc2vec_output"],skip_header=False)
    else:
        print "algo type not recognized. Defaulting to heat."
        return cluster_rep(G,algo="heat")
    
        
    D=distance_nodes(chi)
    #print "D",D
    pur=purity(D,colors,16)
    #print "pur", pur
    if type_analysis=="PCA":
        labels_pred,hom,comp=clustering_representation(chi,colors,nb_clust=kwargs["nb_clust"],type_red="pca", name_addendum=algo,type_label="names",plot=True,savefig=False)
    elif type_analysis=="clust":
        labels_pred,hom,comp=clustering_results(D,colors,nb_clust=kwargs["nb_clust"],type_clust=kwargs["type_clust"],name_addendum=algo, type_label="names",savefig=True)
    else:
        print "type analysis not recognized."
        return np.nan
    
    perf=np.zeros(3+15)
    perf[0]=hom
    perf[1]=comp    
    score=0
    for n in range(nx.number_of_nodes(G)):
        neighbors=np.argsort(D[n,:]).tolist()
        neighbors=neighbors[1:(k+1)]
        success_rate=len([nn for nn in neighbors if colors[nn]==colors[n] ])*1.0/k
        score+=success_rate/N
    perf[2]=score
    perf[3:]=np.mean(pur,0)
    return chi,D,labels_pred,neighbors,np.array(perf)
