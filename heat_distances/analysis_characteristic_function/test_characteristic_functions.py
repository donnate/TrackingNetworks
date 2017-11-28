# -*- coding: utf-8 -*-
"""
Created on Fri May 12 09:12:32 2017

@author: cdonnat
"""


import numpy as np
import networkx as nx 
import pandas as pd
import matplotlib.pyplot as plt
import sklearn as sk
import matplotlib.pyplot as plt
import sys
sys.path.append( '../../GSP/')
from SC import *
from distances_signature import *
from shapes.shapes import *
from heat_diffusion import *
from clustering_via_distances import *
from graph_tools import *
from stats.statistic_checks import *
from purity import *
from characteristic_functions import *
import igraph
from roleX import *

# 1- Start by defining our favorite regular structure
width_basis=36
basis_type="cycle"
#nb_shapes=7
#shape=["fan",8]
#G,colors=build_regular_structure(width_basis,basis_type, nb_shapes,shape, start=0,add_random_edges=0,plot=True,savefig=False)
#A=nx.adjacency_matrix(G).todense()
#nb_clust=len(np.unique(colors))
list_shapes=[["fan",8]]*3+[["star",4]]*3+[["house"]]*3+[["tree",3,2]]*3
G,colors_shape, plugins,colors=build_structure(width_basis,basis_type,list_shapes, start=0,add_random_edges=0,plot=False,savefig=False)
A=nx.adjacency_matrix(G).todense()
nb_clust=len(np.unique(colors))
taus=[1,2,5,10, 25, 50,100,200,300]
heat_print=heat_diffusion(G,taus,diff_type="immediate",type_graph="nx")

##### Start by plottting a few caracteristic function
mode=2
bunch=[0,1,2,10,45,67,90,100,111]
chi=plot_bunch_characteristic_functions(heat_print,mode, bunch)
##### featurize the heatprint
chi=featurize_characteristic_function(heat_print,t=range(1,300,1))
pca=sk.decomposition.PCA(n_components=2)
chi_new=pca.fit_transform(chi)
plot=True
savefig=True
annotate=[True,"color"]
if plot==True:
    sb.plt.figure()
    sb.set_style("white")
    cmap = plt.get_cmap('hot')
    plt.scatter(chi_new[:,0],chi_new[:,1],c=colors,cmap=cmap)
    if annotate[0]==True:
        if annotate[1]=="index":
            lab=range(chi_new.shape[0])
        else:
            lab=colors
            for label,c, x, y in zip(lab,colors, chi_new[:, 0], chi_new[:, 1]):
                #print label,x,y
                plt.annotate(label,xy=(x, y), xytext=(0, 0), textcoords='offset points')
        plt.colorbar()
        plt.title("PCA projection of the wavelet distribution, colored by status indicator")
        if savefig==True:
            plt.savefig("plots_paper/PCA_multiple_shapes_graph_aggregated_chi.png")
purity_chi=purity(D_chi,colors,16)
purity_RoleX=purity(D_roleX,colors,16)           
diff_purity=np.array(purity_chi)-np.array(purity_RoleX)
categories=[(k,v[0]) for k,v in colors_gen_map.iteritems()]
plt.figure()
colors_ref=["red","black","green","blue","salmon","pink","lightblue","orange", "yellow", "purple", "violet","tomato", "teal",\
            "grey","olive", "maroon","lime", "lavender","coral","aqua","sierra", "magenta","yellowgreen","darkgreen"]
for c in np.unique(colors):
    #categories=[(k,v[0]) for k,v in colors_gen_map.iteritems()]
    name=[i for i in range(len(categories)) if categories[i][1]==c]
    name=categories[name[0]][0]
    ind=[i for i in range(len(colors)) if colors[i]==c]
    selected_purity=diff_purity[ind,:]
    plt.plot(np.mean(selected_purity,0),c=colors_ref[c],label=name)
plt.title("Difference in Neighbor Purity levels for the different categories")
plt.legend(loc="upper left")
plt.savefig("plots_paper/ENRON_chi_vs_RoleX_purity_overall.png")
    
