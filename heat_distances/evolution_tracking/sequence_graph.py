# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 12:39:10 2017

@author: cdonnat
"""
import numpy as np
import networkx as nx
import sys
from shapes.shapes import *
from heat_diffusion import *
from distances_between_graphs import *


def make_graph_evolve(G,Time=10, add=4, deleted=7):
    ### Makes a graph evolve through time by randomly removing and adding edges
    ### INPUT:
    ### ===============================================================
    ### G: a nx network
    ### Time: evolution horizon
    ### add: number of new edges created at each time shot
    ### deleted: number of edges deleted at each time shot
    ### OUTPUT:
    ### ===============================================================
    ### NetSeq: a dictionary of nx. graphs
    NetSeq={}
    NetSeq[0]=G
    edge_list=G.edges() 
    for t in range(1,Time+1):
        ## add and randomly remove edges of the graph
         
         #print "initial edge list length", len(edge_list)
         to_be_removed=[edge_list[e] for e in np.random.choice(range(len(edge_list)),deleted,replace=False)]
         #print "to be removed", 
         for item in to_be_removed:
             #print item
             edge_list.remove(item)
         for i in range(add):
             src,dest=np.random.choice(range(G.number_of_nodes()),2)
             already_there=((src,dest) in edge_list)
             while already_there==True:
                 src,dest=np.random.choice(range(G.number_of_nodes()),2)
                 already_there=((src,dest) in edge_list)
             #print src, dest
             edge_list.append((src,dest))
         #print len(edge_list)
         G_Temp=nx.Graph()
         G_Temp.add_nodes_from(G)
         G_Temp.add_edges_from(edge_list)
         NetSeq[t]=G_Temp
    return NetSeq
    
    
    
    
    
    
### define a sequence of evolving graphs
def sequence_evolving_graphs(Time=10,add=4,deleted=7):
    G,colors=build_regular_structure(30,"cycle", 20,["fan", 3], start=0,add_random_edges=0,plot=True,savefig=False)
    NetSeq={}
    NetSeq[0]=G
    for t in range(Time):
        ## add and randomly remove edges of the graph
         edge_list=G.edges()  
         G.remove_edges_from([edge_list[e] for e in np.random.choice(range(len(edge_list)),deleted)])
         for i in range(add):
             src,dest=np.random.choice(range(G.number_of_nodes()),2)
             G.add_edges_from([(src,dest)])
         NetSeq[t]=G
    return NetSeq


### Compute distances
def test(Time,add=2,deleted=1,mode=2):
    NetSeq=sequence_evolving_graphs(Time,add,deleted)
    heat_print={}
    taus=[1, 10, 25, 50]
    for t in range(Time):
        heat_print[t]=heat_diffusion(NetSeq[t],taus=taus,type_graph="nx")
    Distances2=pd.DataFrame(np.zeros((Time,Time)),index=range(Time), columns=range(Time))          
    for it_k in range(Time-1):
        print it_k
        HP1=heat_print[it_k]
        for it_kk in range(it_k+1,Time):
            HP2=heat_print[it_kk]
            _,Distances2.loc[it_k,it_kk]=distances_heatprints_mode(HP1,HP2,mode,type_comp="auc",plot=False,savefig=False)
            Distances2.loc[it_kk,it_k]=Distances2.loc[it_k,it_kk]
    G_time=nx.from_numpy_matrix(Distances2.as_matrix())
    T=nx.minimum_spanning_tree(G_time)
    plt.figure()
    pos=nx.fruchterman_reingold_layout(T) 
    nx.draw_networkx_nodes(T,pos,
                       node_color='r',
                       node_label=Distances2.index,
                       node_size=40,
                   alpha=0.8)
    nx.draw_networkx_edges(T,pos,width=1.0,alpha=0.5)
    nx.draw_networkx_labels(T,pos,{k:Distances2.index[k] for k in range(len( Distances2.index))},font_size=16)
    delta_dist=[Distances2.iloc[t,t+1] for t in range(Time-1)]
    plt.savefig("plots/dMinSpanTreeAdd"+str(add)+"Del"+str(deleted)+"Mode"+str(mode)+".png")
    
    plt.figure()
    plt.plot(range(Time-1),delta_dist)
    plt.title("Distance between consecutive graphs for mode "+str(mode))
    plt.savefig("plots/delta_distAdd"+str(add)+"Del"+str(deleted)+"Mode"+str(mode)+".png")
    return Distances2,T,delta_dist



    
    
        
        
        
