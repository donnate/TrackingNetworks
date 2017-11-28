# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 15:41:49 2017

@author: cdonnat
"""
#### Assess results of the distance matrix via the clustering that it induces on the set of nodes:

import numpy as np
import networkx as nx 
import pandas as pd
import matplotlib.pyplot as plt
import sklearn as sk
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
from sklearn.cluster import DBSCAN
from sklearn.cluster import Birch
from sklearn.cluster import SpectralClustering
from sklearn.cluster import AffinityPropagation
from sklearn.cluster import AgglomerativeClustering
from SC import *
import math
import seaborn as sb


def clustering_results(D,colors, nb_clust=7,type_clust="tsne",type_label="timestamp",name_addendum='',plot=True,savefig=False,filefig="plots/clust.png"):
    ### Clusters the nodes according to a distance matrix D 
    ### INPUT:
    ### ======================================================
    ### D: distance matrix between nodes 
    ### colors: ground truth labels
    ### nb_clust: number of clusters that we want to recover
    ### type_clust: type of clustering : "tsne","kmeans" or "sc"
    ### type_label: type or the color labels: time stamp (if coing from a time series) or names ow
    ### plot, savefig,filefig: additional parameters (for plotting and saving plot results)
    ### OUTPUT:
    ### ======================================================
    ### labels_pred: labels_predicted
    ### hom=homogeneity of the clusters
    ### com=homogeneity of the clusters
    markers = {0:'*',1: '.', 2:',',3: 'o',4: 'v',5: '^',6: '<',7: '>',8: '1',9:'O',10: '+',11:'x',12:'D',13: '|',14: '_',15:'r'}

    N=D.shape[0]
    print "there are ", N, " nodes"
    if type_clust!="sc":
        
        if type_clust=="tsne":
            model = TSNE(n_components=min(7,nb_clust), random_state=0,metric="precomputed")
            np.set_printoptions(suppress=True)
            trans_data=model.fit_transform(D).T
            km=sk.cluster.KMeans(n_clusters=nb_clust)
            km.fit(trans_data.T)
            labels_pred=km.labels_
        else:
            if type_clust=="dbscan":
                model = DBSCAN(eps=0.2,metric="precomputed")
                labels_pred=model.fit_predict(D)
            elif type_clust=="affinity":
                model = AffinityPropagation(damping=.8,preference=-200)
                labels_pred=model.fit_predict(D)
                #labels_pred=[l for l in labels_pred1]
                labels_pred=np.ndarray.flatten(labels_pred)
                print labels_pred
                #reshape((-1,1))
            elif type_clust=="spectral_ref":
                model = SpectralClustering(n_clusters=nb_clust,affinity="precomputed")
                labels_pred=model.fit_predict(D)
            elif type_clust=="birch":
                model=Birch(n_clusters=nb_clust)
                labels_pred=model.fit_predict(D)
            model = TSNE(n_components=min(2,nb_clust), random_state=0,metric="precomputed")
            
            np.set_printoptions(suppress=True)
            trans_data=model.fit_transform(D).T
    
        #plt.subplots_adjust(bottom = 0.1)
        
       
        if plot==True:
            sb.set_context(font_scale=2.0)
            sb.set_style("ticks")
            sb.set_style({"xtick.direction":"in","ytick.direction":"in"})
            plt.figure()
            if type_label=="timestamp":
                plt.scatter(trans_data[0], trans_data[1], c=colors, cmap=plt.cm.rainbow)
                labels = [str(i) for i in range(N)]
            elif type_label=="names":
                #fig, ax = plt.subplots()
                #for xp, yp, c,m in zip(trans_data[0], trans_data[1],labels_pred, colors):
                    #plt.scatter([xp],[yp], color=c*1.0/nb_clust,marker=m)
                plt.scatter(trans_data[0], trans_data[1],c=labels_pred ,cmap=plt.cm.rainbow)
                labels = colors
            else:
                print "label type not recognized"
                labels = [str(i) for i in range(N)]
            
            for label,c, x, y in zip(labels,labels_pred, trans_data[0, :], trans_data[1, :]):
                #print label,x,y
                plt.annotate(label,xy=(x, y), xytext=(0, 0), textcoords='offset points')
                #plt.annotate(label,xy=(x, y))
            plt.title(type_clust+" 2-D projection of the data for method "+name_addendum)
            
            if savefig==True:
                plt.axis("off")
                plt.tight_layout()
                plt.savefig(filefig)
            plt.show()
    
    elif type_clust=="sc":
        U,labels_pred=SC(D, k=min(nb_clust,7),nb_clust=nb_clust,type_sim="sim")
        if plot==True:
            if type_label=="timestamp":
                labels = [str(i) for i in range(N)]
                plt.scatter(U[:,0], U[:,1], c=colors, cmap=plt.cm.rainbow)
            elif type_label=="names":
                labels = [colors[i] for i in range(N)]
                plt.scatter(U[:,0], U[:,1], c=colors, cmap=plt.cm.rainbow)
            else:
                print "label type not recognized"
                labels = [str(i) for i in range(N)]
                plt.scatter(U[:,0], U[:,1], c=colors, cmap=plt.cm.rainbow)
            for label, x, y in zip(labels, U[:, 0], U[:, 1]):
                print label,x,y
                plt.annotate(label,xy=(x, y), xytext=(0, 0), textcoords='offset points')
                #plt.annotate(label,xy=(x, y))
            if savefig==True:
                plt.savefig(filefig)
            plt.show()
    else:
        print "clustering method not recognized"
    hom=sk.metrics.homogeneity_score(colors, labels_pred) 
    comp=sk.metrics.completeness_score(colors, labels_pred) 
    return labels_pred,hom,comp


def clustering_representation(chi,colors, nb_clust=7,type_red="pca",type_clust="km",type_label="timestamp",name_addendum='',plot=True,savefig=False,filefig="plots/clust.png",verbose=False):
    ### Clusters the nodes according to a distance matrix D 
    ### INPUT:
    ### ======================================================
    ### D: distance matrix between nodes 
    ### colors: ground truth labels
    ### nb_clust: number of clusters that we want to recover
    ### type_clust: type of clustering : "tsne","kmeans" or "sc"
    ### type_label: type or the color labels: time stamp (if coing from a time series) or names ow
    ### plot, savefig,filefig: additional parameters (for plotting and saving plot results)
    ### OUTPUT:
    ### ======================================================
    ### labels_pred: labels_predicted
    ### hom=homogeneity of the clusters
    ### com=homogeneity of the clusters
    N=chi.shape[0]
    print "there are ", N, " nodes"
    if type_red=="tsne":
            model = TSNE(n_components=min(7,nb_clust), random_state=0)
            np.set_printoptions(suppress=True)
            trans_data=model.fit_transform(chi)
    else:
            pca=sk.decomposition.PCA(n_components=2)
            trans_data=pca.fit_transform(chi)
    if type_clust=="km":
        km=sk.cluster.KMeans(n_clusters=nb_clust)
        km.fit(trans_data)
        labels_pred=km.labels_
    elif type_clust=="dbscan":
        model = DBSCAN(eps=0.2)
        labels_pred=model.fit_predict(chi)
    elif type_clust=="affinity":
        model = AffinityPropagation(damping=.8,preference=-200)
        labels_pred=model.fit_predict(chi)
                #labels_pred=[l for l in labels_pred1]
        labels_pred=np.ndarray.flatten(labels_pred)

    elif type_clust=="spectral_ref":
        model = SpectralClustering(n_clusters=nb_clust,affinity="precomputed")
        labels_pred=model.fit_predict(trans_data)
    elif type_clust=="birch":
        model=Birch(n_clusters=nb_clust)
        labels_pred=model.fit_predict(chi)
    
        #plt.subplots_adjust(bottom = 0.1)
    if plot==True:
        plt.figure()
        sb.set_context(font_scale=3.0)
        sb.set_style("ticks")
        sb.set_style({"xtick.direction":"in","ytick.direction":"in"})
        if type_label=="timestamp":
            plt.scatter(trans_data[0], trans_data[1], c=colors, cmap=plt.cm.rainbow)
            labels = [str(i) for i in range(N)]
        elif type_label=="names":
                #fig, ax = plt.subplots()
                #for xp, yp, c,m in zip(trans_data[0], trans_data[1],labels_pred, colors):
                    #plt.scatter([xp],[yp], color=c*1.0/nb_clust,marker=m)
            if verbose: print len(labels_pred)
            if verbose: print trans_data[:,0]
            plt.scatter(trans_data[:,0], trans_data[:,1],c=labels_pred,cmap=plt.cm.rainbow,s=100)
            labels = colors
        else:
            print "label type not recognized"
            labels = [str(i) for i in range(N)]
            
        for label,c, x, y in zip(labels,labels_pred, trans_data[:, 0], trans_data[:, 1]):
            #print label,x,y
            plt.annotate(label,xy=(x, y), xytext=(0, 0), textcoords='offset points')
            #plt.annotate(label,xy=(x, y))
        plt.title(type_red+" 2-D projection of the data clustered via "+type_clust + '\n for method '+ name_addendum)
        if savefig==True:
            plt.axis("off")
            plt.tight_layout()
            plt.savefig(filefig)
        plt.show()
    hom=sk.metrics.homogeneity_score(colors, labels_pred) 
    comp=sk.metrics.completeness_score(colors, labels_pred) 
    return labels_pred,hom,comp

def get_nearest_neighbors(D,k,colors,plot=True,savefig=False,filefig="plots/nearest_neighbor_graph.png"):
    ### returns the nearest neighbors graph based on a distance matrix and evaluate the coherence of each node
        ### INPUT:
    ### ======================================================
    ### D: distance matrix between nodes 
    ### colors: ground truth labels
    ### k: number of  neighbors
    ### plot, savefig,filefig: additional parameters (for plotting and saving plot results)
    ### OUTPUT:
    ### ======================================================
    ### G: nearest neighbor graph
    ### coherence= coherence of the nearest neighbors (nb of different types)
    ### alignment= nb of the neighbors of the same type
    N=D.shape[0]
    G=nx.Graph()
    G.add_nodes_from(range(N))
    coherence=[None]*N
    alignment=[None]*N
    for i in range(N):
        nn=np.argsort(D[i,:])
        nn=nn[1:(k+1)]
        G.add_edges_from([(i,l) for l in nn])
        alignment[i]=len([j for j in nn if colors[j]==colors[i]])*1.0/k
    for i in range(N):
        nn=G.neighbors(i)
        nn.append(i)
        coherence[i]=1-len(np.unique([ colors[j] for j in nn ]))*1.0/len(nn)
    if plot==True:
        nx.draw_networkx(G,node_color=colors,cmap="hot") 
        if savefig==True:
                plt.savefig(filefig)
        plt.show()
            
                
                
    return G,coherence, alignment
