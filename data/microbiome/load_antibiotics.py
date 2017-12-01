# -*- coding: utf-8 -*-
"""
Created on Sun May 21 10:17:57 2017

@author: cdonnat
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sklearn as sk
import matplotlib.pyplot as plt
import sys
sys.path.append( '../../../GSP/')



##### load antibiotic data



it=0
for line in open('/Users/cdonnat/Dropbox/GSP/tests_real/microbiome_analysis/data/antibiotics_otu.txt','rb'):
    line=line.strip('\n')
    lst=line.split(" ")
    lst=[ind.strip('"') for ind in lst]
    #print lst
    if it==0:
        index=lst
        #print lst
        
        otu=pd.DataFrame(np.zeros((2582,162)))
        otu.columns=lst
        index=[]
    else:
        index.append(lst[0])
        otu.iloc[it-1,:]=[int(f) for f in lst[1:]]
    it+=1
otu.index=index
print("done otu")
### Shift data
for i in range(162):
    m=np.min(otu.iloc[i,:])
    print m
    otu.iloc[i,:]+=abs(m)

    

it=0
for line in open('/Users/cdonnat/Dropbox/GSP/tests_real/microbiome_analysis/data/antibiotics_taxa.txt','rb'):
    line=line.strip('\n')
    lst=line.split("\t")
    lst=[ind.strip('"') for ind in lst]
    if it==0:
        taxa_col=lst
        taxa=pd.DataFrame(np.zeros((1651,8)))
        taxa.columns=taxa_col
        index_taxa=[]
    else:
        index_taxa.append(lst[0])
        taxa.iloc[it-1,:]=lst[1:]
    it+=1
taxa.index=index_taxa


it=0
index_sample=[]
for line in open('/Users/cdonnat/Dropbox/GSP/tests_real/microbiome_analysis/data/antibiotics_samples.txt','rb'):
    line=line.strip('\n')
    lst=line.split("\t")
    lst=[ind.strip('"') for ind in lst]
    if it==0:
        sample_col=lst
        sample=pd.DataFrame(np.zeros((162,3)))
        print lst
        sample.columns=sample_col
    else:
        index_sample+=[lst[0]]
        sample.iloc[it-1,:]=lst[1:]
    it+=1  
sample.index=index_sample



