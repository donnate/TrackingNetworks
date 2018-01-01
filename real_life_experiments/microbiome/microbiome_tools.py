# -*- coding: utf-8 -*-
"""
Created on Fri May 26 17:16:53 2017

@author: cdonnat
"""
import numpy as np
import pandas as pd 

def delta(A,thres,isAnamed=True):
    try:
        m=A.shape[1]
    except:
        m=1
    B=np.zeros((A.shape[0],m))
    for j in range(m):
            ind=[i for i in range(A.shape[0]) if np.abs(A.as_matrix()[i,j])>thres]
            B[ind,j]=1
    if isAnamed==True:
        B=pd.DataFrame(B,index=A.index)
        if m>1:
            B.columns=A.columns.values
    return B
    
def cooccurences(S,thres):
    gA=pd.DataFrame(np.zeros((S.shape[0],S.shape[0])),index=S.index)
    gA.columns=S.index
    S=delta(S,thres,isAnamed=True)
    for c in S.columns.values:
        gA+= np.matmul((np.reshape(S.loc[:,c],[-1,1])),np.reshape((S.loc[:,c]),[1,-1]))
    gA_compressed_ind=[u for u in S.index if np.sum(gA.loc[:,u])>0]
    gA_compressed=gA.loc[gA_compressed_ind,gA_compressed_ind]
    return gA_compressed
    
    