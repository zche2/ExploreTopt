"""Extract Topt, Method 2.2.1"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math

def findbin(x):
    y = (math.ceil(x)+math.floor(x))/2.
    # y = round(x)
    return y

def running_mean3(ARRAY):
    L = len(ARRAY)
    TRG = np.zeros(L)*np.nan
    for i in range(L):
        if i==0:
            TRG[i] = np.nanmean(ARRAY[i:i+2])
        elif i==L-1:
            TRG[i] = np.nanmean(ARRAY[i-1:i+1])
        else:
            TRG[i] = np.nanmean(ARRAY[i-1:i+2])
    return TRG

def find90(INPUT):
    ARRAY = np.array(INPUT)
    if np.nansum(ARRAY==ARRAY)==0:
        return np.nan
    L = len(ARRAY)
    ARRAY.sort()
    ORDER = np.arange(0,L,1)+0.5
    ORDER = 100*ORDER/L
    # print('order is:',ORDER)
    FIND = np.argwhere(ORDER==90)
    if FIND.shape[0]!=0:
        i = FIND[0,0]
        return ARRAY[i]
    FIND = np.argwhere(ORDER<90)
    i = FIND[-1,0]
    if i==L-1:
        return ARRAY[-1]
    j = i+1
    X = ARRAY[i]+(90-ORDER[i])*(ARRAY[j]-ARRAY[i])/(ORDER[j]-ORDER[i])
    return X

def find_opt_T(T,VI,flag=0,VALID_YR=[0], YR=0):
    ## if Topt can be derived
    if YR not in VALID_YR:
        return np.nan
    ## choose valid data
    if np.nansum(VI==VI)==np.nansum(VI<=0):     #for FLUXNET gpp！
        return np.nan
    elif np.nansum(T==T)==0:
        return np.nan
    else:
        valid = np.logical_and(T==T, VI==VI)
        if ~np.any(valid):
            return np.nan
        T_tmp = T[valid]
        VI_tmp = VI[valid]

        tt = [findbin(i) for i in T_tmp]
        t_array = np.unique(tt)
        t_array.sort()
        L = len(t_array)
        t_bin = np.zeros(L)*np.nan
        vi_bin = np.zeros(L)*np.nan
        for i in range(L):
            valid = tt==t_array[i]
            t_bin[i] = np.nanmean(T_tmp[valid])
            vi90 = np.nanmax(VI_tmp[valid])  # for modelled GPP
            vi_bin[i] = vi90
        # 12月11日更改
        if L>5:
            trg = running_mean3(vi_bin)
        else:
            trg = vi_bin
        index = np.nanargmax(trg)

        if flag==1:
            fig=plt.figure(figsize=(5,5),dpi=150)
            plt.scatter(T_tmp, VI_tmp, s = 30, c='silver')
            plt.plot(t_bin, trg,'k-',linewidth=3)
            plt.scatter(t_bin, vi_bin,c='coral',s = 60, zorder=3)
            global SITE_NAME
            plt.title("%s %d"%(SITE_NAME,YR), fontsize=20)
            plt.xticks(fontsize= 15)
            plt.yticks(fontsize= 15)
            ax=plt.gca()
            ax.spines['bottom'].set_linewidth(1.5)
            ax.spines['left'].set_linewidth(1.5)
            ax.spines['right'].set_linewidth(1.5)
            ax.spines['top'].set_linewidth(1.5)
            plt.show()
            # fig.savefig('path.png'%(SITE_NAME,pic_tt))#;plt.close(fig)

        # if index == len(trg)-1 or index==0:
        if index == len(trg)-1:
            return np.nan
        else:
            return t_bin[index].round(1)

