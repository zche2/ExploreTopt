#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Date: 2023/06/23
Description: 
    Fig 2b
    Topt derived from EC towers
    Tgs derived from EC towers using PELT method (Fang et al., sci data, 2023)
    each PFT calculate Topt~Tgs separately
"""
import numpy as np
from functools import reduce
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm

#====================== def func =========================
def OLS_(x,y,xname = 'Tgs_mean'):
    # linear regression
    x_add = sm.add_constant(x)
    model = sm.OLS(y, x_add)
    results = model.fit()
    coef_df = pd.DataFrame({"params": results.params,
                        "t": round(results.tvalues,3), 
                        "p-values": round(results.pvalues,3) 
                         })
    coef_df[['coef_0.025','coef_0.975']] = results.conf_int()
    # record result
    slope = coef_df.loc[xname,'params']
    slope_err = slope - coef_df.loc[xname,'coef_0.025']
    R2 = results.rsquared
    # signif codes:
    p = coef_df.loc[xname,'p-values']
    if p <= 0.001:
        signif = 'p < 0.001'
    elif p <= 0.01:
        signif = 'p < 0.01'
    elif p <= 0.05:
        signif = 'p < 0.05'
    else :
        signif = 'p = %.2f'%p
    # return slope, R2, signif
    return '%.2f'%slope+'$\pm$'+'%.2f'%slope_err, R2, signif
   
#============================ main ============================ 
if __name__ =='__main__':  

#============================ load files ============================ 
    Tgs_filename = '/Tgsmax_0621_FLUX_NS_PELT.csv'
    Topt_filename = '/Topt_0709_FLUX_check.csv'
    FLUXinfo_filename = '/FLUXNET-data-availability-20220921205702.xlsx'

    Tgs = pd.read_csv(Tgs_filename, index_col = 0)
    Topt = pd.read_csv(Topt_filename, index_col = 0)
    flx = pd.read_excel(FLUXinfo_filename,
                        sheet_name = "merged",index_col = 'site')
    
#===================== process before correlations ======================= 
    # combine some vegtypes with small sample size
    flx.loc[flx['vegtype']=='WSA','vegtype'] = 'SAV'
    flx.loc[flx['vegtype']=='OSH','vegtype'] = 'SHR'
    flx.loc[flx['vegtype']=='CSH','vegtype'] = 'SHR'
    print("total site yrs with Topt found:",np.nansum(Topt==Topt))
    
    # calculate average of each site
    Tgs[['Tgs_mean', 'Tgs_std']] = Tgs.agg((np.mean, np.std), axis=1)
    Topt[['Topt_mean', 'Topt_std']] = Topt.agg((np.mean, np.std), axis=1)
    
    # merge
    dfs = [Tgs[['Tgs_mean', 'Tgs_std']],Topt[['Topt_mean', 'Topt_std']],flx[["vegtype","lat", "lon"]]]
    df = reduce(lambda x, y: pd.merge(x, y, on="site", how="inner"), dfs)
    df = df.dropna(axis=0,how="any") 
 
#==================== loop & plot regression line ====================
    VEG_TYPE = ['ENF', 'EBF', 'DBF', 'MF', 'SHR', 'SAV', 'GRA', 'WET', 'CRO']
    # create a figure
    nrow = 3
    ncol = 3 
    fig, axs = plt.subplots(nrow, ncol, figsize=(10, 10), dpi=300,
                            layout='constrained',sharex=True, sharey=True)
    for nn, ax in enumerate(axs.flat):
        pft = VEG_TYPE[nn]
        temp = df[df["vegtype"]==pft]
        # N: sample saize
        N = len(temp)
        # add text
        slope, R2, signif = OLS_(temp['Tgs_mean'],temp['Topt_mean'])
        str = 'slope = %s\n'%slope+\
            '%s\n'%signif+\
            '$R^{2}$ = %.2f'%R2
        # rg_plot
        sns.regplot(x=temp['Tgs_mean'], y=temp['Topt_mean'], ci=95,ax=ax,
                scatter_kws={"color":"royalblue","alpha":.5,
                             "edgecolors":"none","s":80},
                line_kws={"color":"black","linewidth":2.5})
        
        # generic set
        ax.set(adjustable='box', aspect='equal')
        ax.set_title("%s (N=%d)"%(pft,N),fontdict={'fontsize':20})
        ax.text(13,5,str,fontsize=13)
        ax.set_xlim(3,35)
        ax.set_ylim(3,35)
        ax.set_xlabel('')
        ax.set_ylabel('')
        # axis
        ax.xaxis.set_tick_params(labelsize=18)
        ax.yaxis.set_tick_params(labelsize=18)
      
    
    fig.supxlabel("Mean T$_{gs\\ max}^{air}$ (˚C)", fontsize=25)
    fig.supylabel("T$_{opt}^{eco}$ (˚C)", fontsize=25)
    
    
    
    
    
