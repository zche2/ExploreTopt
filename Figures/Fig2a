#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Date: 2023/06/23
Description: 
    Fig 2a
    Topt derived from EC towers
    Tgs derived from EC towers using PELT method (Fang et al., sci data, 2023)
"""

import pandas as pd
import numpy as np
from functools import reduce
import matplotlib.pyplot as plt
import statsmodels.api as sm

#============================ load files ============================ 
Tgs_filename = '/Tgsmax_0621_FLUX_NS_PELT.csv'
Topt_filename = '/=Topt_0709_FLUX_check.csv'
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
df = df.dropna(axis='index',how="any") 

# print how many sites could be used for analysis
print("sites considered", len(df))

#===================== plot scatter ======================= 
VEG_TYPE = ['ENF', 'EBF', 'DBF', 'MF', 'SHR', 'SAV', 'GRA', 'WET', 'CRO']
MARK_TYPE = ['o', 's', '^', 'D', 'o', 'v', 's', 'D', '^']
COLOR_TYPE = ['forestgreen', 'limegreen', 'springgreen', 'orange', 'khaki', 'lightsalmon', 'palevioletred', 'darkblue','brown']
fig, ax = plt.subplots(figsize=(10, 10),dpi=300)
# plot scatter for each pft
for i in range(9):
    temp = df[df["vegtype"]==VEG_TYPE[i]]
    ax.errorbar(temp["Tgs_mean"],temp["Topt_mean"],yerr=temp["Topt_std"],xerr=temp["Tgs_std"],
                ecolor='k',elinewidth=0.5,marker=MARK_TYPE[i],
                mfc=COLOR_TYPE[i],mec='k',mew=1,ms=12,alpha=0.6,
                capsize=3,capthick=2.5,linestyle="none",label=VEG_TYPE[i])
# plot 1:1 line
flag = np.arange(df['Tgs_mean'].min(),df['Tgs_mean'].max()+5,5)
ax.plot([0,50],[0,50],linestyle=':',color='silver',linewidth=1.5) 
# sets
ax.set_xlim((2,39))
ax.set_ylim((2,39))   
ax.set_xlabel("Mean T$_{gs\\ max}^{air}$ (˚C)", fontsize=25)
ax.set_ylabel("T$_{opt}^{eco}$ (˚C)", fontsize=25)
ax.xaxis.set_tick_params(labelsize=18)
ax.yaxis.set_tick_params(labelsize=18)
ax.legend(loc='upper left', ncol=2,
          frameon=False, markerscale=0.9,
          fontsize=15)
axe=plt.gca();#获得坐标轴的句柄
axe.spines['bottom'].set_linewidth(1.5);###设置底部坐标轴的粗细
axe.spines['left'].set_linewidth(1.5);####设置左边坐标轴的粗细
ax.spines['right'].set_visible(False);###设置右边坐标轴的粗细
ax.spines['top'].set_visible(False);####设置上部坐标轴的粗细

#================= plot regression line =================== 
X = sm.add_constant(df['Tgs_mean'])
model = sm.OLS(df['Topt_mean'], X)
results = model.fit()
R2 = round(results.rsquared,2)
coef_df = pd.DataFrame({"params": results.params,   # coefficients
                        "t": round(results.tvalues,3),       # t-value
                        "p-values": round(results.pvalues,4) # p-value
                         })
coef_df[['coef_0.025','coef_0.975']] = results.conf_int() # confidence interval
print(coef_df)

# predict Y
Topt_predict = results.fittedvalues
Y_predict = results.predict(sm.add_constant(flag))
# RMSE
rmse = ((df['Topt_mean']-Topt_predict)**2).mean()**0.5
print('R2=%r\nrmse=%.2f'%(R2,rmse))
ax.plot(flag, Y_predict,
        linestyle='--', color='k',
        linewidth=3, zorder=3) 
# add text
slope = coef_df.loc['Tgs_mean','params']
slope_err = slope - coef_df.loc['Tgs_mean','coef_0.025']
str = 'slope = %.2f $\pm$ %.2f\n'%(slope, slope_err)+\
    '$R^{2}$ = %.2f\n'%R2+\
    'RMSE = %.2f\n'%(rmse)
ax.text(22,5,str,fontsize = 24)
fig.show()

