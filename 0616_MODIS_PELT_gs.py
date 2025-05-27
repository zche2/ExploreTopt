"""
date: 2023/06/16
description: 用Fang(2023)的方法计算生长季，写起来要简便一些，因为不再需要划分半球
"""

import numpy as np
import matplotlib.pyplot as plt
from dateutil import rrule
import datetime as dt

from scipy import interpolate    # cubic spline interpolation
from scipy.signal import savgol_filter    # savgol filter

from changepy import pelt    # Change Point Detection: PELT
from changepy.costs import normal_var
## see: https://github.com/ruipgil/changepy


## 定义求单个点的GS的函数：

def sos_eos(NIRv, dd, SG_filter=True):
    """
    input: NIRv
    output: series of SOS and EOS
    """
    ## process
    tck = interpolate.splrep(dd, NIRv)
    d_num = np.arange(1,5475.1,1)
    NIRv_interp = interpolate.splev(d_num, tck)
    if SG_filter:
        NIRv_SGfilter = savgol_filter(NIRv_interp, 23*7, 3)
    else:
        NIRv_SGfilter = NIRv_interp
    
    var = np.std(NIRv_SGfilter)
    index = pelt(normal_var(NIRv_SGfilter,var), len(NIRv_SGfilter), pen=0.5)
    
    ## mean value between the two change pts
    me = np.zeros(len(index))
    cpts = np.zeros(len(index)+1, dtype=int)   # cpts is day NO.
    cpts[:-1] = index
    cpts[-1] = len(tdmax)-1

    for i in range(len(me)):
        start = cpts[i]
        end = cpts[i+1]
        me[i] = np.nanmean(NIRv_SGfilter[start:end])  # 两个change pts之间的算数平均

    ## identifying potential bottoms
    # diff1
    diff1 = np.diff(me)
    # sign
    sign = np.sign(diff1)
    # diff2
    diff2 = np.diff(sign)
    # find index where diff2 = 2
    loca = np.argwhere(diff2==2)+1
    bottom_idx = cpts[loca]
    bottom_idx = np.append(bottom_idx, len(NIRv_SGfilter)-1)
    bottom_val = me[loca]
    bottom_val = np.append(bottom_val, me[-1])

    # find peak and rule out false GS
    peak_idx = np.zeros(len(bottom_idx)-1, dtype=int)
    peak_val = np.zeros(len(bottom_idx)-1)
    maxgpp = np.nanmax(me)   # max gpp of the whole time series
    mingpp = np.nanmin(me)   # min gpp of the whole ts
    overall_ampl = maxgpp-mingpp
    for i in range(len(peak_idx)):
        start = loca[i,0]
        if i==len(peak_idx)-1:
            end = len(me)-1
        else:
            end = loca[i+1,0]
        maxgpp_rg = np.nanmax(me[start:end])
        mingpp_rg = min(me[start], me[end])
        maxgpp_idx = np.nanargmax(me[start:end])
        if maxgpp_rg > 0.25*maxgpp and (maxgpp_rg-mingpp_rg) > 0.25*overall_ampl:
        #### in two case the peaks are found:
        #### 1) the cycle peak larger than 25% of the overall peak
        #### 2) the cycle amplitude larger than 25% of the...
            peak_idx[i] = cpts[maxgpp_idx+start]
            peak_val[i] = maxgpp_rg
        else:
            peak_idx[i] = -1       # 对于不合理的peak, 将其下标置于-1
            peak_val[i] = np.nan

    ## SOS and EOS
    sos = []
    eos = []
    for i in range(len(peak_idx)):
        peak_tmp_idx = peak_idx[i]
        if peak_tmp_idx==-1:
            continue
        peak_tmp_val = peak_val[i]
        # find SOS
        btm_tmp_idx = bottom_idx[i]
        btm_tmp_val = bottom_val[i]
        amplitude = btm_tmp_val + (peak_tmp_val-btm_tmp_val)*0.25
        ts = NIRv_SGfilter[btm_tmp_idx:peak_tmp_idx]
        sos_tmp = np.argwhere(ts>=amplitude)[0]+btm_tmp_idx
        sos.append(sos_tmp[0])
        # find eos
        btm_tmp_idx = bottom_idx[i+1]
        btm_tmp_val = bottom_val[i+1]
        amplitude = btm_tmp_val + (peak_tmp_val-btm_tmp_val)*0.25
        ts = NIRv_SGfilter[peak_tmp_idx:btm_tmp_idx]
        eos_tmp = np.argwhere(ts>=amplitude)[-1]+peak_tmp_idx
        eos.append(eos_tmp[0])

    return sos,eos

def cal_Tgs(Tdmax, NIRv, dd, SG_filter=True):
    """
    input: daily Tmax, NIRv
    output: Tgs
    """
    valid_t = np.nansum(Tdmax==Tdmax)
    valid_vi = np.nansum(NIRv==NIRv)
    if valid_t<0.5*len(Tdmax) or valid_vi<0.5*len(NIRv):
        return np.nan
    idx = NIRv==NIRv
    sos,eos = sos_eos(NIRv[idx], dd[idx], SG_filter=SG_filter)
    Tgs = np.nanmean([np.nanmean(Tdmax[sos[i]:eos[i]]) for i in range(len(sos))])
    return Tgs
    
    
if __name__ =='__main__':  
    ## MODIS的记录日期
    year = np.arange(2001,2015.1,1)
    date = ['01-01','01-17','02-02','02-18','03-06','03-22','04-07','04-23','05-09','05-25','06-10','06-26',
           '07-12','07-28','08-13','08-29','09-14','09-30','10-16','11-01','11-17','12-03','12-19']

    dd = []
    dstart = dt.datetime.strptime('2000-12-31', "%Y-%m-%d").date()

    for i in year:
        for j in date:
            temp = '%d'%i+'-'+j

            dtime = dt.datetime.strptime(temp, "%Y-%m-%d").date()
            temp_minus = rrule.rrule(rrule.DAILY, dtstart = dstart, until = dtime).count()
            dd.append(temp_minus - 1)

            
    
    ipath1 = '/scratch/hezichang/data/MODIS/'
    ipath2 = '/scratch/hezichang/data/CRUNCEP/'
    opath = '/scratch/hezichang/data/MODIS/'

    ndvi = np.load(ipath1+'NDVI_16d_halfdeg_2001_2015.npy')
    nir = np.load(ipath1+'NIR_16d_halfdeg_2001_2015.npy')
    NIRv = ndvi*nir
    tdmax = np.load(ipath2+'Tmax_d_halfdeg_2001_2015.npy')
    
    # mask ndvi<0.1
    me_ndvi = np.nanmean(ndvi, axis=0)
    invalid = me_ndvi<0.1
    print(invalid.shape)
    NIRv[:,invalid] = np.nan
    tdmax[:,invalid] = np.nan
    invalid = invalid[:300,:]
    
    # calculate Tgs
    Tgs = np.zeros((100,720))
    for i in range(200,300):
        for j in range(720):
            temp = cal_Tgs(tdmax[:,i,j], NIRv[:,i,j], np.array(dd))
            print('(%d, %d): %.1f'%(i,j,temp))
            Tgs[i,j] = temp
            
            
    '''
    Tgs = np.array([cal_Tgs(tdmax[:,i,j], NIRv[:,i,j], np.array(dd)) for i in range(300) for j in range(720)])
    Tgs = Tgs.reshape(300,720)
    
    
    Tgs[invalid] = np.nan
    '''
    print(np.nansum(Tgs==Tgs))
    np.save(file=opath+'Tgs_halfdeg_01-15_PELT_200_300.npy', arr=Tgs)
    
    
    
    
    
    
    
    
    
    

