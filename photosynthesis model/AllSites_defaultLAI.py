"""
Date: 2024/01/29
Desciption: 
    Calculate Topt with a model including canopy
    Use site-specific LAI, otherwise, use default LAI, see ppt#5 p31
"""

import os
import datetime
import pandas as pd
import numpy as np
import const
from Topt_func_Canopy import Topt_4_Canopy

# need to change for every test
TEST_NAME = 'Canopy_defaultLAI/'

# mkdir to save the files
PATH_NAME = "/path2FLUXNET_output_Topt/" + TEST_NAME
SAVE_NAME = PATH_NAME + "%s.csv"
if not os.path.exists(PATH_NAME):
    os.makedirs(PATH_NAME)

#  variable used in FvCB model
# VPD: hPa, PA: kPa, but the model requires：kPa, hPa.
keys = ["TA_F", "VPD_F", "PA_F", "SW_IN_F", "CO2_F_MDS"]        
QC_keys = QC_keys = ["TA_F_QC", "VPD_F_QC", "PA_F_QC", "SW_IN_F_QC", "CO2_F_MDS_QC"]

# sites in south hemisphere
name_S = ['AU-DaP','AU-DaS','AU-How','AU-Stp','AU-Tum','BR-Sa1']
half_yr = datetime.timedelta(days=182)

# params for photosynthesis used
dft_params = {
        'aEV': 42600, 'bEV': 1.14,
        'aSV': 645.13, 'bSV': -0.38,
        'aSJ': 658.8, 'bSJ': -0.52, 'bSJ_home': -0.32,
        'Vcmax25': 55,
        'ARJV': 2.56, 'BRJV': -0.0202, 'BRJV_home': -0.0173
    }

# load site PFT and LAI
INFO_PATH = '/site_info_LAI_202401.xlsx'
info = pd.read_excel(INFO_PATH, index_col=0)

# loop for each site
Chara = [30,26]
PATH = ['/FLUXNET_subset/',
           '/AMERIFLUX_subset/']
SET = ['FLX','AMF']
    
for c,path,sett in zip(Chara, PATH, SET):
    print(c,path)
    allpath=os.listdir(path)
    allpath.sort()
    for ipath in allpath:
        if ipath[0:3]==sett:   # FLX/AMF
            allfile = os.listdir(path+ipath)
            for ifile in allfile:
                if len(ifile)>30 and ifile[c]=='H':    
                    SITE_NAME = ifile[4:10]
                    FILE_NAME = path+ipath+'/'+ifile
                    print(FILE_NAME)

                    # read csv and 
                    site = pd.read_csv(FILE_NAME,engine='python')   # HH fluxnet file
                    # add timestep
                    site["timestamp"] = pd.to_datetime(site["TIMESTAMP_START"],format ="%Y%m%d%H%M")

                    ## south hemisphere
                    if SITE_NAME in name_S:
                        site["timestamp"] = site["timestamp"]+half_yr

                     ## add time
                    site["year"] = site["timestamp"].dt.year
                    site["month"] = site["timestamp"].dt.month
                    site["yearmonth"] = site["month"]+site["year"]*100
                    site["doy"] = site["timestamp"].dt.dayofyear
                    site["ymd"] = site["yearmonth"]*1000+site["doy"]
                    site["ymd_hm"] = site["ymd"]*100+site["timestamp"].dt.hour+site["timestamp"].dt.minute/60.

                    ## check data quality
                    YMD_qc = site[QC_keys+["year","yearmonth","ymd", "ymd_hm"]]
                    YMD_qc.loc[YMD_qc["CO2_F_MDS_QC"]<0,["CO2_F_MDS_QC"]] = 3
                    YMD_qc["QC_tt"] = YMD_qc[QC_keys].max(axis=1)
                    YMD_qc.loc[YMD_qc["QC_tt"]<=1,["QC_tt"]] = 0
                    YMD_qc.loc[YMD_qc["QC_tt"]>1,["QC_tt"]] = 1

                    ## valid hour
                    valid_HM = YMD_qc[YMD_qc["QC_tt"]<=0]["ymd_hm"].tolist()   # QC_tt的可能取值：0，1（invalid）

                    ## agg from hm into yr
                    Y_qc = YMD_qc[["year","yearmonth","ymd","QC_tt"]].groupby(by="year").mean().reset_index()
                    valid_yr = Y_qc[Y_qc["QC_tt"]<=0.3]["year"].tolist()

                    if len(valid_yr)==0:
                        print("none")
                    else:
                        flx_keys = site.loc[site["ymd_hm"].isin(valid_HM),["year","month","yearmonth",
                                                                           "doy","ymd","timestamp","GPP_NT_VUT_REF"]+keys]

                        ## Thm
                        flx_dmax = flx_keys[["timestamp","year","yearmonth","doy","ymd","TA_F"]].groupby(["ymd"]).agg({
                            "TA_F":np.max,"year":np.mean, "yearmonth":np.mean,"timestamp":np.min,"doy":np.mean}).reset_index()

                        flx_ym = flx_dmax[flx_dmax["year"].isin(valid_yr)].groupby(["yearmonth"]).mean().reset_index()
                        flx_ym = flx_ym.groupby(["year"]).max().reset_index()
                        Thm = flx_ym["TA_F"].mean()

                        ## Tgs_30d
                        flx_dave = flx_keys[["timestamp","year","yearmonth","doy","ymd","TA_F"]].groupby(["ymd"]).agg({
                            "TA_F":np.mean,"year":np.mean, "yearmonth":np.mean, "doy":np.mean, "timestamp":np.min}).reset_index()
                        # average
                        flx_keys["TA_30d"] = np.nan
                        d30 = datetime.timedelta(days=30)
                        for i in range(len(flx_dmax)):
                            uptime = flx_dmax["timestamp"][i]
                            downtime = flx_dmax["timestamp"][i]-d30
                            tempTA = flx_dave[
                                np.logical_and(flx_dave["timestamp"]<uptime,flx_dave["timestamp"]>=downtime)]["TA_F"].mean()
                            flx_keys.loc[flx_keys["ymd"]==flx_dmax["ymd"][i],["TA_30d"]] = tempTA

                        ## 3 highest GPP months
                        GPPkeys = ["GPP_NT_VUT_REF","NEE_VUT_REF_QC"]
                        GPPqc_ymd = site[GPPkeys+["year","yearmonth","ymd"]]
                        GPPqc_ymd.loc[GPPqc_ymd["NEE_VUT_REF_QC"]<=1,["NEE_VUT_REF_QC"]] = 0
                        GPPqc_ymd.loc[GPPqc_ymd["NEE_VUT_REF_QC"]>=2,["NEE_VUT_REF_QC"]] = 1
                        # agg into daily
                        GPPqc_ymd = GPPqc_ymd.groupby(["ymd"]).agg({
                            "year":np.mean, "yearmonth":np.mean,
                            "GPP_NT_VUT_REF":np.max,"NEE_VUT_REF_QC":np.mean}).reset_index()

                        # agg into from hm to month, find GPP_valid_m
                        GPPqc_ym = GPPqc_ymd.groupby(["yearmonth"]).mean().reset_index()
                        GPPvalid_mo = GPPqc_ym[GPPqc_ym["NEE_VUT_REF_QC"]<=0.3]["yearmonth"].tolist()  

                        GPPqc_y = GPPqc_ym.loc[GPPqc_ym["yearmonth"].isin(GPPvalid_mo),
                                               ["yearmonth","year","GPP_NT_VUT_REF"]].groupby(by="year").apply(
                            lambda x:x.nlargest(3,"GPP_NT_VUT_REF"))
                        mark_mo = GPPqc_y[GPPqc_y["year"].isin(valid_yr)]["yearmonth"].tolist()

                        # subset to mark_mo
                        flx_trg = flx_keys.loc[flx_keys["yearmonth"].isin(mark_mo)]
                        flx_trg["Iabs_apx"] = flx_trg["SW_IN_F"]*2.3

                        # get LAI
                        if SITE_NAME in info.index:
                            LAI = info.loc[SITE_NAME,'Max. LAI']
                            if np.isnan(LAI):
                                LAI = const.PFT_LAI[info.loc[SITE_NAME,'vegtype']]
                            print('LAI is %r'%LAI)
                        else:
                            print('site not found')
                            continue
                        
                        # initialize & calculate An
                        this_site = Topt_4_Canopy(flx_trg, Thm, dft_params)
                        flx_trg["GPP_model"] = this_site.NetPhoto_Canopy(LAI)

                        # save as csv
                        flx_trg[["timestamp","year","yearmonth","doy",
                                  "ymd","TA_F","TA_30d","GPP_model","GPP_NT_VUT_REF"]].to_csv(SAVE_NAME%SITE_NAME)



