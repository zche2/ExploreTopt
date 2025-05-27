import json
import pandas as pd
import numpy as np
import os
import csv
from ToptCalculation import find_opt_T

with open("Json_file_to_modelling_outputs","r", encoding='UTF-8') as f:
    mydict = json.load(f)

path      = mydict["path"]
SAVE_NAME = mydict["SAVE_NAME"]

allpath=os.listdir(path)
allpath.sort()

# prepare a year list
all_yr = list(np.arange(1992,2022))
final = [["site"]+all_yr]

for ipath in allpath:
    SITE_NAME = ipath[0:6]
    FILE_NAME = path+'/'+ipath
    # read csv, already QC
    site = pd.read_csv(FILE_NAME, index_col=0) 
    # record valid yr
    valid_yr = site["year"].unique()
    print(valid_yr, type(valid_yr))
    # aggregate into ymd, cal Tdmax and GPP_acc
    site.loc[site["GPP_model"]<0,'GPP_model'] = np.nan
    site = site[["year","yearmonth","ymd","TA_F","GPP_model"]].groupby(["ymd"]).agg({
        "TA_F":np.max,"year":np.mean, "yearmonth":np.mean,"GPP_model":np.sum}).reset_index()

    # cal Topt
    if len(valid_yr)==0:
        Topt = np.zeros(30) - 9999
    else:
        print(SITE_NAME)
        Topt = np.zeros(len(all_yr)) - 9999
        for i,yr in enumerate(all_yr):
            if yr not in valid_yr:
                continue
            else:
                print(i,yr)
                subset = site.loc[site['year']==yr]
                Topt[i] = find_opt_T(subset['TA_F'], subset['GPP_model'],0, valid_yr, yr)
    
    Topt = Topt.tolist()
    final.append([SITE_NAME]+Topt)

# save as new csv
for i in range(len(final)):
    with open(SAVE_NAME,'a+', newline='') as f:      
        csv_writer = csv.writer(f)
        csv_writer.writerow(final[i])


