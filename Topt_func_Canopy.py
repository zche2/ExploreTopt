''' Photosynthesis Model'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from const import const


class Topt_4_Canopy():
    def __init__(self, df, MAT, params = None):
        # Ca, t2m, Iabs, pb, VPD, iv, TGROWTH, THOME of all timestamp
        self.df = df
        self.Ca = df['CO2_F_MDS']
        self.t2m = df['TA_F']
        self.I0 = df['Iabs_apx']
        self.pb = df['PA_F'] * 10  # unit: hPa
        self.VPD = df['VPD_F'] * .1  # unit: kPa
        self.Tgs = df['TA_30d']
        self.Thome = MAT  # a constant temperature
        self.params = params or {}
    
    def Arrhenius(self, E):
        trg = np.exp((self.t2m-const.tref)*E/((const.tref+const.tp_00)*const.R*(self.t2m+const.tp_00)))
        return trg
    
    def Arrhenius_modified_s(self, E, D, S):
        trg = np.exp((self.t2m-const.tref)*E/((const.tref+const.tp_00)*const.R*(self.t2m+const.tp_00)))
        fz = 1+np.exp(((const.tref+const.tp_00)*S-D)/((const.tref+const.tp_00)*const.R))
        fm = 1+np.exp(((self.t2m+const.tp_00)*S-D)/((self.t2m+const.tp_00)*const.R))
        return trg*fz/fm
    
    def Assimilation_curve(self, Iabs):
        '''
        Modelling net carbon assimilation with EC tower meteorological data as forcings
        This model was adapted from Yin (2009) and the Fortran version of ORCHIDEE
        Partially comparable to a big leaf model, somehow representing a leaf scale simulation (?)
        [KEY UPDATE] parameters representing thermal adjustment (Kumarathunge et al.,2018)
        p.s. C3 plants
        '''
        O = self.pb/const.pb_std*0.209*1e6     # converted to ubar, 0.209 molar fraction

        # leaf boundary layer conductance: gb
        gb_h2o = const.gb_ref*44.6*(const.tp_00/(self.t2m+const.tp_00))*(self.pb/const.pb_std)
        gb_co2 = gb_h2o/const.ratio_H2O_to_CO2

        # params using Arr_mod_s function
        E_Vcmax = self.params['aEV'] + self.params['bEV'] * self.Tgs * 1000      # Unit: KJ/mol
        S_Vcmax = self.params['aSV'] + self.params['bSV'] * self.Tgs        
        Vcmax = self.params['Vcmax25'] * self.Arrhenius_modified_s(
            E_Vcmax, const.D_Vcmax, S_Vcmax
        )   
        S_Jmax = self.params['aSJ'] + self.params['bSJ'] * self.Tgs + self.params['bSJ_home'] * self.Thome      
        JVr = self.params['ARJV'] + self.params['BRJV'] * self.Tgs + self.params['BRJV_home'] * self.Thome     
        Jmax = JVr * self.params['Vcmax25'] * self.Arrhenius_modified_s(
            const.E_Jmax, const.D_Jmax, S_Jmax
        )
        gm = const.gm25 * self.Arrhenius_modified_s(
            const.E_gm, const.D_gm, const.S_gm
        )

        # params using Arr function
        KmC = const.KmC25 * self.Arrhenius(const.E_Kmc)
        KmO = const.KmO25 * self.Arrhenius(const.E_Kmo)
        gamma_star = const.gamma_star25 * self.Arrhenius(const.E_gamma_star)
        fvpd = 1 / (1 / (const.a1 - const.b1 * self.VPD) - 1)
        Rd = 0.01 * self.params['Vcmax25'] * self.Arrhenius(const.E_Rd)

        prbl = (const.alpha_LL * Iabs + Jmax)**2 - 4*const.theta * Jmax * const.alpha_LL * Iabs
        J = (const.alpha_LL * Iabs + Jmax - prbl**.5) / (2*const.theta)

        for limit_photo in range(2):
            # controlled by Rubisco functioning (carbon assimilation)
            if limit_photo==0:
               x1 = Vcmax
               x2 = KmC*(1+O/KmO)       # It should be O not Oi (comment from Vuichard)
            # otherwise, by Rubisco regeneration (light limitation)
            else:
               x1 = J/4
               x2 = 2*gamma_star
               
            a = const.g0*(x2+gamma_star) + (const.g0/gm+fvpd)*(x1-Rd)
            b = self.Ca*(x1-Rd) - gamma_star*x1 - Rd*x2
            c = self.Ca + x2 + (1/gm+1/gb_co2)*(x1-Rd)
            d = x2 + gamma_star + (x1-Rd)/gm
            m = 1/gm + (const.g0/gm+fvpd)*(1/gm+1/gb_co2)
            
            p = -(d + (x1-Rd)/gm + a*(1/gm+1/gb_co2)
                  + c*(const.g0/gm+fvpd))/m
            q = (d*(x1-Rd) + a*c + b*(const.g0/gm+fvpd))/m
            r = -(a*b)/m
            
            QQ = (p**2 - 3*q)/9
            UU = (2*p**3 - 9*p*q + 27*r)/54
            
            idx_1 = (QQ>0) & (abs(UU/QQ**1.5)<=1)
            PSI = np.arccos(UU[idx_1]/QQ[idx_1]**1.5)
            A1_tmp = -1 * Rd
            A1_tmp[idx_1] = -2*QQ[idx_1]**0.5*np.cos(PSI/3) - p[idx_1]/3
            
            # if change result
            if limit_photo==0:
                A1 = A1_tmp
            else:
                idx_2 = A1_tmp<A1
                A1[idx_2] = A1_tmp[idx_2]

        return A1

    def NetPhoto_Canopy(self, LAI0):
        '''
        Calculate hourly An using Assimilation_curve
        Integrate along the canopy with given LAI and coefficients for PAR (K and rho)
        Using 5-point Gaussian Integration Method
        '''
        LAI = np.array([.0469101, .2307534, .500000000, .7692465, .9530899]) * LAI0
        WGT = np.array([.1184635, .2393144, .284444444, .2393144, .1184635])
        
        # light absorbed by each layer, I0 = self.I0
        Iabs_cnpy =  const.K_PAR * (1 - const.rho_PAR) * np.exp(-const.K_PAR * LAI)
        
        # net CO2 assimilation by every layer with LAI
        A1 = self.Assimilation_curve(self.I0 * Iabs_cnpy[0])
        A2 = self.Assimilation_curve(self.I0 * Iabs_cnpy[1])
        A3 = self.Assimilation_curve(self.I0 * Iabs_cnpy[2])
        A4 = self.Assimilation_curve(self.I0 * Iabs_cnpy[3])
        A5 = self.Assimilation_curve(self.I0 * Iabs_cnpy[4])

        # weighted average
        self.df['An'] = LAI0 * (WGT[0] * A1
                                 + WGT[1] * A2
                                 + WGT[2] * A3
                                 + WGT[3] * A4
                                 + WGT[4] * A5)
        
        return self.df['An']
    
    def NetPhoto_Canopy_desLight(self, LAI0, Iabs_):
        '''
        By given Iabs, calculate An
        '''
        LAI = np.array([.0469101, .2307534, .500000000, .7692465, .9530899]) * LAI0
        WGT = np.array([.1184635, .2393144, .284444444, .2393144, .1184635])
        
        # light absorbed by each layer, I0 = self.I0
        Iabs_cnpy =  const.K_PAR * (1 - const.rho_PAR) * np.exp(-const.K_PAR * LAI)
        
        # net CO2 assimilation by every layer with LAI
        Iabs = pd.Series(data = np.zeros(len(self.I0)) + Iabs_)
        A1 = self.Assimilation_curve(Iabs * Iabs_cnpy[0])
        A2 = self.Assimilation_curve(Iabs * Iabs_cnpy[1])
        A3 = self.Assimilation_curve(Iabs * Iabs_cnpy[2])
        A4 = self.Assimilation_curve(Iabs * Iabs_cnpy[3])
        A5 = self.Assimilation_curve(Iabs * Iabs_cnpy[4])

        # weighted average
        self.df['An_%d'%Iabs_] = LAI0 * (WGT[0] * A1
                                 + WGT[1] * A2
                                 + WGT[2] * A3
                                 + WGT[3] * A4
                                 + WGT[4] * A5)
        
        return self.df['An_%d'%Iabs_]
    
    def NetPhoto_Canopy_varyLight(self, LAI0, light_min, light_max, intervals):
        '''
        given a series of light intensity, calculate An separately.
        '''
        light_lev = np.linspace(light_min, light_max, num=intervals)
        for i in light_lev:
            self.NetPhoto_Canopy_desLight(LAI0, Iabs_=i)

            
