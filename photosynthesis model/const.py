class _const:
    class ConstError(TypeError): pass
    class ConstCaseError(ConstError): pass

    def __setattr__(self, name, value):
        if name in self.__dict__:
            raise self.ConstError("can't change const %s" % name)
        # if not name.isupper():
        #     raise self.ConstCaseError('const name "%s" is not all uppercase' % name)
        self.__dict__[name] = value

const = _const()
#
const.R = 8.314463
const.gb_ref = 1./25. 
const.pb_std = 1013
const.tp_00 = 273.15
const.ratio_H2O_to_CO2 = 1.6
#
const.E_Kmc = 79430.
const.E_Kmo = 36380.
const.E_Sco = -24460.
const.E_gamma_star = 37830.
const.E_Jmax = 40710.
const.D_Vcmax = 200000.
const.D_Jmax = 200000.
const.E_gm = 49600.
const.S_gm = 1400.
const.D_gm = 437400.
const.E_Rd = 46390.
#
const.KmC25 = 404.9
const.KmO25 = 278400
const.Sco25 = 2800
const.gm25 = 0.4  #### mol m-2 s-1 bar-1
const.gamma_star25 = 42.75
const.a1 = 0.85
const.b1 = .14
const.g0 = 0.00625
const.h_protons = 4.
const.theta = 0.85
const.alpha_LL = 0.24
# tref
const.tref = 25
# light
const.K_PAR = 0.715   # extinction coefficient for PAR (sphere)
const.rho_PAR = 0.057   # canopy reflection coefficient for PAR (sphere)
# PFT-specific LAI
PFT_LAI = {'GRA': 4, 'CRO': 4,
           'CSH': 4, 'OSH': 4,
           'SAV': 4, 'WSA': 4,
           'ENF': 5, 'EBF': 7,
           'DBF': 5, 'DNF': 5,
           'MF':  5, 'WET': 4,
           }
