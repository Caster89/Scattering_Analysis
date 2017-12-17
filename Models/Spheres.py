from config import _AvlbUnits, _UnitsSymbols,_UnitsConv, _AvlbSASFit, _AvlbSASFitDic, _AvlbSASFitDicInv, _lmfitModels, _lmfitModelFunctions, _lmfitDistFunctions
import numpy as np
import scipy as sp
from scipy import signal
from scipy import interpolate
from scipy.integrate import dblquad, tplquad
import logging
import lmfit
from Form_Factors import*
from lmfit import minimize, Parameter, report_fit
try:
    import gmpy2
    from gmpy2 import mpz,mpq,mpfr,mpc
except:
    gmpy2 = None
    from decimal import Decimal
    print 'The Schultz fitting function is optimized by using the GMPY2 module to\
deal with the large numbers required. The module was not found, so Decimal will\
be used instead, but the calculations will be slower.'

class SingleGaussianSphereModel(lmfit.Model):
    def __init__(self, *args, **kwargs):
        super(SingleGaussianSphereModel, self).__init__(single_gauss_spheres,*args,**kwargs)


    def guess(self, data, **kwargs):
        pass



class DoubleGaussianSphereModel(lmfit.Model):
    def __init__(self, *args, **kwargs):
        super(DoubleGaussianSphereModel, self).__init__(double_gauss_spheres, *args, **kwargs)

    def guess(self, data, **kwargs):
        pass


class SingleSchultzSphereModel(lmfit.Model):
    def __init__(self, *args, **kwargs):
        super(SingleSchultzSphereModel, self).__init__(single_schultz_spheres, *args, **kwargs)

    def guess(self, data, **kwargs):
        pass


class DoubleSchultzSphereModel(lmfit.Model):
    def __init__(self, *args, **kwargs):
        super(DoubleSchultzSphereModel, self).__init__(double_schultz_spheres, *args, **kwargs)

    def guess(self, data, **kwargs):
        pass