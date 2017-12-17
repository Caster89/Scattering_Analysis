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

class MonodisperseCubeModel(lmfit.Model):
    def __init__(self, *args, **kwargs):
        super(MonodisperseCubeModel, self).__init__(monodisperse_cube,*args,**kwargs)


    def guess(self, data, **kwargs):
        pass

class SingleGaussianCubeModel(lmfit.Model):
    def __init__(self, *args, **kwargs):
        super(SingleGaussCubeModel, self).__init__(single_gaussian_cube,*args,**kwargs)


    def guess(self, data, **kwargs):
        pass