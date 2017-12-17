from config import _AvlbUnits, _UnitsSymbols,_UnitsConv, _AvlbSASFit, _AvlbSASFitDic, _AvlbSASFitDicInv, _lmfitModels, _lmfitModelFunctions, _lmfitDistFunctions
import numpy as np
import scipy as sp
from scipy import signal
from scipy import interpolate
from scipy.integrate import dblquad, tplquad
import logging
import lmfit
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


def single_gauss_distribution(x, R_av, sigma, I0):
    """single_gauss_distribution: returnd the PDF for a normal distribution
    given the necessary parameters.
        Args:
            x (:obj: numpy.array): the list of values for which the PDF has to be
                calculated
            R_av (int): the mean of the distribution
            sigma (int): the dispersion of the distribution
            I0 (int): the intensity of the distribution
        Returns
            The vector with the propability densities with the same size as x.
    """
    return 10**(I0)*np.exp(-(x-R_av)**2/(2*sigma**2))/np.sqrt(2*sigma**2*np.pi)

def double_gauss_distribution(x,R1_av, sigma1, R2_av, sigma2, I0, ratio):
    """double_gauss_distribution: returns the PDF for normal distribution
    given the necessary parameters.
        Args
            x (:obj: numpy.array): the list of values or which the PDF has to be
                calculated
            R_av1 (int): the average of the first distribution.
            sigm1 (int): the dispersion of the first distribution
            R_av2 (int): the average of the second distribution.
            Z2 (int): the dispersion of the second distribution
            I0 (int): the intensity of the distribution
            ratio (int): the ratio between the two distributions
        Returns
            The vector with the PDF function with the same size as x
    """
    return ratio*single_gauss_distribution(x, R1_av, sigma1, I0) + (1-ratio)*single_gauss_distribution(x, R2_av, sigma2, I0)
# z = 1 / sigma^2
# \text{SZ}(R,N,R_a,k) =  \frac{N}{R_a} \left(\frac{R}{R_a}\right)^{k-1} \frac{k^k\exp(-kR/R_a)}{\Gamma(k)}
#http://sasfit.ingobressler.net/manual/Schultz-Zimm


def single_schultz_distribution(x, R_av, Z, I0):
    """single_schultz_distribution: returns the PDF for the schultz function
    given the necessary parameters.
        Args
            x (:obj: numpy.array): the list of values or which the PDF has to be
                calculated
            R_av (int): the average of the distribution.
            Z (int): the dispersion of the distribution, as defined for a Flory
                Schultz distribution (Z = 1/sigma^2)
            I0 (int): the intensity of the distribution
        Returns
            The vector with the probabilty densities with the same size as x
    """
    if gmpy2 is None:
        x = np.array([Decimal(xx) for xx in x])
        R_av = Decimal(R_av)
        Z = Decimal(Z)
        shape=Z+Decimal(1.)
        I0 = Decimal(I0)
        gamma = Decimal(sp.special.gamma(float(shape)))
        #scale=1
        #gamma=bins**(shape-1)*(np.exp(-bins/scale) /(sp.special.gamma(shape)))
        returnVal = Decimal(10.)**I0*(shape/R_av)**shape * x**Z * np.exp(-(shape/R_av)*x) / gamma
    else:
        x = np.array([gmpy2.mpfr(xx) for xx in x])
        R_av = gmpy2.mpfr(R_av)
        Z = gmpy2.mpfr(Z)
        shape=Z+gmpy2.mpfr(1.)
        I0 = gmpy2.mpfr(I0)
        gamma = gmpy2.mpfr(sp.special.gamma(float(shape)))
        #scale=1
        #gamma=bins**(shape-1)*(np.exp(-bins/scale) /(sp.special.gamma(shape)))
        returnVal = 10.**I0*(shape/R_av)**shape * x**Z * np.array([gmpy2.exp(-(shape/R_av)*xx) for xx in x]) / gamma
    returnVal = np.array([float(rr) for rr in returnVal])
    #print 'Single Schultz returned a vector of {} with {}nan values'.format(returnVal.shape,np.sum(np.isnan(returnVal)))
    #print 'given R_av:{} Z:{} I0:{}'.format(R_av,Z,I0)
    return returnVal
    #10**(I0)*(shape**shape*(x/R_av)**Z*np.exp(-shape*x/R_av)/(R_av*sp.special.gamma(shape)))

def double_schultz_distribution(x, R1_av, Z1, R2_av, Z2, I0, ratio):
    """double_schultz_distribution: returns the PDF for the schultz function
    given the necessary parameters.
        Args
            x (:obj: numpy.array): the list of values or which the PDF has to be
                calculated
            R_av1 (int): the average of the first distribution.
            Z1 (int): the dispersion of the first distribution, as defined for a Flory
                Schultz distribution (Z = 1/sigma^2)
            R_av2 (int): the average of the second distribution.
            Z2 (int): the dispersion of the second distribution, as defined for a Flory
                Schultz distribution (Z = 1/sigma^2)
            I0 (int): the intensity of the distribution
            ratio (int): the ratio between the two distributions
        Returns
            The vector with the PDF function with the same size as x
    """
    #return I0

    returnVal =  ratio*single_schultz_distribution(x,R1_av,Z1,I0) + (1-ratio)*single_schultz_distribution(x,R2_av,Z2,I0)
    #print 'Single Schultz returned a vector of {} with {}nan values'.format(returnVal.shape,np.sum(np.isnan(returnVal)))

    return returnVal


