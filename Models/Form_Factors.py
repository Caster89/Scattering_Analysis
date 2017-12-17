from config import _AvlbUnits, _UnitsSymbols,_UnitsConv, _AvlbSASFit, _AvlbSASFitDic, _AvlbSASFitDicInv, _lmfitModels, _lmfitModelFunctions, _lmfitDistFunctions
import numpy as np
import scipy as sp
from scipy import signal
from scipy import interpolate
from scipy.integrate import dblquad, tplquad
import logging
import lmfit
from Distributions import *
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

"""
The distriutions presented here can be found in Polydispersity analysis of scattering
data from self-assembled systems. Phy. Rev. A,45, 2428-2438.
DOI: 10.1103/PhysRevA.45.2428
"""
def single_gauss_spheres(q,R_av = 1,sigma = 1,I0 = 1,bckg=0):
    """sing_gauss_spheres: calculates the scattering pattern of an assembly of
    spheres which have a Gaussian number density size distribution.
        Args
            q (numpy.array): the array containg the list of q-values for which to
                calculate the scattering
            R_av (int): the mean of the size distribution. Defaults to 1
            sigma (int): the dispersion of the distribution. Defaults to 1
            I0 (int): the prefactor which includes information on the scattering
                length density (SLD) and the concentration of particles. Defaults
                to 1
            bckg (int): the background value to use in case the background is
                not perfectly subtracted. Defaults to 0.
        Returns
            the scattering curve which has the same size as q
    """
    P_q=(4.*np.pi/q**3.)**2.*((1.+q**2.*(R_av**2.+sigma**2.))/2.+\
        (((-1.+q**2*R_av**2)/2.-3./2.*q**2.*sigma**2.-2.*q**4.*sigma**4.)*\
        np.cos(2.*q*R_av)-q*R_av*(1.+2.*q**2.*sigma**2.)*\
        np.sin(2.*q*R_av))*np.exp(-2.*q**2.*sigma**2))

    return np.array(10**(I0)*P_q+bckg)


def double_gauss_spheres(q,R1_av = 1,sigma1 = 1, R2_av = 1, sigma2 = 1, I0 = 1,ratio=0.5, bckg = 0):
    """double_gauss_spheres: calculates the scattering pattern of an assembly of
    spheres which have a bimodal Gaussian size distribution.
        Args
            q (numpy.array): the array containg the list of q-values for which to
                calculate the scattering
            R_av1 (int): the mean of the size distribution of the first
                peak. Defaults to 1
            sigma1 (int): the dispersion of the first peak. Defaults to 1
            R_av2 (int): the mean of the size distribution of the second
                peak. Defaults to 1
            sigma2 (int): the dispersion of the second peak. Defaults to 1
            I0 (int): the prefactor which includes information on the scattering
                length density (SLD) and the concentration of particles. Defaults
                to 1
            ratio (int): the ratio between the first and the second peak. Defaults
                to 0.5
            bckg (int): the background value to use in case the background is
                not perfectly subtracted. Defaults to 0.
        Returns
            the scattering curve which has the same size as q
    """
    return np.array(ratio*single_gauss_spheres(q,R1_av, sigma1,I0,0)+(1-ratio)*single_gauss_spheres(q,R2_av, sigma2,I0,0)+bckg)


def single_schultz_spheres(q, R_av = 1, Z = 50, I0 = 1, bckg = 0 ):
    """sing_schultz_spheres: calculates the scattering pattern of an assembly of
    spheres which have a Schultz-Zimm size distribution. Devimal is used to
    ensure that the for vey monodisperse distributions (Z>171) the values are not
    rounded off to inf. The integrated function is taken from 'Analysis of small
    angle neutron scattering spectra from pplydisperse interacting colloids',
    DOI: 10.1063/1.446055
    sigma = R_av/(Z+1)^0.5
    ?The Z parameter is defined as z = 1 / sigma^2.?
    The definition was taken from:
    ttp://sasfit.ingobressler.net/manual/Schultz-Zimm
    """
    if gmpy2 is None:
        aD = np.array([Decimal((Z+1.)/(qq*R_av)) for qq in q])
        """
        numpy trigonometric functions do not support Decimal, therefore the
        numpy array is created on the spot using float numbers and transforming
        them to Decimal after the calculation
        """
        a = (Z+1.)/(q*R_av)

        p1 = Decimal(8. * np.pi**2 * R_av**6 * (Z+1)**(-6)) * aD**Decimal(Z+7.)
        G11 = aD**Decimal(-(Z+1.)) - (Decimal(4.)+aD**2)**(Decimal(-(Z+1.)/2.)) *\
        np.array([Decimal(np.cos((Z+1) * np.arctan(2./aa))) for aa in a])

        G12 = Decimal((Z+2.)*(Z+1.)) * (aD**Decimal(-(Z+3.)) + (Decimal(4.) + aD**Decimal(2.))**Decimal(-(Z+3.)/2.) *\
        np.array([Decimal(np.cos((Z+3)*np.arctan(2./aa))) for aa in a]))

        G13 = Decimal(2.*(Z+1.)) * (Decimal(4.) + aD**2.)**Decimal(-(Z+2.)/2) *\
        np.array([Decimal(np.sin((Z+2.)*np.arctan(2./aa))) for aa in a])
        G1 = G11+G12-G13
        returnVal = Decimal(10**I0)*p1*G1+Decimal(bckg)
    else:
        a = np.array([mpfr((Z+1.)/(qq*R_av)) for qq in q])
        a2 = a**2
        a2_1 = (mpfr(4.)+a2)
        R_av = mpfr(R_av)
        Z = mpfr(Z)
        I0 = mpfr(I0)
        bckg = mpfr(bckg)
        """
        numpy trigonometric functions do not support Decimal, therefore the
        numpy array is created on the spot using float numbers and transforming
        them to Decimal after the calculation
        """
        p1 = 8. * np.pi**2 * R_av**6 * (Z+1)**(-6) * a**(Z+7.)
        #G11 = a**-(Z+1.) - (4.+a**2)**(-(Z+1.)/2.) *\
        #np.array([gmpy2.cos((Z+1) * gmpy2.atan(2./aa)) for aa in a])
        G11 = a**-(Z+1.) - a2_1**(-(Z+1.)/2.) *\
        np.array([gmpy2.cos((Z+1) * gmpy2.atan(2./aa)) for aa in a])

        G12 = (Z+2.)*(Z+1.) * (a**-(Z+3.) + a2_1**(-(Z+3.)/2.) *\
        np.array([gmpy2.cos((Z+3)*gmpy2.atan(2./aa)) for aa in a]))

        G13 = 2.*(Z+1.) * a2_1**(-(Z+2.)/2) *\
        np.array([gmpy2.sin((Z+2.)*gmpy2.atan(2./aa)) for aa in a])

        G1 = G11+G12-G13
        returnVal = 10**I0*p1*G1+bckg
    returnVal = np.array(returnVal.astype(np.float64))
    #print 'Single_schultz calculated with:\nR_av:{} Z:{} I0:{}'.format(R_av, Z, I0)
    #print 'length is:{}, of which nan: {}'.format(len(returnVal), np.sum(np.isnan(returnVal)))
    return returnVal

def single_schultz_spheres_old(q,R_av = 1,Z = 1, I0 = 1, bckg = 0):
    """sing_schultz_spheres: calculates the scattering pattern of an assembly of
    spheres which have a Flory schultz size distribution. the Z parameter is
    defined as z = 1 / sigma^2. THe definistion was taken forom:
    ttp://sasfit.ingobressler.net/manual/Schultz-Zimm
        Args
            q (numpy.array): the array containg the list of q-values for which to
                calculate the scattering
            R_av (int): the mean of the size distribution. Defaults to 1
            Z (int): the dispersion of the distribution. For a Flory-Schultz
                distribution the Z parameter is defined as Z = 1/sigma^2.
                Defaults to 1
            I0 (int): the prefactor which includes information on the scattering
                length density (SLD) and the concentration of particles. Defaults
                to 1
            bckg (int): the background value to use in case the background is
                not perfectly subtracted. Defaults to 0.
        Returns
            the scattering curve which has the same size as q
    """

    a = (Z+1.)/(q*R_av)

    P_q = 8.*np.pi**2*R_av**6*(Z-1.)**(-6.)*a**(Z+7.)*(a**(-(Z+1.))- \
        (4.+a**2)**(-(Z+1.)/2)*np.cos((Z+1.)*np.arctan(2/a)) + \
        (Z+2.)*(Z+1.)*(a**(-Z-3.)+(4+a**2)**((-Z-3.)/2.)*np.cos((Z+3.)*np.arctan(2./a))) - \
        2.*(Z+1.)*(4.+a**2.)**(-(Z+2.)/2.)*np.sin((Z+2.)*np.arctan(2./a)))
    return np.nan_to_num(10**I0*P_q+bckg)


def double_schultz_spheres(q, R1_av = 1, Z1 = 1, R2_av = 1,Z2 = 1, I0 = 1, ratio = 0.5, bckg = 0):
    """double_schultz_spheres: calculates the scattering pattern of an assembly of
    spheres which have a bimodal Flory Schultz distribution.
        Args
            q (numpy.array): the array containg the list of q-values for which to
                calculate the scattering
            R_av1 (int): the mean of the size distribution of the first
                peak. Defaults to 1
            Z1 (int): the dispersion of the first distribution. For a Flory-Schultz
                distribution the Z parameter is defined as Z = 1/sigma^2.
                Defaults to 1
            R_av2 (int): the mean of the size distribution of the second
                peak. Defaults to 1
            Z2 (int): the dispersion of the second distribution. For a Flory-Schultz
                distribution the Z parameter is defined as Z = 1/sigma^2.
                Defaults to 1
            I0 (int): the pre-factor which includes information on the scattering
                length density (SLD) and the concentration of particles. Defaults
                to 1
            ratio (int): the ratio between the first and the second peak. Defaults
                to 0.5
            bckg (int): the background value to use in case the background is
                not perfectly subtracted. Defaults to 0.
        Returns
            the scattering curve which has the same size as q
    """
    return np.nan_to_num(ratio*single_schultz_spheres(q,R1_av,Z1,I0,0)+(1-ratio)*single_schultz_spheres(q,R2_av,Z2,I0,0)+bckg)


def monodisperse_cube(q, L=1, I0=1, bckg = 0):
    """
    http://www.sasview.org/sasview/user/models/model_functions.html#rectangularprismmodel
    :param q: the wavevector, vna be aither a number of a numpy array
    :param L: The side of the cube
    :param I0: The prefactor in front of the form factor
    :param bckg: The constant background to sum
    :return: The complete, integrated form factor for a cube
    """
    def FF(theta, phi):
        A = q*L/2.*np.cos(theta)
        B = q*L/2.*np.sin(theta)*np.sin(phi)
        C = q*L/2.*np.sin(theta)*np.cos(phi)
        return np.sinc(A)*np.sinc(B)+np.sinc(C)
    return 10**I0*dblquad(FF, 0, np.pi/2., lambda x: 0, lambda x: np.pi/2.0)[0]+bckg


def single_gaussian_cube(q, L_av=1, sigma=1, I0=1, bckg = 0):
    def FF(theta,phi,L):
        A = q*L/2.*np.cos(theta)
        B = q*L/2.*np.sin(theta)*np.sin(phi)
        C = q*L/2.*np.sin(theta)*np.cos(phi)
        return single_gauss_distribution(L,L_av,sigma,1)*np.sinc(A)*np.sinc(B)+np.sinc(C)

    l_min = max(0,L_av-4*(L_av*sigma))
    l_max = L_av+4*(L_av*sigma)
    return 10**I0*tplquad(FF, 0, np.pi/2., lambda x: 0, lambda x: np.pi/2.0,)[0]+bckg




