from config import *
import numpy as np
import scipy as sp
from lmfit import minimize, Parameter, report_fit

def single_gauss_spheres(q,R_av = 1,sigma = 1,I0 = 1,bckg=0):
    """sing_gauss_spheres: calculates the scattering pattern of an assembly of
    spheres which have a Gaussian size distribution.
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
    return 10**(I0)*P_q+bckg

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
    return ratio*single_gauss_spheres(q,R1_av, sigma1,I0,0)+(1-ratio)*single_gauss_spheres(q,R2_av, sigma2,I0,0)+bckg
   #return I0*(GaussSpheres(q,R_av1, sigma1,ratio,0)+GaussSpheres(q,R_av2, sigma2,1-ratio,0))

def single_schultz_spheres(q,R_av = 1,Z = 1, I0 = 1, bckg = 0):
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
    return I0*P_q+bckg

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
    return ratio*single_schultz_spheres(q,R1_av,Z1,I0,0)+(1-ratio)*single_schultz_spheres(q,R2_av,Z2,I0,0)+bckg

def single_gauss_distribution(x, R_av, sigma, I0):
    """single_gauss_distribution: returnd the PDF for a normal distribution
    given the necessary parameters.
        Args:
            x (:obj: numpy.array): the list of values for which the PDF has to be
                calculated
            R_av (int): the mean of the distribution
            sigma (int): the dispersion of the distribution
            I0 (int): the intensity of hte distribution
        Returns
            The vector with hte propability densities with the same size as x.
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
    shape=Z+1
   #scale=1
   #gamma=bins**(shape-1)*(np.exp(-bins/scale) /(sp.special.gamma(shape)))

    return 10**(I0)*(shape**shape*(x/R_av)**Z*np.exp(-shape*x/R_av)/(R_av*sp.special.gamma(shape)))

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
    return ratio*single_schultz_distribution(x,R1_av,Z1,ratio) + (1-ratio)*single_schultz_distribution(x,R2_av,Z2,I0)