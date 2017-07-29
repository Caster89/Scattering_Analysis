from config import _AvlbUnits, _UnitsSymbols,_UnitsConv, _AvlbSASFit, _AvlbSASFitDic, _AvlbSASFitDicInv, _lmfitModels, _lmfitModelFunctions, _lmfitDistFunctions
import numpy as np
import scipy as sp
from scipy import signal
from scipy import interpolate
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

'''
The distriutions presented here can be found in Polydispersity analysis of scattering
data from self-assembled systems. Phy. Rev. A,45, 2428-2438.
DOI: 10.1103/PhysRevA.45.2428
'''
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
    #return I0*(GaussSpheres(q,R_av1, sigma1,ratio,0)+GaussSpheres(q,R_av2, sigma2,1-ratio,0))
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
    print 'Single_schultz calculated with:\nR_av:{} Z:{} I0:{}'.format(R_av, Z, I0)
    print 'length is:{}, of which nan: {}'.format(len(returnVal), np.sum(np.isnan(returnVal)))
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
def guess_from_dist(xs,ys,fitType = None, goodness = False, max_mode = 'wrap', verbose = True):
    """guess_from_dist: returns the parameter suggestions for a distribution curve.
        Args:
            xs (list): a list of the x coordinates of the curve
            ys (list): a list with the y coordinates of the curve
            fitType (str): The distribution with which to fit the data provided.
                if None then all the available distributions are tried and the
                one returning the lowest redchi is chosen. The redchi of double
                distributions is "weighed" by 0.9, in order to avoid adding
                unecessary complications to the model. Defaults to None
            goodness (bool): whether the function should return the redchi.
                Defaults to False
            max_mode (string): How the edges of the vector are treated. 'wrap'
            (wrap around) or 'clip' (treat overflow as the same as the last (or
            first) element). Default is 'wrap'
        Returns:
            paramSugg: a dictionary of the values of the parameters found
            fitType: returns the fitType used, nevessary if None was passed in
                order to determine the best fitting function
            redchi: the reduced Chi of the fit
    """

    if verbose:
        logging.basicConfig(format='%(levelname)s:%(message)s', level=10)
    else:
        logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

    maxIdx = sp.signal.argrelextrema(ys, np.greater, order = 3,\
        mode = max_mode )[0]

    numbPeaks = len(maxIdx)

    maxValues = ys[maxIdx]
    maxPos = xs[maxIdx]
    #The peaks are ordered not only by intensity, but by weighing the relative
    #volume of the particles.
    maxIdxOrd = (-maxValues*maxPos**3).argsort()[:numbPeaks]
    maxIdx = maxIdx[maxIdxOrd]
    if verbose:
        print 'Fitting_Models:guess_from_dist(): The indixes of the maxima found are: {}'.format(maxIdx)
        print 'Fitting_Models:guess_from_dist(): The maxima intensities are: {}'.format(ys[maxIdx])

    if fitType is None:
        chosenFit = 'Sing_Gauss'
        redchi = np.inf
        paramSugg = {}
        for currFit in _lmfitDistFunctions.keys():
            tempParamSugg, _, tempchi = guess_from_dist(xs,ys,fitType = currFit, goodness = True, verbose = verbose)
            if chosenFit.split('_')[0] == currFit.split('_')[0]:
                if tempchi<redchi:
                    redchi = tempchi
                    chosenFit = currFit
                    paramSugg = tempParamSugg
            elif 'Sing' in chosenFit:
                if tempchi<0.9*redchi:
                    redchi = tempchi
                    chosenFit = currFit
                    paramSugg = tempParamSugg
            else:
                if 0.9*tempchi<redchi:
                    redchi = tempchi
                    chosenFit = currFit
                    paramSugg = tempParamSugg
        if goodness:
            return [paramSugg, chosenFit, redchi]
        else:
            return [paramSugg,chosenFit]


    if numbPeaks == 0:
        if 'Sing' in fitType:
            maxIdx = np.array([5])
        elif 'Double' in fitType:
            maxIdx = np.array([5,5])
    elif numbPeaks == 1:
        maxIdx = np.array([maxIdx[0]])
        if 'Double' in fitType:
            maxIdx = np.array([maxIdx, maxIdx])
    elif numbPeaks > 1:
        if 'Sing' in fitType:
            maxIdx = np.array([maxIdx[0]])
        else:
            maxIdx = np.array(maxIdx[:2])



    maxValues = np.array([ys[maxIdx]]).flatten()

    if verbose:
        print 'Fitting_Models:guess_from_dist():After manipulation the maxima indexes are:{}'.format(maxIdx)
        print 'Fitting_Models:guess_from_dist(): The peaks intensities array is of type: {} with values: {}'.format(type(maxValues),maxValues)
        print 'Trying to use function: {}'.format(_lmfitDistFunctions[fitType])
    maxPos = np.array([xs[maxIdx]]).flatten()
    #print'Fitting_Models:guess_from_dist(): The maxima array is of type: {}\nwith radius values: {}'.format(type(maxPos),maxPos)
    #logging.info('The maxima considered are: {}'.format(maxPos))

    model = lmfit.Model(globals()[_lmfitDistFunctions[fitType]])
    params = model.make_params()
    params['I0'].set(value = np.log10(maxValues[0]), min = -20)
    if len(maxValues)==1:
        params['R_av'].set(value = maxPos[0], min = 0.3)
        if 'Schultz' in fitType:
            params['Z'].set(value = 40, min = 0)
        else:
            params['sigma'].set(value = maxPos[0]*0.1, min = 0)
    else:
        params['ratio'].set(value= 0.5, max = 1.0, min = 0.0)
        params['R1_av'].set(value = maxPos[0], min = 0.3)
        params['R2_av'].set(value = maxPos[1], min = 0.3)
        if 'Schultz' in fitType:
            params['Z1'].set(value = 40, min = 0)
            params['Z2'].set(value = 40, min = 0)
        else:
            params['sigma1'].set(value = maxPos[0]*0.1, min = 0)
            params['sigma2'].set(value = maxPos[1]*0.1, min = 0)
    #print 'Xs has {} elements (of which {}nan) while ys has {} (of which {}nan)'.format(xs.shape, np.sum(np.isnan(xs)), ys.shape, np.sum(np.isnan(ys)))
    #print params
    #logging.debug('Xs:{} \n ys {}'.format(xs,ys))
    if verbose:
        print 'Starting parameters:\n{}'.format(params.valuesdict())
    result = model.fit(ys, params = params, x = xs, fit_kws = {'nan_policy' : 'omit'})
    print 'guess from dist with dist {} found:\n{}'.format(fitType,result.params.valuesdict())
    if goodness:
        return [result.params.valuesdict(), fitType, result.redchi]
    else:
        return [result.params.valuesdict(),fitType]
