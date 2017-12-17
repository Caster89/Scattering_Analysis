from lmfit import Parameter, Parameters
import numpy as np
import logging
from config import _lmfitParamsOrder

logger = logging.getLogger(__name__)

def create_parameters(parameter_list):
    """Create the parameters variable for the lmfit models to use. The information needed is passed as a dictionary
    (possible OrderedDict) of dictionaries with all the parameters in the right order

    :param parameter_list: (dict) The keys are the names of the parameters while each value is a dictionary itself
        which contains the value, min, max and vary parameters to initialize the parameter. All these values have
        default values.
    :return: params: (Parameters) The parameters which can be used to fit an lmfit model.
    """
    logger.debug('Creating the parameters for the fit starting from the dictionary:\n{}'.format(parameter_list))
    params = Parameters()
    try:
        for p in _lmfitParamsOrder:
            if parameter_list.get(p, None):
                logger.debug('Creating parameter:{}\n'.format(p),parameter_list[p])
                param = parameter_list.get(p)
                params.add(p,
                           value=param.get('value',1),
                           min=param.get('min', -np.inf),
                           max=param.get('max', np.inf),
                           vary=param.get('vary', True))
        logger.debug('Created the parameters:\n{}'.format(params))
        return params
    except Exception as exc:
        logger.error('There was an error with the initialization of the Parameters:\n{exc}'.format(exc=exc))
        raise


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
    # The peaks are ordered not only by intensity, but by weighing the relative volume of the particles.
    maxIdxOrd = (-maxValues*maxPos**3).argsort()[:numbPeaks]
    maxIdx = maxIdx[maxIdxOrd]
    if verbose:
        print 'Fitting_Models:guess_from_dist(): The indixes of the maxima found are: {}'.format(maxIdx)
        print 'Fitting_Models:guess_from_dist(): Corresponding to radii of: {}'.format(xs[maxIdx])
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
    #print 'guess from dist with dist {} found:\n{}'.format(fitType,result.params.valuesdict())
    if goodness:
        return [result.params.valuesdict(), fitType, result.redchi]
    else:
        return [result.params.valuesdict(),fitType]

