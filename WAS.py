from config import _AvlbUnits, _UnitsConv, _AvlbSASFit, _lmfitModels, _lmfitModelFunctions, _lmfitDistFunctions 
from .Scattering_Object import ScatteringObject
import Fitting_Models
import Expected_Maximization as EM
import numpy as np
import scipy as sp
import scipy.signal
import numpy.matlib
import logging
#import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Qt5Agg")
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D
try :
	import lmfit
except ImportError:
	print "Could not import the lmfit package. Some fittings will be disabled."
	lmfitAvlb = False
else:
	lmfitAvlb = True
	import Fitting_Models


	 
class WideAngleScattering(ScatteringObject):
	
	def __init__(self, fname, **kwargs):
		self.PVMFit = None
		super(WideAngleScattering, self).__init__(fname, **kwargs)
		
	def plot_data(self, qRange=[0,np.inf],yShift=0, ax = None, figSize = (6,6), xUnits = 'nm', yUnits ='m',\
					**kwargs):
		"""plot_data: Plots the data. If the axis is passed then the data is plotted
			on the given axis, if not a figure of size figSize is created and the axis places in
			it. The units on the x and y axis can be selected. The keyword arguments are used
			to set the parameters of the plots.
				Args:
					qRange (list): the q-range over which the data should be plotted. Defaults to
						[0, np.inf]
					yShift (int): indicated by how much the data should be shifted vertically. This
						allows the plot_data function to be used in time resolved data. The data is
						added to shift. Defaults to 0
					ax (:obj: plt.axes): the axis on which to plot the data. If it is not provided
						then one will be created. Defaults to None
					figSize (tuple): used to set the size of the figure in case the axis is not
						provided. Defaults to (6,6)
					xUnits (string): the units used on the x-axis. They have to be defined in the
						_AvlbUnits in the config.py file. Defaults to nm
					yUnits (string): the units used on the xy-axis. They have to be defined in the
						_AvlbUnits in the config.py file. Defaults to m
					**kwargs: all the arguments used to set the values of the plot.
				Returns:
					ax (:obj: plt.axes): the axis on which the data was plotted
					fig (:obj: lt.fig): the figure on which the data was plotted
			"""
		for k in self.plot_dict:
			if k in kwargs:
				self.plot_dict[k] = kwargs[k]
		if ax is None:
			fig = Figure(figsize = figSize)
			ax = fig.add_subplot(111)
		tempq = self.q[np.logical_and(self.q>min(qRange), self.q<max(qRange))]
		tempI = self.Iabs[np.logical_and(self.q>min(qRange),self.q<max(qRange))]
		
		qConvCoeff = 1
		IConvCoeff = 1
		if (xUnits in _AvlbUnits) and (self.qUnits in _AvlbUnits):
			qConvCoeff = _UnitsConv[self.qUnits]/_UnitsConv[xUnits]	

		if (yUnits in _AvlbUnits) and (self.IUnits in _AvlbUnits):
			IConvCoeff = _UnitsConv[self.IUnits]/_UnitsConv[yUnits]	

		 
		ax.plot(tempq*qConvCoeff,tempI*IConvCoeff+(yShift),color = self.plot_dict['color'],\
				linestyle = self.plot_dict['linestyle'],linewidth = self.plot_dict['linewidth'],\
				marker = self.plot_dict['marker'], ms = self.plot_dict['ms'],\
				markeredgecolor = self.plot_dict['mec'])
		#Set the x an y limits
		ax.set_ylim(min(tempI*IConvCoeff+(yShift)),max(tempI*IConvCoeff*10**(yShift)))
		ax.set_xlim(min(tempq*qConvCoeff),max(tempq*qConvCoeff))
		 
		#Set the labels for the 2 axis and the tick parameters
		ax.set_xlabel(r'Wavevector ({}$^{{-1}}$)'.format(xUnits), fontsize = kwargs.get('lableSize',20))
		ax.set_ylabel(r'Intensity ({}$^{{-1}}$)'.format(yUnits), fontsize = kwargs.get('lableSize',20))
		ax.tick_params(labelsize = kwargs.get('labelSize',18),size = kwargs.get('tickSize',6),	
		 				width = kwargs.get('tickWidth',3) )
		#implement a method to place the major and minor tick marks
		
		#plt.show()
		#plt.pause(0.001)
		return [ax]
		 
	def fit_peaks(self, fitType = 'PVM', **kwargs):
		"""fit_peaks: fits the wide-angle peaks from the scattering.
			Args:
				fitType (String): the type of distribution used to fit
					the peaks. Defaults to PVM
				**kwargs(dict): dictionary of possible parameters to pass
					to the function
						numbPeaks (int): the number of peaks to be fitted
						peakWidth (int): the average width of the peak in
							number of points. Defaults to 1/10th of the
							number of point in vector q.
		"""
		
		numbPeaks = kwargs.get('numbPeaks',0)
		qRange = kwargs.get('qRange', np.array([0,np.inf]))
		tempq = self.q[np.logical_and(self.q>min(qRange), self.q<max(qRange))]
		tempI = self.Iabs[np.logical_and(self.q>min(qRange),self.q<max(qRange))]
		
		
		#Use scipy to find the indices of the local maxima. Order tells how
		#many point on either side to consider
		order = max(int(len(tempq)/kwargs.get('peakWidth',20)),10)
		maxima = sp.signal.argrelextrema(tempI, np.greater, order = order)[0]
		if numbPeaks == 0:
			numbPeaks = len(maxima)
		else:
			#If a set number of peaks was given all the local maxima are take
			#and the numbPeaks having the largest value are chosen
			if numbPeaks >len(maxima):
				print "The number of requested peaks ({}) is superior to the number".format(numbPeaks)
				print "of relevant peaks found ({}). The fit might result \"overfit\".".format(len(maxima))
			#Use scipy local maxima to find the position of the local maxima. The points are
			#ordered by intensity and the first numbPeaks are selected. The extremes (in
			#the left extreme) are not considered maximums. This gives an approximate position
			# for the position of the peak, which can be passed to the create_MPV_model.
			else:
				#The intensities of all the peaks found and stored in maxValues
				maxValues = tempI[maxima]
				#The first 'numbPeaks' indexes of the peaks . argsort returns an
				#increasing array. -maxValues is used to obtain an increasing order
				maxIdx = (-maxValues).argsort()[:numbPeaks]
				#The position of the first numbPeaks, ordered by intensity
				maxima = maxima[maxIdx]
				#order the peaks by their position
				maxima.sort()
				
		peaksInfo = []
		
		for m in xrange(numbPeaks):
			if m < len(maxima):
				peakDict = {'Center': {'value': tempq[maxima[m]],\
									'min': 0.9*tempq[maxima[m]],\
									'max': 1.1*tempq[maxima[m]],\
									'vary': True}}

			else:
				peakDict = {'Center': {'value': 1,\
									'min': 0,\
									'max': None,\
									'vary': True}}
			peaksInfo.append(peakDict.copy())
		if fitType == 'PVM':
			model, params = self.create_MPV_model(numbPeaks, bckg = kwargs.get('bckg','Cubic'),\
												  PeaksInfo = peaksInfo, q = tempq, I = tempI)

			result = model.fit(tempI, params = params, x=tempq)

			if kwargs.get('plot', False):

				fig = Figure(figsize = (12,6))
				ax = fig.add_subplot(111)
				#ax2 = fig.add_subplot(122)
				ax.plot(tempq, tempI, 'k', linestyle = 'None', marker = 'o')

				ax.plot(tempq, result.best_fit,'r')
				#plt.show()
			self.PVMFit = {'qRange': qRange, 'ChiSqr': result.chisqr, 'RedChi': result.redchi,\
					   'Residuals': result.residual, 'Params': {}, 'bckgParams':{}, 'NumbPeaks': numbPeaks}
		
			for p in result.params:
				#Save the background paramters separatly from the distribution parameters
				if 'bckg' in p:
					self.PVMFit['bckgParams'][p] = result.params[p].value
				else:
					#The full width half maximum is not important
					if 'fwhm' not in p:
						self.PVMFit['Params'][p] = result.params[p].value
				
	def plot_fit(self, peakType = 'PVM', ax = None, qRange = None, **kwargs):
		"""plot_fit: plots the fit of the peaks if available. If the axis is provided
		it plots in the given axis if not a new one is created anc used. If the qRange
		is provided it will plot the fit over the gven q range, even though the curve outside
		of the fitted range will be plotted using a dashed line. If None is provided then it
		will only plot withing the fitted range
		"""
		if peakType == 'PVM':
			if self.PVMFit is None:
				print 'The data has not been fitted yet using a Pseudo-Voigt model.'
				return 0
		else:
			print 'The only available fit model for the moment is the Pseudo-Voigt model (PVMM).'
			print '{} is not a valid model'.format(peakType)
			return 0
		
		if ax is None:
			fig = Figure(figsize = kwargs.get('figSize',(6,6)))
			ax = fig.add_subplot(111)
		if qRange is None:
			qRange = self.PVMFit['qRange']
		
		tempq = self.q[np.logical_and(self.q > min(qRange), self.q < max(qRange))]
		tempI = self.Iabs[np.logical_and(self.q > min(qRange), self.q < max(qRange))]
		ax.plot(tempq, tempI, color = self.fit_plot_dict['color'], linestyle = self.fit_plot_dict['linestyle'],linewidth = self.fit_plot_dict['linewidth'],\
						marker = self.fit_plot_dict['marker'], ms = self.fit_plot_dict['ms'],\
						markeredgecolor = self.fit_plot_dict['mec'])
		
		bckgDegree = len(self.PVMFit['bckgParams'])-1
		bckgModel = lmfit.models.PolynomialModel(degree = bckgDegree, prefix = 'bckg_')
		bckgParams = bckgModel.make_params()
		for p in bckgParams:
			bckgParams[p].set(value = self.PVMFit['bckgParams'][p])
		if peakType == 'PVM':
			numbPeaks = int(len(self.PVMFit['Params'])/4)
			model = lmfit.models.PseudoVoigtModel(prefix = 'pvm0_')
			for x in range(1,numbPeaks):
				model += lmfit.models.PseudoVoigtModel(prefix = 'pvm{}_'.format(x))
			
			#model += bckgModel
			params = model.make_params()
			for p in params:
				if 'fwhm' not in p:
					params[p].set(value = self.PVMFit['Params'][p])
		
		model += bckgModel
		params += bckgParams
		if qRange == self.PVMFit['qRange']:
			ax.plot(tempq, model.eval(params = params, x = tempq), color = self.fit_plot_dict['fit_color'],\
					linestyle = self.fit_plot_dict['fit_linestyle'],linewidth = self.fit_plot_dict['fit_linewidth'],\
						marker = self.fit_plot_dict['fit_marker'], ms = self.fit_plot_dict['fit_ms'],\
						markeredgecolor = self.fit_plot_dict['fit_mec'])
		else:
			qFit = self.q[np.logical_and(self.q > min(self.PVMFit['qRange']), self.q < max(self.PVMFit['qRange']))]
			ax.plot(tempq, model.eval(params = params, x = tempq),color = self.fit_plot_dict['fit_color'],\
					linestyle = self.fit_plot_dict['fit_linestyle'],linewidth = self.fit_plot_dict['fit_linewidth'],\
						marker = self.fit_plot_dict['fit_marker'], ms = self.fit_plot_dict['fit_ms'],\
						markeredgecolor = self.fit_plot_dict['fit_mec'])
			ax.plot(qFit, model.eval(params = params, x = qFit),color = self.fit_plot_dict['fit_color'],\
					linestyle = self.fit_plot_dict['fit_linestyle'],linewidth = self.fit_plot_dict['fit_linewidth'],\
						marker = self.fit_plot_dict['fit_marker'], ms = self.fit_plot_dict['fit_ms'],\
						markeredgecolor = self.fit_plot_dict['fit_mec'])
		return ax
	
	def get_fit_params(self, fitType, param = None):
		'''get_fit_params returns the parameter or parameters requested. If a parameter
		is specified in the param argument it will return that argument. If the name of
		the argument is preceded by the peak number (e.g. 'pvm0_center') only one parameter
		will be returned. If the general name is provided (e.g. center) then the parameters
		for all the peaks are returned. If no parameter is specified then all the parameters
		are returned. Passing 'bckg' as parameter will return the coefficients of the backg.
			Args:
			
			fitType (String): the peak distribution used for the fit.
			param (String): the parameter required. If None returns
				all the parameters used for the fit. Defaults to None
		'''
		if fitType == 'PVM':
			if self.PVMFit is None:
				return None
			if param  is None:
				return self.PVMFit
			else:
				if param.startswith('pvm'):
					return {param: self.PVMFit['Params'].get(param,None)}
				elif param == 'bckg':
					return self.PVMFit['bckParams']
				else:
					result = {}
					for i in range(self.PVMFit['NumbPeaks']):
						result['pvm{}_{}'.format(i, param)] = self.PVMFit['Params'].get('pvm{}_{}'.format(i, param), None)
					return result
	
	def get_peak_position(self, fitType, origPos = 1):
		centers = []
		origPos = np.array(origPos)
		if self.PVMFit is not None:
			for i in xrange(self.PVMFit):
				centers.append(self.PVMFit['Params']['pvm{}_center'.format(i)])
			centers = np.array(centers)
			centers = np.sort(centers)
			
			if origPos != 1:
				if isinstance(origPos, int):
					centers /= origPos
				elif isinstance(origPos, np.array):
					centers = centers/np.sort(origPos)
				else:
					print 'origPos should be either int or a list'
		return centers

	def get_peak_intensity(self, fitType):
		intensity = []
		if self.PVMFit is not None:
			for i in xrange(self.PVMFit):
				intensity.append(self.PVMFit['Params']['pvm{}_center'.format(i)])
			intensity = np.array(intensity)
			
		return intensity
	
	def get_peak_intensity(self, fitType, K = 0.9):
		tau = []
		if self.PVMFit is not None:
			for i in xrange(self.PVMFit):
				mu = self.PVMFit['Params']['pvm{}_center'.format(i)]
				fwhm = self.PVMFit['Params']['pvm{}_sigma'.format(i)]**2
				tempTau = K*2*np.pi
				tau.append(tempTau)
			tau = np.array(tau)
			
		return tau
	
	def create_MPV_model(self, NumbPeaks, **kwargs):
		"""create_MPV_model: used to create a model with multiple pseudo-Voigt peaks
		and a background.
			Args:
				BckgParams (list of dict): a list of dict with the necessary information for the
					background curve
						value (int): the value of the coefficient of te curve
						min (int): the minimum  value the coefficient can have
						max (int): the maximum value the coefficient can have
						vary (bool): whether the coefficient should vary
			Returns:
				The complete model and parameters needed to fit the WAXS curve
				
		"""
		bckgPower = {'Linear': 1, '2-Quadratic': 2, 'Cubic' : 3, '4-Quadratic': 4}
		
		#bckgMP = 
		bckg = lmfit.models.PolynomialModel( 3, prefix = 'bckg_')
		bckgParams = bckg.make_params()
		if kwargs.get('BckgParams', None) is not None:
			bgparams = kwargs.get('BckgParams')
			for bgp in bgparams:
				bckgParams['bckg_'+bgp].set(value = bgp.get('value', 1), min = bgp.get('min', -np.inf),\
								max = bgp.get('max', np.inf), vary = bgp.get('vary', True))
		else:
			for p in bckgParams:
				bckgParams[p].set(value = 0)
		
		model = lmfit.models.PseudoVoigtModel(prefix = 'pvm0_')
		params = model.make_params()
		peaksInfo = kwargs.get('PeaksInfo', None)
		if peaksInfo is not None:
				#Set-up the parameters using the information provided
				params = self._create_PV_params(params,prefix = 'pvm0_', paramDict = peaksInfo[0])
		for np in xrange(NumbPeaks-1):
			tempModel = lmfit.models.PseudoVoigtModel(prefix = 'pvm{}_'.format(np+1))
			tempParams = tempModel.make_params()
			if kwargs.get('PeaksInfo', None) is not None:
				#Set-up the parameters using the information provided
				tempParams = self._create_PV_params(tempParams,prefix = 'pvm{}_'.format(np+1),\
													paramDict = peaksInfo[np+1])

			model += tempModel
			params += tempParams
		model += bckg
		params += bckgParams
		return [model,params]
		
	def _create_PV_params(self, params, prefix, paramDict):
		"""Helper function used to create the set the values for the
		parameter. Implemented to the create_MPV_model is not too
		cluttered.
			Args:
				params (lmfit.Parameters): lmfit parameters object used
					to fit the data. All the parameters should exist
				prefix (string): the characteristic string of the parameter
					used when multiple models are concatenated
				paramDict (dict of dict): dictionary of dictionaries containing
					all the values for the parameters. Any value not found
					will be set at a default value.
			Returns:
				the parameter variable with all the values initiated
		"""
		if paramDict.get('Center',None) is not None:
			params['{}center'.format(prefix)].set(value = paramDict['Center'].get('value',1),\
												   min = paramDict['Center'].get('min',0),\
												   max = paramDict['Center'].get('max',np.inf),\
												   vary = paramDict['Center'].get('vary',True))
		else:
			params['{}center'.format(prefix)].set(value = 1, min = 0, max = np.inf, vary = True)
		
		if paramDict.get('Sigma',None) is not None:
			params['{}sigma'.format(prefix)].set(value = paramDict['Sigma'].get('value',1),\
												   min = paramDict['Sigma'].get('min',0),\
												   max = paramDict['Sigma'].get('max',np.inf),\
												   vary = paramDict['Sigma'].get('vary',True))
		else:
			params['{}sigma'.format(prefix)].set(value = 1, min = 0, max = np.inf, vary = True)
			
		if paramDict.get('Amplitude',None) is not None:
			params['{}amplitude'.format(prefix)].set(value = paramDict['Amplitude'].get('value',1),\
												   min = paramDict['Amplitude'].get('min',0),\
												   max = paramDict['Amplitude'].get('max',np.inf),\
												   vary = paramDict['Amplitude'].get('vary',True))
		else:
			params['{}amplitude'.format(prefix)].set(value = 1, min = 0, max = np.inf, vary = True)
			
		if paramDict.get('Fraction',None) is not None:
			params['{}fraction'.format(prefix)].set(value = paramDict['Fraction'].get('value',0.5),\
												   min = paramDict['Fraction'].get('min',0),\
												   max = paramDict['Fraction'].get('max',1),\
												   vary = paramDict['Fraction'].get('vary',True))
		else:
			params['{}fraction'.format(prefix)].set(value = 0.5, min = 0, max = 1, vary = True)
			
		return params