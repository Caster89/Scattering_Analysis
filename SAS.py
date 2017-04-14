from config import _AvlbUnits, _UnitsConv, _AvlbSASFit, _lmfitModels, _lmfitModelFunctions, _lmfitDistFunctions 
from Scattering_Object import ScatteringObject
import Fitting_Models
import Expected_Maximization as EM
import numpy as np
import numpy.matlib
import logging
#import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Qt5Agg")
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import scipy as sp
from scipy import signal
from scipy import interpolate
from mpl_toolkits.mplot3d import Axes3D
try :
	import lmfit
except ImportError:
	print "Could not import the lmfit package. Some fittings will be disabled."
	lmfitAvlb = False
else:
	lmfitAvlb = True
	import Fitting_Models


	 
class SmallAngleScattering(ScatteringObject):
	"""SmallAngleScattering: this class is a child of ScatteringObject used to store
	data arriving from small angle scattering measurements. It provides a method
	to plot the data in log log scale as well as to fit it. There are several
	different fitting methods available. The different methods can be divided in
	2 main categories: those using the lmfit library and the one based on the
	expected maximization method. The former are not available if the lmfit module
	was not correctly loaded.
	TODO:
		generalize the fitting/potting, and organize it better in order to have a skimmer
		fitting/plotting function.
		Implement the Porod and Guinier.
	"""
	def __init__(self,fname,**kwargs):

		self.singleGauss = None
		self.doubleGauss = None
		self.singleSchultz = None
		self.doubleSchultz = None
		self.EMFit = None
		self.fitResults = {f : None for f in _AvlbSASFit}
		self.porod_invariant = None
		self.guinier_regime = None
		
		super(SmallAngleScattering, self).__init__(fname, **kwargs)
	
	def plot_data(self, qRange=[0,np.inf],yShift=1, ax = None, figSize = (6,6), xUnits = 'nm', yUnits ='m',\
					**kwargs):
		"""plot_data: Plots the data in a log-log plot. If the axis is passed then the data is plotted
		on the given axis, if not a figure of size figSize is created and the axis places in
		it. The units on the x and y axis can be selected. The keyword arguments are used
		to set the parameters of the plots.
			Args:
				qRange (list): the q-range over which the data should be plotted. Defaults to
					[0, np.inf]
				yShift (int): indicated by how much the data should be shifted vertically. This
					allows the plot_data function to be used in time resolved data. The data is
					multiplied by shift because it is plotted in loglog. Defaults to 1
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
				fig (:obj: plt.fig): the figure on which the data was plotted
		"""
		for k in self.plot_dict:
			if k in kwargs:
				self.plot_dict[k] = kwargs[k]
		#plt.ion()
		if ax is None:
			fig = Figure(figsize = figSize)
			ax = fig.add_subplot(111)
		tempq = self.q[np.logical_and(self.q>min(qRange), self.q<max(qRange))]
		tempI = self.Iabs[np.logical_and(self.q>min(qRange),self.q<max(qRange))]
		if self.Ierr is not None:
			tempE = self.Ierr[np.logical_and(self.q>min(qRange),self.q<max(qRange))]
		
		qConvCoeff = 1
		IConvCoeff = 1
		if (xUnits in _AvlbUnits) and (self.qUnits in _AvlbUnits):
			qConvCoeff = _UnitsConv[self.qUnits]/_UnitsConv[xUnits]	

		if (yUnits in _AvlbUnits) and (self.IUnits in _AvlbUnits):
			IConvCoeff = _UnitsConv[self.IUnits]/_UnitsConv[yUnits]
				

		if self.Ierr is not None and kwargs.get('plot_error', False):
			#This part was modified, instead of asking the values every time, a plot dictionary was created
			#to store all the properties, and is updated in case new parameters are passed to the plotting
			#function.
			ax.errorbar(tempq*qConvCoeff,tempI*IConvCoeff*yShift,yerr = tempE*IConvCoeff*yShift, **self.plot_dict)
						
					#color = kwargs.get('color','k'),\
		 			#linestyle = kwargs.get('linestyle', 'None'),marker = kwargs.get('marker','o'),\
		 			#markeredgecolor = kwargs.get('color','k'), ms = kwargs.get('ms',2))
			ax.set_xscale("log", nonposx = "clip")
			ax.set_yscale("log", nonposx = "clip")
		else:
			ax.loglog(tempq*qConvCoeff,tempI*IConvCoeff*yShift, **self.plot_dict)
					  
					  #color = kwargs.get('color','k'),\
						#linestyle = kwargs.get('linestyle', 'None'),marker = kwargs.get('marker','o'),\
						#markeredgecolor = kwargs.get('color','k'), ms = kwargs.get('ms',2))
		#Set the x an y limits
		ax.set_ylim(min(tempI*IConvCoeff*yShift),max(tempI*IConvCoeff*yShift))
		ax.set_xlim(min(tempq*qConvCoeff),max(tempq*qConvCoeff))
		 
		#Set the labels for the 2 axis and the tick parameters
		ax.set_xlabel(r'Wavevector ({}$^{{-1}}$)'.format(xUnits), fontsize = kwargs.get('lableSize',20))
		ax.set_ylabel(r'Intensity ({}$^{{-1}}$)'.format(yUnits), fontsize = kwargs.get('lableSize',20))
		ax.tick_params(labelsize = kwargs.get('labelSize',18),size = kwargs.get('tickSize',6),	
		 				width = kwargs.get('tickWidth',3) )
		#implement a method to place the major and minor tick marks
		
		#plt.show()
		#plt.pause(0.001)
		return ax
		 
	def plot_slopes(self, slope, frequency = 10,slope_color ='k', **kwargs):
		"""plot_slopes: plots the data along with a series of parallel lines with given slope.
		this can be useful when trying ot detrmine the slope of a particlare scattering curve.
			Args:
				slope (int): the power of the polynomial to plot x^(slope). As the slope is
					(almost) always decreasing, an value given will be taken as negative
				frequency (int): the distance between subsequent lines. The nth line will
					therefore be frequency^(n)*x^slope. Defaults to 10
				slope_color (char): a valid character to define the color with which the
					sloped lines should be drawn. Defaults to 'k' (black)
				**kwargs (dict): any keywords which can be used with the plot_data function
			Retruns:
				the axis and the figure of the plot
		"""
		slope = -abs(slope)
		#plt.ion()
		if kwargs.get('ax', None) is None:
			ax, _ = self.plot_data(**kwargs)
		else:
			ax = kwargs.get('ax')
		
		qRange = kwargs.get('qRange',[0,np.inf])
		tempq = self.q[np.logical_and(self.q>min(qRange), self.q<max(qRange))]
		tempI = self.Iabs[np.logical_and(self.q>min(qRange),self.q<max(qRange))]
		qExt = ax.get_xlim()
		IExt = ax.get_ylim()
		lowLeft = [min(tempq), min(IExt)]
		upRight = [max(tempq),max(IExt)]
		baseLine = tempq**(slope)*IExt[0]/tempq[0]**(slope)
		n=0
		visible = True
		while visible:
			tempLine = baseLine*frequency**(n)
			if tempLine[-1]>=IExt[1]:
				visible = False
			ax.loglog(tempq,tempLine,slope_color,linestyle = '--')
			n+=1
		#plt.draw()
		#plt.ioff()
		return[ax,fig]
	
	def fit_data(self, fitType = "Sing_Gauss",**kwargs):
		"""fit_data: wrapper function used to decide which fitting algorithm to use
		and to set it up.
			Args:
				fitType (string): the method used for fitting the data. Has to be defined in _AvlbSASFit
					contained in config.py. Defaults to Sing_Gauss.
				**kwargs (dict):
					The keywords to pass to the defired fitting method.
					plot (bool): whether to plot or not the fit. Defaults to False
					paramSugg (dict): dictionary passed on to the create_params function. It containes
						any suggestion used to create the parameters
					
		"""
		if kwargs.get('verbose', False):
			logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
		else:
			logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
		#Check to make sure that if the fitting method required used the lmfit
		#package, it is in fact available
		if fitType in _lmfitModels:
			if not lmfitAvlb:
				print "Impossible to fit using {} model because lmfit was not important.".format(type)
				return 0
			
			params = self.create_params(fitType, kwargs.get('paramSugg', {}))
			self.lmfit_fit(fitType, params, **kwargs)
			
		if fitType == "EM":
			self.em_fit( **kwargs)
		
		#if kwargs.get('plot',False):
			#self.plot_fit(type)
			
		if fitType not in _AvlbSASFit:
			print "{} is not an implemented fit. The available options are:".format(fitType)
			for ft in _AvlbSASFit:
				print "\t - {}".format(ft)
			return 0
		
	def create_params(self, fitType, paramSugg = {}):
		"""create_params: creates the Parameter object for the lmfit package based on the
		type of fit which is chosen between the available fits.
			Args:
				fitType (string): the method used for fitting. Has to be defined in _lmfitModels
				paramSugg (dict): dictionary containing all the default values for the parameters.
					these depend on the fit used.
			Returns:
				params (:obj: Parameters): a Parameters object from the lmfit package with all
					the parameters needed to fit using the selected model.
					
		"""
		params = lmfit.Parameters()
		if fitType == "Sing_Gauss":
			params.add("R_av", value = paramSugg.get('R_av',10), min = paramSugg.get('R_min',0),\
					 max = paramSugg.get('R_max',None),vary = paramSugg.get('R_vary',True))
			params.add("sigma", value = paramSugg.get('sigma',0.5), min = paramSugg.get('sigma_min',0),\
					 max = paramSugg.get('sigma_max',None),vary = paramSugg.get('sigma_vary',True))
			params.add("I0", value = paramSugg.get('I0',1), min = paramSugg.get('I0_min',None),\
					 max = paramSugg.get('I0_max',None),vary = paramSugg.get('I0_vary',True))
			params.add("bckg", value = paramSugg.get('bckg',0), min = paramSugg.get('bckg_min',None),\
					 max = paramSugg.get('bckg_max',None),vary = paramSugg.get('bckg_vary',False))
					 
		elif fitType == "Double_Gauss":
			params.add("R1_av", value = paramSugg.get('R1_av',10), min = paramSugg.get('R1_min',0),\
					 max = paramSugg.get('R1_max',None),vary = paramSugg.get('R1_vary',True))
			params.add("sigma1", value = paramSugg.get('sigma1',0.5), min = paramSugg.get('sigma1_min',0),\
					 max = paramSugg.get('sigma1_max',None),vary = paramSugg.get('sigma1_vary',True))
			params.add("R2_av", value = paramSugg.get('R2_av',10), min = paramSugg.get('R2_min',0),\
					 max = paramSugg.get('R2_max',None),vary = paramSugg.get('R2_vary',True))
			params.add("sigma2", value = paramSugg.get('sigma2',0.5), min = paramSugg.get('sigma2_min',0),\
					 max = paramSugg.get('sigma2_max',None),vary = paramSugg.get('sigma2_vary',True))
			params.add("ratio", value = paramSugg.get('ratio',0.5), min = paramSugg.get('ratio_min',0),\
					 max = paramSugg.get('ratio_max',1),vary = paramSugg.get('ratio_vary',True))
			params.add("I0", value = paramSugg.get('I0',1), min = paramSugg.get('I0_min',None),\
					 max = paramSugg.get('I0_max',None),vary = paramSugg.get('I0_vary',True))
			params.add("bckg", value = paramSugg.get('bckg',0), min = paramSugg.get('bckg_min',None),\
					 max = paramSugg.get('bckg_max',None),vary = paramSugg.get('bckg_vary',False))
					 
		elif fitType == "Sing_Schultz":
			params.add("R_av", value = paramSugg.get('R_av',10), min = paramSugg.get('R_min',0),\
					 max = paramSugg.get('R_max',None),vary = paramSugg.get('R_vary',True))
			params.add("Z", value = paramSugg.get('Z',10), min = paramSugg.get('Z_min',0),\
					 max = paramSugg.get('Z_max',None),vary = paramSugg.get('Z_vary',True))
			params.add("I0", value = paramSugg.get('I0',1), min = paramSugg.get('I0_min',None),\
					 max = paramSugg.get('I0_max',None),vary = paramSugg.get('I0_vary',True))
			params.add("bckg", value = paramSugg.get('bckg',0), min = paramSugg.get('bckg_min',None),\
					 max = paramSugg.get('bckg_max',None),vary = paramSugg.get('bckg_vary',False))
					 
		elif fitType == "Double_Schultz":
			params.add("R1_av", value = paramSugg.get('R1_av',10), min = paramSugg.get('R1_min',0),\
					 max = paramSugg.get('R1_max',None),vary = paramSugg.get('R1_vary',True))
			params.add("Z1", value = paramSugg.get('Z1',10), min = paramSugg.get('Z1_min',0),\
					 max = paramSugg.get('Z1_max',None),vary = paramSugg.get('Z1_vary',True))
			params.add("R2_av", value = paramSugg.get('R2_av',10), min = paramSugg.get('R2_min',0),\
					 max = paramSugg.get('R2_max',None),vary = paramSugg.get('R2_vary',True))
			params.add("Z2", value = paramSugg.get('Z2',10), min = paramSugg.get('Z2_min',0),\
					 max = paramSugg.get('Z2_max',None),vary = paramSugg.get('Z2_vary',True))
			params.add("ratio", value = paramSugg.get('ratio',0.5), min = paramSugg.get('ratio_min',0),\
					 max = paramSugg.get('ratio_max',1),vary = paramSugg.get('ratio_vary',True))
			params.add("I0", value = paramSugg.get('I0',1), min = paramSugg.get('I0_min',None),\
					 max = paramSugg.get('I0_max',None),vary = paramSugg.get('I0_vary',True))
			params.add("bckg", value = paramSugg.get('bckg',0), min = paramSugg.get('bckg_min',None),\
					 max = paramSugg.get('bckg_max',None),vary = paramSugg.get('bckg_vary',False))
					 
		return params
	
	def lmfit_fit(self, fitType, params, **kwargs):
		"""lmfit_fit: fits the data to one of the model types. kwargs should contain 2
		qRanges, a first one over which the fit is done over the wide qrange where the 
		values are normalized by 1/I. A second at low-q range for fitting the I0. If no
		ranges are provided the whole q range will be used for the first fit and 10% of 
		range will be used for the second fit
			Args:
				fitType (string): the type of distribution used to fit the data
				params (:obj: Parameters): The Parameter object from lmfit with all the
					parameters provided.
				**kwargs (dict): the set of parameters used for the fit.
					qRange (list): the first range of q-values used for fitting with
						weights = 1/I. Defaults to [0,np.inf] 
					RedqRange (list): The second range of q-values used for fitting the I0.
						Defaults to [0,0.1*max(q)]
					plotFit (bool): decides whether or not to plot the fitted data
					fit_kws (dict): contains the dictionary to pass over to the scipy function
						for the fit.
		"""
		if kwargs.get('verbose', False):
			logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
		else:
			logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
			
		model = lmfit.Model(getattr(Fitting_Models, _lmfitModelFunctions[fitType]),\
							independent_variables = ['q'])
		
		#Initial fit over the wider qRange
		qRange = kwargs.get('qRange',np.array([0,np.inf]) )

		tempq = self.q[np.logical_and(self.q>min(qRange), self.q<max(qRange)) ]
		tempI = self.Iabs[np.logical_and(self.q>min(qRange), self.q<max(qRange)) ]

		
		if "Schultz" in fitType:
			print " Due to the form of the Schultz distribution the data has to be \
					shifted to m^(-1). Therefore the intensity fit might be off. Not\
					sure, will have to check."
			convFact = _UnitsConv['m']/_UnitsConv[self.qUnits]
			tempq = tempq*convFact
		
		result = model.fit(tempI, params, q=tempq, weights = 1., fit_kws = kwargs.get('fit_kws',{}))
		logging.debug('_____________FIRST FIT____________\n The fit succeded: {}\n Lmfit message: {}\n The scipy return code is: {}\n The scipy message is: {}\n _______________________'.format(result.success,\
																																									   result.message,\
																																									   result.ier,\
																																									   result.lmdif_message))
		
		#Second fit done over the reduced range in order to correctly fit the
		#intensity
		RedqRange = kwargs.get('RedqRange',[0,0.1*max(self.q)])
		
		tempq = self.q[np.logical_and(self.q>min(RedqRange), self.q<max(RedqRange)) ]
		tempI = self.Iabs[np.logical_and(self.q>min(RedqRange), self.q<max(RedqRange)) ]

		if "Schultz" in fitType:
			convFact = _UnitsConv['m']/_UnitsConv[self.qUnits]
			tempq = tempq*convFact
		#All parameters except I0 are set to constant
		for p in result.params:
			if p is not "I0":
				result.params[p].set(vary = False)
		result = model.fit(tempI, result.params, q = tempq, fit_kws = kwargs.get('fit_kws',{}))
		
		logging.debug('_____________SECOND FIT____________\n The fit succeded: {}\n Lmfit message: {}\n The scipy return code is: {}\n The scipy message is: {}\n _______________________'.format(result.success,\
																																									   result.message,\
																																									   result.ier,\
																																									   result.lmdif_message))
		
		if kwargs.get('plot', False):
			tempq = self.q[np.logical_and(self.q>min(qRange), self.q<max(qRange)) ]
			tempI = self.Iabs[np.logical_and(self.q>min(qRange), self.q<max(qRange)) ]
			logging.info('Plotting still needs to be implemented for pyqt5')
			#plt.loglog(tempq,tempI)
			#plt.loglog(tempq, model.eval(params = result.params, q=tempq))
			#plt.show()
		
		
		#The results of the fit are stored in the correct dictionary
		self.fitResults[fitType] = {p:result.params[p].value for p in result.params}
		self.fitResults[fitType]['qRange'] = kwargs.get('qRange',[0,np.inf])
		self.fitResults[fitType]['redchi'] = result.redchi
		if fitType == 'Sing_Gauss':
			self.singleGauss = {}
			for p in result.params:
				self.singleGauss[p] = result.params[p].value
			self.singleGauss['qRange'] = kwargs.get('qRange',[0,np.inf])
			self.singleGauss['redchi'] = result.redchi
			#self.singleGauss = {'R_av': result.params['R_av'].value, 'sigma': result.params['sigma'].value,\
		 							#'I0': result.params['I0'].value, 'bckg': result.params['bckg'].value, \
		 							#'range': kwargs.get('qRange',[0,np.inf])}
		elif fitType == 'Double_Gauss':
			self.doubleGauss = {}
			for p in result.params:
				self.doubleGauss[p] = result.params[p].value
			self.doubleGauss['qRange'] = kwargs.get('qRange',[0,np.inf])
			#self.doubleGauss = {'R_av1': result.params['R_av1'].value, 'sigma1': result.params['sigma1'].value,\
									#'R_av2': result.params['R_av2'].value, 'sigma2': result.params['sigma2'].value,\
		 							#'I0': result.params['I0'].value, 'bckg': result.params['bckg'].value,\
		 							#'range': kwargs.get('qRange',[0,np.inf])}	
		if fitType == 'Sing_Schultz':
			self.singleSchultz = {}
			for p in result.params:
				self.singleSchultz[p] = result.params[p].value
			self.singleSchultz['qRange'] = kwargs.get('qRange',[0,np.inf])
			self.singleSchultz['ConvFact'] = convFact
			#self.singleSchultz = {'R_av': result.params['R_av'].value, 'Z': result.params['Z'].value,\
		 							#'I0': result.params['I0'].value, 'bckg': result.params['bckg'].value,\
		 							#'range': kwargs.get('qRange', [0, np.inf]), 'ConvFact': convFact}
		elif fitType == 'Double_Schultz':
			self.doubleSchultz = {}
			for p in result.params:
				self.doubleSchultz[p] = result.params[p].value
			self.doubleSchultz['qRange'] = kwargs.get('qRange',[0,np.inf])
			self.doubleSchultz['ConvFact'] = convFact
			#self.singleSchultz = {'R_av1': result.params['R_av1'].value, 'Z1': result.params['Z1'].value,\
									#'R_av2': result.params['R_av2'].value, 'Z2': result.params['Z2'].value,\
		 							#'I0': result.params['I0'].value, 'bckg': result.params['bckg'].value,\
		 							#'range': kwargs.get('qRange',[0,np.inf]),'ConvFact': convFact}
		
		if kwargs.get('plotFit',False):
			self.plot_fit(fitType)
		
	def em_fit(self, **kwargs):
		if kwargs.get('verbose', False):
			logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
		else:
			logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
			
		convFact = _UnitsConv[kwargs.get('units','nm')]/_UnitsConv[self.qUnits]
		logging.debug('orignial number of points (q-values): {}'.format(len(self.q)))
		qRange = kwargs.get('qRange',np.array([0,np.inf]) )
		logging.debug("the given qRange is: {}, while the extremes of the q vector are: {}, {}".format(qRange,self.q[0],self.q[-1]))

		tempq = self.q[np.logical_and(self.q>min(qRange), self.q<max(qRange)) ]
		tempI = self.Iabs[np.logical_and(self.q>min(qRange), self.q<max(qRange)) ]
		logging.debug('After applying qRange the number of points is : '.format(len(tempq)))
		xk, Rvec, H = EM.expected_maximization_method(tempq*convFact, tempI, Rmin = kwargs.get('Rmin',0.1),\
							Rmax = kwargs.get('Rmax',0), eps = kwargs.get('eps', 1e-2),\
							k=kwargs.get('k',3), numbElements=kwargs.get('numbElements',100),\
							maxIter = kwargs.get('maxIter',10000), verbose = kwargs.get('verbose',False))
		self.EMFit = {'Rvec': Rvec, 'H': H, 'xk': xk, 'units': kwargs.get('units','nm'), 'qRange': qRange}
		self.fitResults['EM'] = {'Rvec': Rvec, 'H': H, 'xk': xk, 'units': kwargs.get('units','nm'), 'qRange': qRange}
	
	def lmfit_from_em(self, qRange, **kwargs):
		'''Uses the results from the EM fit as a starting point for the lmfit method.
		looks for at most 2 peaks and uses gaussian distributions. Implemetatio for the Schultz
		distribution still have to be made
		'''
		if kwargs.get('verbose', False):
			logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
		else:
			logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
		
		qRange = kwargs.get('qRange',np.array([0,np.inf]) )
		
		if self.fitResults['EM'] is None:
			logging.error('Expected maximization has not been performed on the sample. It will be performed using the default arguments')
			self.em_fit(qRange = qRange)
		
		Rvec = self.fitResults['EM']['Rvec']
		xk = self.fitResults['EM']['xk'].flatten()
		
		order = 1
		#max(int(len(Rvec)/kwargs.get('peakWidth',5)),5)
		maxima = sp.signal.argrelextrema(xk, np.greater, order = order)[0]
		logging.debug('The maxima found are: {}'.format(maxima))
		numbPeaks = len(maxima)
		
		if numbPeaks > 1:
			#The intensities of all the peaks found and stored in maxValues
			maxValues = xk[maxima]
			#The first 'numbPeaks' indexes of the peaks . argsort returns an
			#increasing array. -maxValues is used to obtain an increasing order
			maxIdx = (-maxValues).argsort()[:numbPeaks]
			#The position of the first numbPeaks, ordered by intensity, from highest th lowes
			maxima = maxima[maxIdx]
			if len(maxima)>2:
				maxima = maxima[:2]
		logging.debug('THe maxima considered are: {}'.format(maxima))
		model = lmfit.models.GaussianModel(prefix = 'gauss1_')
		params = model.make_params()
		params['gauss1_center'].set(value = Rvec[maxima[0]])
		params['gauss1_amplitude'].set(value = xk[maxima[0]])
		if len(maxima) > 1:
			model += lmfit.models.GaussianModel(prefix = 'gauss2_')
			params += model.make_params()
			params['gauss2_center'].set(value = Rvec[maxima[1]])
			params['gauss2_amplitude'].set(value = xk[maxima[1]])
		
		result = model.fit(xk, params = params, x = Rvec)
		
		if len(maxima) == 1:
			logging.debug('Performing Single Gaussian Fit')
			paramSugg = {'R_av': result.params['gauss1_center'].value, 'sigma': result.params['gauss1_sigma'].value}
			self.fit_data(fitType = 'Single_Gauss', qRange = qRange, patamSugg=paramSugg,verbose = kwargs.get('verbose',False),fit_kws = kwargs.get('fit_kws',{})) 
			
		if len(maxima) == 2:
			logging.debug('Performing Double Gaussian Fit')
			paramSugg = {'R1_av': result.params['gauss1_center'].value, 'sigma1': result.params['gauss1_sigma'].value}
			paramSugg = {'R2_av': result.params['gauss2_center'].value, 'sigma2': result.params['gauss2_sigma'].value}
			self.fit_data(fitType = 'Double_Gauss', qRange = qRange, patamSugg=paramSugg,verbose = kwargs.get('verbose',False),fit_kws = kwargs.get('fit_kws',{})) 
			
	
	def plot_fit(self, fitType,ax = None, fig = None, **kwargs):
		"""plot_fit: used to plot the results of the fit. Creates 2 axis, on the left the
		data plus the theoretical scattering curve (in solid line in the fitted range and as 
		dashed line outside of the range). On the right plots the particle distribution as
		 a histogram. If the two axis are provided it used the given axes.
			Args:
				fitType (string): the type of fit to plot
				ax (list of :obj: plt.axes): the two axis used to plot the curves
					and the distribution. If not provided they will be created.
					Defaults to None
				**kwargs (dict): the different settings used for the plotting:
					qUnits (string): the units used on the x-axis. they must be defined
						in _AvlbUnits, contained in the config.py file. Defaults to nm
					radRange (list): the range in which to plot the distribution. If
						it doesn't exist the method will try to find a range in which the
						whole peak is visible
			Returns:
				the handles to both axes and the figure
					
		"""
		if kwargs.get('verbose', False):
			print 'logging set to debug'
			logging.basicConfig(format='%(levelname)s:%(message)s', level = logging.DEBUG)
		else:
			print 'logging set to info'
			logging.basicConfig(format='%(levelname)s:%(message)s', level = logging.INFO)
		
		qConvFact = _UnitsConv[kwargs.get('qUnits','nm')]/_UnitsConv[self.qUnits]
		IConvFact = _UnitsConv[kwargs.get('IUnits','m')]/_UnitsConv[self.IUnits]
		
		for k in self.fit_plot_dict:
			if k in kwargs:
				self.fit_plot_dict[k] = kwargs[k]
		#This first set of elif serves to prepare the model and params objects if the
		#required fit was made using lmfit functions. It also serves the purpose
		#of defining the range over which the data was fit, in order to differentiate
		#where the fit was made from where it was extrapolated, and to check that the
		#fit was actually done.
		if fitType in _lmfitModels:
			model = lmfit.Model(getattr(Fitting_Models, _lmfitModelFunctions[fitType]),\
									independent_variables = ['q'])
			distModel = lmfit.Model(getattr(Fitting_Models, _lmfitDistFunctions[fitType]),\
									independent_variables = ['x'])
			if fitType == "Sing_Gauss":
				if self.singleGauss is None:
					logging.error("No fit using a single gaussian has been performed yet!")
					return 0
				else:
					params = self.create_params(fitType)
					for p in params:
						params[p].set(value = self.singleGauss[p])
					RedqRange = self.singleGauss['qRange']
			if fitType == "Double_Gauss":
				if self.fitResults['Double_Gauss'] is None:
					logging.error("No fit using a double gaussian has been performed yet!")
					return 0
				else:
					params = self.create_params(fitType)
					for p in params:
						params[p].set(value = self.fitResults['Double_Gauss'][p])
						#params[p].set(value = self.doubleGauss[p])
					logging.debug("The parameters are:\n {}".format(params))
					print "The parameters are:\n {}".format(params)
					#model = lmfit.Model(getattr(Scattering_Models,fitType), independent_variables = ['q'])
					RedqRange = self.doubleGauss['qRange']
			if fitType == "Sing_Schultz":
				if self.singleSchultz is None:
					logging.error("No fit using a single schultz has been performed yet!")
					return 0
				else:
					params = self.create_params(fitType)
					for p in params:
						params[p].set(value = self.singleSchultz[p])
					#model = lmfit.Model(getattr(Scattering_Models,fitType), independent_variables = ['q'])
					RedqRange = self.singleSchultz['qRange']
			if fitType == "Double_Schultz":
				if self.doubleSchultz is None:
					logging.error("No fit using a double schultz has been performed yet!")
					return 0
				else:
					params = self.create_params(fitType)
					for p in params:
						params[p].set(value = self.doubleSchultz[p])
					#model = lmfit.Model(getattr(Scattering_Models,fitType), independent_variables = ['q'])
					RedqRange = self.doubleSchultz['qRange']
		if fitType == "EM":
			if self.EMFit is None:
				logging.error("No expected maximization has been performed yet!")
				return 0
			else:
				RedqRange = self.EMFit['qRange']
				distModel = 'EM'
				
		#The axes are set-up, the first one to plot the scattering curve+fit, the second 
		#to plot the distribution
		if ax is None or len(ax) != 2:
			fig = Figure(figsize=kwargs.get('figsize',(22,10)))
			ax=[]
			gs = matplotlib.gridspec.GridSpec(1,3)
			ax.append(fig.add_subplot(gs[0,0:2]))
			ax.append(fig.add_subplot(gs[0,2]))
			#ax.append(plt.subplot2grid((1,3), (0,0), rowspan = 1, colspan = 2))
			#ax.append(plt.subplot2grid((1,3), (0,2), rowspan = 1, colspan = 1))
		
		#If the error is available the experimental data is plotted with error
		#bars, if not as normal scatter points 
		
		if self.Ierr is not None and kwargs.get('plot_error',False):			
			ax[0].errorbar(self.q*qConvFact, self.Iabs*IConvFact, yerr=self.Ierr*IConvFact,\
							color = self.fit_plot_dict['color'], linestyle = self.fit_plot_dict['linestyle'],\
							marker = self.fit_plot_dict['marker'], ms = self.fit_plot_dict['ms'],\
							markeredgecolor = self.fit_plot_dict['mec'])
			#color = kwargs.get('expColor','k'), linestyle = kwargs.get('expLinestyle','None'),\
			#marker = kwargs.get('marker','o'), ms = kwargs.get('ms',2), \
			#markeredgecolor = kwargs.get('expColor','k'))
			ax[0].set_xscale("log", nonposx = "clip")
			ax[0].set_yscale("log", nonposx = "clip")
		else:
			ax[0].loglog(self.q*qConvFact, self.Iabs*IConvFact,\
						color = self.fit_plot_dict['color'], linestyle = self.fit_plot_dict['linestyle'],linewidth = self.fit_plot_dict['linewidth'],\
						marker = self.fit_plot_dict['marker'], ms = self.fit_plot_dict['ms'],\
						markeredgecolor = self.fit_plot_dict['mec'])
		
		#The two Schultz distributions require particular attention as they had to be 
		#scaled to m in order to correctly fit.
		tempq = self.q[np.logical_and(self.q>min(RedqRange), self.q<max(RedqRange))]
		
		#If the radii range in which to plot the distribution is not given it has to be calculated.
		#This is done using the _find_distribution_limits function.
		rLimits = kwargs.get('radRange',None)
		if rLimits is None:
			if fitType in _lmfitModels:
				rLimits = self._find_distribution_limits(fitType, distModel,params)
			else:
				rLimits = self._find_distribution_limits(fitType, distModel)
		logging.debug('Calculated limits for radii: {}'.format(rLimits))
		radii = np.arange(min(rLimits), max(rLimits), (max(rLimits)-min(rLimits))/2000)
		
		if fitType in _lmfitModels:
			if fitType == "Sing_Schultz":
				ax[0].loglog(self.q*qConvFact, model.eval(params = params, q=self.q*self.singleSchultz['ConvFact']),\
							color = self.fit_plot_dict['fit_color'], linestyle = '--', linewidth = self.fit_plot_dict['fit_linewidth'],\
							marker = self.fit_plot_dict['fit_marker'], ms = self.fit_plot_dict['fit_ms'],\
							markeredgecolor = self.fit_plot_dict['fit_mec'])
				ax[0].loglog(tempq*qConvFact, model.eval(params = params, q=tempq*self.singleSchultz['ConvFact']),\
							color = self.fit_plot_dict['fit_color'], linestyle = self.fit_plot_dict['fit_linestyle'], linewidth = self.fit_plot_dict['fit_linewidth'],\
							marker = self.fit_plot_dict['fit_marker'], ms = self.fit_plot_dict['fit_ms'],\
							markeredgecolor = self.fit_plot_dict['fit_mec'])
				
				ax[1].plot(radii*self.singleSchultz['ConvFact']/qConvFact, distModel.eval(params = params, x = radii),\
						   color = self.fit_plot_dict['fit_color'], linestyle = self.fit_plot_dict['fit_linestyle'], linewidth = self.fit_plot_dict['fit_linewidth'],\
							marker = self.fit_plot_dict['fit_marker'], ms = self.fit_plot_dict['fit_ms'],\
							markeredgecolor = self.fit_plot_dict['fit_mec'])
			elif fitType == "Double_Schultz":
				ax[0].loglog(self.q*qConvFact, model.eval(params = params, q=self.q*self.doubleSchultz['ConvFact']),\
							color = self.fit_plot_dict['fit_color'], linestyle = '--', linewidth = self.fit_plot_dict['fit_linewidth'],\
							marker = self.fit_plot_dict['fit_marker'], ms = self.fit_plot_dict['fit_ms'],\
							markeredgecolor = self.fit_plot_dict['fit_mec'])
				ax[0].loglog(tempq*qConvFact, model.eval(params = params, q=tempq*self.doubleSchultz['ConvFact']),\
							color = self.fit_plot_dict['fit_color'], linestyle = self.fit_plot_dict['fit_linestyle'], linewidth = self.fit_plot_dict['fit_linewidth'],\
							marker = self.fit_plot_dict['fit_marker'], ms = self.fit_plot_dict['fit_ms'],\
							markeredgecolor = self.fit_plot_dict['fit_mec'])
				
				ax[1].plot(radii*self.doubleSchultz['ConvFact']/qConvFact, distModel.eval(params = params, x = radii),\
						   color = self.fit_plot_dict['fit_color'], linestyle = self.fit_plot_dict['fit_linestyle'], linewidth = self.fit_plot_dict['fit_linewidth'],\
							marker = self.fit_plot_dict['fit_marker'], ms = self.fit_plot_dict['_fitms'],\
							markeredgecolor = self.fit_plot_dict['mec'])
			else:
				logging.debug('Plotting {}'.format(fitType))
				ax[0].loglog(self.q*qConvFact, model.eval(params = params, q=self.q),\
							color = self.fit_plot_dict['fit_color'], linestyle = '--', linewidth = self.fit_plot_dict['fit_linewidth'],\
							marker = self.fit_plot_dict['fit_marker'], ms = self.fit_plot_dict['fit_ms'],\
							markeredgecolor = self.fit_plot_dict['fit_mec'])
				ax[0].loglog(tempq*qConvFact, model.eval(params = params, q = tempq),\
							color = self.fit_plot_dict['fit_color'], linestyle = self.fit_plot_dict['fit_linestyle'], linewidth = self.fit_plot_dict['fit_linewidth'],\
							marker = self.fit_plot_dict['fit_marker'], ms = self.fit_plot_dict['fit_ms'],\
							markeredgecolor = self.fit_plot_dict['fit_mec'])
				
				ax[1].plot(radii/qConvFact, distModel.eval(params = params, x = radii),\
							color = self.fit_plot_dict['fit_color'], linestyle = self.fit_plot_dict['fit_linestyle'], linewidth = self.fit_plot_dict['fit_linewidth'],\
							marker = self.fit_plot_dict['fit_marker'], ms = self.fit_plot_dict['fit_ms'],\
							markeredgecolor = self.fit_plot_dict['fit_mec'])
		elif fitType == 'EM':
			#self._plot_em_result(ax, RedqRange, **kwargs)
			#print 'tempq shape', tempq.shape
			#expI = np.dot(self.EMFit['H'], self.EMFit['xk']).squeeze()
			#print 'H shape', expI.shape
			ax[0].loglog(tempq*qConvFact, np.dot(self.EMFit['H'], self.EMFit['xk'])*IConvFact,\
						color = self.fit_plot_dict['fit_color'], linestyle = self.fit_plot_dict['fit_linestyle'], linewidth = self.fit_plot_dict['fit_linewidth'],\
						marker = self.fit_plot_dict['fit_marker'], ms = self.fit_plot_dict['fit_ms'],\
						markeredgecolor = self.fit_plot_dict['fit_mec'])
			#The H matrix has to be recreated to calculate the whole q range
			Vvec = self.fitResults['EM']['Rvec']**3*4/3*np.pi
			FF = EM.sphere_form_factor(self.q,self.fitResults['EM']['Rvec'])
			H = FF*np.matlib.repmat(Vvec,len(self.q),1)
			
			ax[0].loglog(self.q*qConvFact, np.dot(H, self.fitResults['EM']['xk'])*IConvFact,\
						color = self.fit_plot_dict['fit_color'], linestyle = '--', linewidth = self.fit_plot_dict['fit_linewidth'],\
						marker = self.fit_plot_dict['fit_marker'], ms = self.fit_plot_dict['fit_ms'],\
						markeredgecolor = self.fit_plot_dict['fit_mec'])
			
			Volume = np.array(4./3.*self.fitResults['EM']['Rvec']**3*np.pi)
			#np.array([x/v for x,v in zip(xk,Volume)])
			sizeDist = np.array([x/v for x,v in zip(self.fitResults['EM']['xk'],Volume)])
			numbFrac = sizeDist
			#numbFrac = self.EMFit['xk'].flatten()/self.EMFit['xk'].flatten().sum()
			totVolume = self.fitResults['EM']['xk'].sum()
			#totVolume = np.sum(Volume*self.EMFit['xk'][:,0])
			#numbDist = np.array([x/v for x,v in zip(xk,Volume)])

			#volFrac = (Volume*self.EMFit['xk'][:,0])/totVolume
			volFrac = (self.EMFit['xk'][:,0])/totVolume
			print kwargs
			if kwargs.get('volumeDist', False):
				
				ax[1].plot(self.EMFit['Rvec']/qConvFact, volFrac,\
							color = self.fit_plot_dict['fit_color'], linestyle = self.fit_plot_dict['fit_linestyle'], linewidth = self.fit_plot_dict['fit_linewidth'],\
							marker = self.fit_plot_dict['fit_marker'], ms = self.fit_plot_dict['fit_ms'],\
							markeredgecolor = self.fit_plot_dict['fit_mec'])
			else:
				ax[1].plot(self.EMFit['Rvec']/qConvFact, numbFrac,\
							color = self.fit_plot_dict['fit_color'], linestyle = self.fit_plot_dict['fit_linestyle'], linewidth = self.fit_plot_dict['fit_linewidth'],\
							marker = self.fit_plot_dict['fit_marker'], ms = self.fit_plot_dict['fit_ms'],\
							markeredgecolor = self.fit_plot_dict['fit_mec'])
		else:
			logging.info("Plotting of other fitting methods has not been implemented yet")
			return 0
		#plt.show()
		
		return (ax,fig)
				
	def _find_distribution_limits(self, fitType, model, params = None):
		"""_find_distribution_limits: determins the limits of the distribution plot
		if they are not given. This is done by taking the peak (for bimodal dist.
		the left peak is used for the left limit while the right peak is used for
		the right limit), and calculating the value of the distribution.
		from 0 to the peak. The left boundary will be taken as the value at which the
		distribution reaches 1/n the value at the peak (n=100?), or 0. The right boundary
		will be taken by calculating the distribution between the rightmost peak and a
		position simmetrically opposit the left peak. The limit will be selected as the
		value at which the distribution reaches 1/n the value at the peak (n=100?) or the
		symmetric position.
			Args:
				fitType (string): the type of fit to plot
				model (lmfit.Model): the lmfit model used to calculate the distribution
				params (lmfit.Parameters): the lmfit parameters object which contains
					the parameters of the distribution
			Returns:
				Numpy.array witht he two values in increasing order
		"""
		if fitType in _lmfitModels:
			if params is None:
				print "To calculate the radii limits for an lmfit model you need to provide the prameters"
				return 0
		
			if "Sing" in fitType:
				R_av = params['R_av'].value
			elif "Double" in fitType:
				R_av = min(params['R1_av'].value, params['R2_av'].value)
				print 'The largest peak is at: ', R_av
			xs = np.arange(0,R_av,R_av/100.)
			#print xs
			ys = model.eval(params = params, x = xs)
			#print ys
			#print 'xs: ',xs
			#print 'params: ',params
			#print 'ys: ',ys
			
			lowerThan = xs[ys<ys[-1]/100.]
			print 'lower than: ', lowerThan
			if lowerThan.size == 0:
				lLimit = 0
			else:
				lLimit = lowerThan[-1]
			if "Double" in fitType:
				R_av = max(params['R1_av'].value, params['R2_av'].value)
				print "The second peak is at: ", R_av
			xs = np.arange(R_av,2*R_av, R_av/100.)
			ys = model.eval(params = params, x = xs)
			lowerThan = xs[ys<ys[0]/100.]
			
			#fig1 = plt.figure()
			#ax1 = fig1.add_subplot(111)
			#ax1.plot(xs,ys)
			#ax1.plot([xs[0],xs[-1]],[ys[0]/100.,ys[0]/100.])
			#plt.show()
			
			if lowerThan.size == 0:
				rLimit = 0
			else:
				rLimit = lowerThan[0]
		elif fitType == 'EM':
			#Should be changed so that the maxima are found, the first and secon
			lLimit = self.EMFit['Rvec'][0]
			rLimit = self.EMFit['Rvec'][-1]
			
		
		return np.array([lLimit,rLimit])
	
	def fit_guinier_regime(self, error_lim = 0.01, min_q = 3,skip_initial = 0, plot_guinier = False):
		"""guinier_regime: fit the data in the guinier regime using the linear relationship
			Args:
				error_lim (long): acceptable relative error between the linear fit and the
					data. the value has to be between 0 and 1, 0 excluded. Defaults to 0.01
				min_q (int): minimum number of points needed to validate the Guinier regime
					fit. Has to be greater than 3. Defaults to 3
				skip_initial (int): number of initial points to ignore
				plot_guinier (bool): whether to plot the guinier fit or not. Defaults to
					False
		"""
		if error_lim == 0 or error_lim > 1:
			print 'The provided error limit {}, is not a valid value. The error limit has\
			to be included within (0,1]'.format(error_lim)
			return 0
		if min_q < 3:
			print 'In order for the fit to be valid there have to be at least 3 points\
			over which to fit. The given number of points is {}'.format(min_q)
			return 0
			
		q_2 = self.q[skip_initial:]**2
		I_log = np.log(self.Iabs[skip_initial:])
		linMod = lmfit.models.LinearModel()
		curr_q = min_q
		#guinierRegime = True
		advance = True
		while advance:
			params = linMod.guess(I_log[:curr_q], x = q_2[:curr_q])
			result = linMod.fit(I_log[:curr_q],params = params, x = q_2[:curr_q])
			Rg = np.sqrt(-result.params['slope'].value*3)
			error = np.mean(result.residual/I_log[:curr_q])
			guinierRegime = Rg*self.q[(skip_initial+curr_q)]<1
			if error>error_lim or not guinierRegime:
				advance = False
				curr_q -=1
			else:
				curr_q +=1
			
		if curr_q > min_q:
			#Recalculate the fit of the last valid q
			params = linMod.guess(I_log[:curr_q], x = q_2[:curr_q])
			result = linMod.fit(I_log[:curr_q],params = params, x = q_2[:curr_q])
			Rg = np.sqrt(-result.params['slope'].value*3)
		if plot_guinier:
			fig = Figure(figsize = (6,6))
			ax = fig.add_subplot(111)
			ax.plot(q_2[:(2*curr_q)], I_log[:(2*curr_q)], 'k', marker = 'o', linestyle = 'None')
			ax.plot(q_2[:curr_q], result.eval(x = q_2[:curr_q]), 'r', linewidth = 2)
			ax.plot(q_2[:(2*curr_q)], result.eval(x = q_2[:(2*curr_q)]), 'r', linestyle = '--',linewidth = 2)
			
		self.guinier_regime = [result, curr_q, self.q[curr_q], guinierRegime]
		return {'best_fit': result,'last_q_pos': skip_initial+curr_q, \
				'last_q': self.q[(skip_initial+curr_q)],'within_Guinier': guinierRegime,\
				'initial_points_skipped': skip_initial, 'Rg': Rg}

	def find_porod_invariant(self, drho, units = 'm', qRange = [0,np.inf], precision = 0.1, postExtrap = None,\
						preExtrap = None, plot_porod = False):
		"""porod_invariant: calculate the porod invariant for the given data. The porod
		invariant is calculated as the integral of I*q^2. The data is integrated, and
		then extrapolated to infinity by giving a constant slope to the signal of the
		type I(q) = q^(slope). The complete equation is:
		Q = pi**2*phi*(1-phi)*(drho)**2
		for simplicity, as the systems have a low volume fraction, the 1-phi is simplified
		to 1. therefore
		phi = Q/(2*np.pi**2*drho**2)
		
			Args:
				drho (int): the difference in scattering length density between the two
					phases
				units (string): the units used for the drho. The intensity will be converted
					to the right units. Defaults to 'm'
				qRange (list): the q range to use for the experimental data. If None is give
					than the whole range is considered. Defaults to [0,np.inf]
				precision (int): the precision used to calculate the spline to interpolate
					the data points. Defaults to 0.01
				postExtrap (list of dict): list of regimes to use to extrapolate to large
					q values. Each dictionary must have a slope (automatically considered
					negative), and the two extremes between which to integrate within the
					regeme. If the last regeme is steeper than -3 then it will be used to
					infinity. If not -4 will be used. If None it will interpolate using
					-4 at the end of the experimenta data. Defaults to None
				preExtrap (int): the slope to use to interpolare at low q. If None is give
					then the interpolation will be constant. Defaults to None
		"""
		tempq = self.q[np.logical_and(self.q>min(qRange), self.q<max(qRange))]
		tempI = self.Iabs[np.logical_and(self.q>min(qRange), self.q<max(qRange))]
		Iq2 = tempI * tempq**2
		
		if plot_porod:
			fig = Figure(figsize = (6,6))
			ax = fig.add_subplot(111)
		#s is the precision used to calculate the spline
		s = np.average((Iq2*precision)**2)
		f_exp = sp.interpolate.UnivariateSpline(tempq, Iq2, k = 3, s = s)
		#experimental integral
		Q = f_exp.integral(tempq[0],tempq[-1])
		
		if plot_porod:
			ax.plot(tempq,Iq2)
		#the general extrapolation function is of the type:
		# I*q^2 = a*x^(slope)*x**2
		extrap = lambda x,a,slope: a*x**(slope)*x**2
		#the integral is going to be like the function but with
		#and extra *x
		#extrapInt = lambda x,a,slope: a*x**(slope)*x**2*x
		if len(postExtrap) == 0:
			#first the coefficient to align the data is calculated
			I_ratio = Iq2[-1]/extrap(tempq[-1],1,-4)
			Q += sp.integrate.quad(extrap,tempq[-1],np.inf, args = (I_ratio, -4))[0]
			if plot_porod:
				qs = np.array([tempq[-1],2*tempq[-1]])
				ax.plot(qs, estrap(qs,I_ratio,-4))
		else:
			#prevSlope
			last_q = tempq[-1]
			last_I = Iq2[-1]
			for i in xrange(len(postExtrap)):
				#get values form the dictionary
				curr_slope = -abs(postExtrap[i].get('slope',-4))
				max_q = postExtrap[i].get('max_q', np.inf)
				
				if np.isinf(max_q) and curr_slope > -3:
					print 'The final regime cannot decrease slower than -3 (given{}), the integral will\
					not converge.'. format(curr_slope)
					return 0
				I_ratio = last_I/extrap(last_q,1,curr_slope)
				Q += sp.integrate.quad(extrap,last_q,max_q, args = (I_ratio, curr_slope))[0]
				
				if plot_porod:
					if np.isinf(max_q):
						qs = np.array([last_q,2*last_q])
					else:
						qs = np.array([last_q,max_q])
					ax.plot(qs, extrap(qs,I_ratio,curr_slope))
				
				last_q = max_q
				last_I = extrap(last_q,I_ratio,curr_slope)
				
				
				if np.isinf(max_q):
					break
			#If the extrapolation list did not end with a regime up to infinity
			#then the signal is extended with slopw -4
			if not np.isinf(last_q):
				I_ratio = last_I/extrap(last_q,1,-4)
				Q += sp.integrate.quad(extrap,last_q,np.inf, args = (I_ratio, -4))[0]
				if plot_porod:
					qs = np.array([last_q,2*last_q])
					ax.plot(qs, estrap(qs,I_ratio,-4))
			
		if preExtrap is None:
			#if no pre-regeme is given then the data is integrated as
			#a rectangle.
			Q += Iq2[0]*tempq[0]
		else:
			I_ratio = Iq2[0]/extrap(tempq[0],1,preExtrap)
			Q += sp.integrate.quad(extrap,0,tempq[0], args = (I_ratio, preExtrap))[0]
		
		self.porod_invariant = Q/(2*np.pi**2*drho**2)
		return Q/(2*np.pi**2*drho**2)
	
	def get_fit_params(self, fitType, param = None, **kwargs):
		'''get_fit_params returns a dictionary with the various fitting parameters as keys and
		the fit as values.
			Args:
				fitType (String): the fit type for which the parameters are required.
					can be either Sing_Gauss, Double_Gauss, Sing_Schultz, Double_Schultz
					or EM
				param (list: String): if a parameter name is specified it will return the specific
					parameter.
			Return: returns the dictionary containing the fitting parameters and their values
		'''
		if kwargs.get('verbose', False):
			logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
		else:
			logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
			
		if fitType == 'Sing_Gauss':
			if self.singleGauss is None:
				return None
			elif param is None:
				return self.singleGauss
			else:
				if isinstance(param, list):
					return {p: self.singleGauss.get(p, None) for p in params}
				elif isinstance(param, str):
					return {param: self.singleGauss.get(param, None)}
				else:
					logging.error('The format in which the parameter was passed is not recognizable: \n {} \n Please provide either a string or a list of strings'.format(param))
					
		elif fitType == 'Double_Gauss':
			if self.doubleGauss is None:
				return None
			elif param is None:
				return self.doubleGauss
			else:
				if isinstance(param, list):
					return {p: self.doubleGauss.get(p, None) for p in params}
				elif isinstance(param, str):
					return {param: self.doubleGauss.get(param, None)}
				else:
					logging.error('The format in which the parameter was passed is not recognizable: \n {} \n Please provide either a string or a list of strings'.format(param))
					
		elif fitType == 'Sing_Schultz':
			if self.singleSchultz is None:
				return None
			elif param is None:
				return self.singleSchultz
			else:
				if isinstance(param, list):
					return {p: self.singleSchultz.get(p, None) for p in params}
				elif isinstance(param, str):
					return {param: self.singleSchultz.get(param, None)}
				else:
					logging.error('The format in which the parameter was passed is not recognizable: \n {} \n Please provide either a string or a list of strings'.format(param))

		elif fitType == 'Double_Schultz':
			if self.doubleSchultz is None:
				return None
			elif param is None:
				return self.doubleSchultz
			else:
				if isinstance(param, list):
					return {p: self.doubleSchultz.get(p, None) for p in params}
				elif isinstance(param, str):
					return {param: self.doubleSchultz.get(param, None)}
				else:
					logging.error('The format in which the parameter was passed is not recognizable: \n {} \n Please provide either a string or a list of strings'.format(param))

		elif fitType == 'EM':
			if self.EMFit is None:
				return None
			elif param is None:
				return self.EMFit
			else:
				if isinstance(param, list):
					return {p: self.EMFitz.get(p, None) for p in params}
				elif isinstance(param, str):
					return {param: self.EMFit.get(param, None)}
				else:
					logging.error('The format in which the parameter was passed is not recognizable: \n {} \n Please provide either a string or a list of strings'.format(param))

		else:
			print 'Fit Type not recognized'
	
	#def fit_avlb(self, fitType):
		#if fitType in _AvlbSASFit:
		