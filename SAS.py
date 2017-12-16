from config import _AvlbUnits, _UnitsSymbols,_UnitsConv, _AvlbSASFit, _AvlbSASFitDic, _AvlbSASFitDicInv, _lmfitModels, _lmfitModelFunctions, _lmfitDistFunctions
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
from matplotlib.ticker import ScalarFormatter
from matplotlib.figure import Figure
from matplotlib import rc
rc('text', usetex=True)
import scipy as sp
from scipy import signal
from scipy import interpolate
from collections import OrderedDict
from mpl_toolkits.mplot3d import Axes3D
try :
	import lmfit
except ImportError:
	print "Could not import the lmfit package. Some fittings will be disabled."
	lmfitAvlb = False
else:
	lmfitAvlb = True
	import Fitting_Models

try:
    import gmpy2
    from gmpy2 import mpz,mpq,mpfr,mpc
except:
    gmpy2 = None
    from decimal import Decimal
    print 'The Schultz fitting function is optimized by using the GMPY2 module to\
deal with the large numbers required. The module was not found, so Decimal will\
be used instead, but the calculations will be slower.'



class SmallAngleScattering(ScatteringObject):
	"""SmallAngleScattering: this class is a child of ScatteringObject used to store
	data arriving from small angle scattering measurements. It provides a method
	to plot the data in log log scale as well as to fit it. There are several
	different fitting methods available. The different methods can be divided in
	2 main categories: those using the lmfit library and the one based on the
	expected maximization method. The former are not available if the lmfit module
	was not correctly loaded. It is also possible to first fit using the EM
	method, and later use it as a starting point for the lmfit fit.
	"""
	def __init__(self,fname = None,**kwargs):
		'''The results of the fit are store in a dictionary, which is initaiate
		with None values using the available fits described in _AvlbSASFit. Two
		other dictionaries are created ot store the result from the porod
		invariant and the guinier regime
		'''
		self.fitResults = {f : None for f in _AvlbSASFit}
		self.porod_invariant = None
		self.guinier_regime = None

		super(SmallAngleScattering, self).__init__(fname, **kwargs)

	def plot_data(self, qRange=[0,np.inf],yShift=1, ax = None, fig = None, figSize = (6,6),\
	 			  xUnits = None, yUnits =None, plot_dic={}, axLabelSize = 20,\
				  labelSize = 18, tickSize = 6, tickWidth = 3, **kwargs):
		"""plot_data: Plots the data in a log-log plot. If the axis is passed
		then the data is plotted on the given axis, if not a figure of size
		figSize is created and the axis places in it. The units on the x and y
		axis can be selected. The keyword arguments are used to set the
		parameters of the plots.
			Args:
				-qRange (list): the q-range over which the data should be
					plotted. Defaults to [0, np.inf]
				-yShift (int): indicated by how much the data should be shifted
					vertically. This allows the plot_data function to be used in
					time resolved data. The data is multiplied by shift because
					it is plotted in loglog. Defaults to 1
				-ax (:obj: plt.axes): the axis on which to plot the data. If it
					is not provided then one will be created. Defaults to None
				-fig (:obj:plt.figure): the figure in which the axis is placed
				-figSize (tuple): used to set the size of the figure in case the
					axis is not provided. Defaults to (6,6)
				-xUnits (string): the units used on the x-axis. They have to be
					defined in the _AvlbUnits in the config.py file. Defaults to
					nm
				-yUnits (string): the units used on the xy-axis. They have to be
					defined in the _AvlbUnits in the config.py file. Defaults to
					m
				-plot_dic(dict): any keyword used by matplotlib in plot can be
					placed here to influence the appearance of the plot
				**kwargs: all the arguments used to set the values of the plot.
					-axLabelSize (int): the size of the font used for the lable
						of the x and y axis ('Wavelength' and 'Intensity')
					-labelSize (int): the size of the numbers along the x and y
						axis.
					-tickSize (int): the length of the major tickmarks on the x
						and y axis
					-tickWidth (int): the width of the major tickmarks on the x
						and y axis
			Returns:
				ax (:obj: plt.axes): the axis on which the data was plotted
				fig (:obj: plt.fig): the figure on which the data was plotted
		"""

		if kwargs.get('verbose', False):
			logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
		else:
			logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
		#Any keywords which might alter the plotting parameters are stored
		#in the dictionary. This also allows the next plot to retain the
		#same plotting style.
		for k in self.plot_dict:
			if k in plot_dic:
				self.plot_dict[k] = plot_dic[k]
		#plt.ion()
		if ax is None:
			fig = Figure(figsize = figSize)
			ax = fig.add_subplot(111)
		else:
			if fig is None:
				fig = ax.get_figure()
		mask = np.logical_and(self.q>min(qRange), self.q<max(qRange))
		tempq = self.q[mask]
		tempI = self.Iabs[mask]
		if self.Ierr is not None:
			tempE = self.Ierr[mask]

		#conversion coefficients are calculated in order to plot the graph
		#in the correct scale
		qConvCoeff = 1
		IConvCoeff = 1
		if (xUnits in _AvlbUnits) and (self.qUnits in _AvlbUnits):
			qConvCoeff = _UnitsConv[self.qUnits]/_UnitsConv[xUnits]
		else:
			xUnits = self.qUnits

		if (yUnits in _AvlbUnits) and (self.IUnits in _AvlbUnits):
			IConvCoeff = _UnitsConv[self.IUnits]/_UnitsConv[yUnits]
		else:
			yUnits = self.IUnits
		logging.debug('SAS.plot_data():\nqConvCoeff: {}\nIConvCoeff: {}'.format(qConvCoeff,IConvCoeff))

		if (self.Ierr is not None) and (kwargs.get('plot_error', False)):
			'''
			This part was modified, instead of asking the values every time, a
			plot dictionary was created to store all the properties, and is
			updated in case new parameters are passed to the plotting function.
			'''
			ax.errorbar(tempq*qConvCoeff,tempI*IConvCoeff*yShift,yerr = tempE*IConvCoeff*yShift, **self.plot_dict)
					#color = kwargs.get('color','k'),\
		 			#linestyle = kwargs.get('linestyle', 'None'),marker = kwargs.get('marker','o'),\
		 			#markeredgecolor = kwargs.get('color','k'), ms = kwargs.get('ms',2))
			ax.set_xscale("log", nonposx = "clip")
			ax.set_yscale("log", nonposx = "clip")
		else:
			logging.debug('SAS.plot_data():\nplotting q: {}\nplotting I: {}'.format(tempq*qConvCoeff,tempI[10:]*IConvCoeff*yShift))
			ax.loglog(tempq*qConvCoeff,tempI*IConvCoeff*yShift, **self.plot_dict)
					  #color = kwargs.get('color','k'),\
						#linestyle = kwargs.get('linestyle', 'None'),marker = kwargs.get('marker','o'),\
						#markeredgecolor = kwargs.get('color','k'), ms = kwargs.get('ms',2))

		#Set the x an y limits
		#yLimits = [min(tempI[np.isfinite(tempI)]*IConvCoeff*yShift),max(tempI[np.isfinite(tempI)]*IConvCoeff*yShift)]
		ax.set_ylim(min(tempI[np.isfinite(tempI)]*IConvCoeff*yShift),max(tempI[np.isfinite(tempI)]*IConvCoeff*yShift))
		ax.set_xlim(min(tempq[np.isfinite(tempI)]*qConvCoeff),max(tempq[np.isfinite(tempI)]*qConvCoeff))

		#Set the labels for the 2 axis and the tick parameters
		if xUnits in _AvlbUnits:
			ax.set_xlabel(r'Wavevector ({}$^{{-1}}$)'.format(_UnitsSymbols[xUnits]), fontsize = kwargs.get('lableSize',20))
		else:
			ax.set_xlabel(r'Wavevector', fontsize = axLabelSize)
		if yUnits in _AvlbUnits:
			ax.set_ylabel(r'Intensity ({}$^{{-1}}$)'.format(_UnitsSymbols[yUnits]), fontsize = kwargs.get('lableSize',20))
		else:
			ax.set_ylabel(r'Intensity', fontsize = axLabelSize)
		ax.tick_params(labelsize = labelSize, size = tickSize,
		 				width = tickWidth)
		#Format The Scientific notation
		#ax.ticklabel_format(style='sci',scilimits=(-2,2),axis='y')
		#implement a method to place the major and minor tick marks
		return (ax, fig)

	def plot_slopes(self, slope, frequency = 10,slope_color ='k', ax = None, **kwargs):
		"""plot_slopes: plots the data along with a series of parallel lines with given slope.
		this can be useful when trying ot detrmine the slope of a particlare scattering curve.
			Args:
				slope (int): the power of the polynomial to plot x^(slope). As the slope is
					(almost) always decreasing, any value given will be taken as negative
				frequency (int): the distance between subsequent lines. The nth line will
					therefore be frequency^(n)*x^slope. Defaults to 10
				slope_color (char): a valid character to define the color with which the
					sloped lines should be drawn. Defaults to 'k' (black)
				ax (obj:axis): the axis on which to draw the slopes
				**kwargs (dict): any keywords which can be used with the plot_data function
			Retuns:
				the axis and the figure of the plot
		"""
		slope = -abs(slope)

		#If no axis is given it is supposed that the slopes have to be drawn
		#over the data from the object, therefore the data will be drawn and the
		#lines added
		if ax is None:
			ax, fig = self.plot_data(**kwargs)
		##qRange could be substituted by ax.limits
		qRange = kwargs.get('qRange',[0,np.inf])
		mask = np.logical_and(self.q>min(qRange), self.q<max(qRange))
		tempq = self.q[mask]
		tempI = self.Iabs[mask]
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
		return [ax,fig]

	def fit_data(self, fitType = "Sing_Gauss",**kwargs):
		"""fit_data: wrapper function used to decide which fitting algorithm to
		use and to set it up.
			Args:
				-fitType (string): the method used for fitting the data. Has to
					be defined in _AvlbSASFit contained in config.py. Defaults
					to Sing_Gauss.
				**kwargs (dict): The keywords to pass to the desired fitting method.
					-plot (bool): whether to plot or not the fit.
						Defaults to False
					-paramSugg (dict): dictionary passed on to the create_params
						function. It containes any suggestion used to create the
						parameters

		"""

		if kwargs.get('verbose', False):
			logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
		else:
			logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
		#Check to make sure that if the fitting method required used the lmfit
		#package, it is in fact available
		if fitType in _lmfitModels:
			if not lmfitAvlb:
				logging.error("Impossible to fit using {} model because lmfit was not importent.".format(fitType))
				return 0

			params = self.create_params(fitType, kwargs.get('paramSugg', {}))
			self.lmfit_fit(fitType, params, **kwargs)

		if fitType == "EM":
			self.em_fit( **kwargs)

		if fitType == "Dist_from_EM":
			self.lmfit_from_em(**kwargs)

		#if kwargs.get('plot',False):
			#self.plot_fit(type)

		if fitType not in _AvlbSASFit:
			print "{} is not an implemented fit. The available options are:".format(fitType)
			for ft in _AvlbSASFit:
				print "\t - {}".format(ft)
			return 0

	def create_params(self, fitType, paramSugg = {}):
		"""create_params: creates the Parameter object for the lmfit package
		based on the type of fit which is chosen between the available fits.
			Args:
				fitType (string): the method used for fitting. Has to be defined
					in _lmfitModels
				paramSugg (dict): dictionary containing all the default values
					for the parameters. These depend on the fit used.
			Returns:
				params (:obj: Parameters): a Parameters object from the lmfit
				 	package with all the parameters needed to fit using the
					selected model.

		"""
		params = lmfit.Parameters()
		if fitType == "Sing_Gauss":
			params.add("R_av", value = paramSugg.get('R_av',10),\
						min = paramSugg.get('R_min',0),\
					 	max = paramSugg.get('R_max',None),\
						vary = paramSugg.get('R_vary',True))
			params.add("sigma", value = paramSugg.get('sigma',0.5),\
			 			min = paramSugg.get('sigma_min',0),\
					 	max = paramSugg.get('sigma_max',None),\
						vary = paramSugg.get('sigma_vary',True))


		elif fitType == "Double_Gauss":
			params.add("R1_av", value = paramSugg.get('R1_av',5),\
						min = paramSugg.get('R1_min',0),\
					 	max = paramSugg.get('R1_max',None),\
						vary = paramSugg.get('R1_vary',True))
			params.add("sigma1", value = paramSugg.get('sigma1',0.5),
			 			min = paramSugg.get('sigma1_min',0),\
					 	max = paramSugg.get('sigma1_max',None),\
						vary = paramSugg.get('sigma1_vary',True))
			params.add("R2_av", value = paramSugg.get('R2_av',10),\
						min = paramSugg.get('R2_min',0),\
					 	max = paramSugg.get('R2_max',None),\
						vary = paramSugg.get('R2_vary',True))
			params.add("sigma2", value = paramSugg.get('sigma2',0.5),\
						min = paramSugg.get('sigma2_min',0),\
					 	max = paramSugg.get('sigma2_max',None),\
						vary = paramSugg.get('sigma2_vary',True))
			params.add("ratio", value = paramSugg.get('ratio',0.5),\
						min = paramSugg.get('ratio_min',0),\
					 	max = paramSugg.get('ratio_max',1),\
						vary = paramSugg.get('ratio_vary',True))


		elif fitType == "Sing_Schultz":
			params.add("R_av", value = paramSugg.get('R_av',10),\
						min = paramSugg.get('R_min',0),\
					 	max = paramSugg.get('R_max',None),\
						vary = paramSugg.get('R_vary',True))
			params.add("Z", value = paramSugg.get('Z',10),\
						min = paramSugg.get('Z_min',0),\
					 	max = paramSugg.get('Z_max',None),
						vary = paramSugg.get('Z_vary',True))


		elif fitType == "Double_Schultz":
			params.add("R1_av", value = paramSugg.get('R1_av',5),\
						min = paramSugg.get('R1_min',0),\
					 	max = paramSugg.get('R1_max',None),\
						vary = paramSugg.get('R1_vary',True))
			params.add("Z1", value = paramSugg.get('Z1',10),\
						min = paramSugg.get('Z1_min',0),\
					 	max = paramSugg.get('Z1_max',None),\
						vary = paramSugg.get('Z1_vary',True))
			params.add("R2_av", value = paramSugg.get('R2_av',10),\
						min = paramSugg.get('R2_min',0),\
					 	max = paramSugg.get('R2_max',None),\
						vary = paramSugg.get('R2_vary',True))
			params.add("Z2", value = paramSugg.get('Z2',10),\
						min = paramSugg.get('Z2_min',0),\
					 	max = paramSugg.get('Z2_max',None),\
						vary = paramSugg.get('Z2_vary',True))
			params.add("ratio", value = paramSugg.get('ratio',0.5),\
						min = paramSugg.get('ratio_min',0),\
					 	max = paramSugg.get('ratio_max',1),\
						vary = paramSugg.get('ratio_vary',True))

		params.add("I0", value = paramSugg.get('I0',1),\
					min = paramSugg.get('I0_min',None),\
					max = paramSugg.get('I0_max',None),\
					vary = paramSugg.get('I0_vary',True))
		params.add("bckg", value = paramSugg.get('bckg',0),\
					min = paramSugg.get('bckg_min',None),\
					max = paramSugg.get('bckg_max',None),\
					vary = paramSugg.get('bckg_vary',False))

		return params

	def lmfit_fit(self, fitType, params, qRange = np.array([0,np.inf]),\
	 			RedqRange = None, plotFit = False, fit_kws = {},store_history=True, **kwargs):
		"""lmfit_fit: fits the data to one of the model types. kwargs should
		contain 2 qRanges, a first one over which the fit is done over the wide
		qrange where the values are normalized by 1/I. A second at low-q range
		for fitting the I0. If no ranges are provided the whole q range will be
		used for the first fit and 10% of range will be used for the second fit
			Args:
				-fitType (string): the type of distribution used to fit the data
				-params (:obj: Parameters): The Parameter object from lmfit with
					all the parameters provided.
				-qRange (list): the first range of q-values used for fitting with
					weights = 1/I. Defaults to [0,np.inf]
				-RedqRange (list): The second range of q-values used for fitting
					the I0. Defaults to [0,0.1*max(q)].
				-plotFit (bool): decides whether or not to plot the fitted data.
					Defaults to False
				-fit_kws (dict): contains the dictionary to pass over to the
					scipy function for the fit. Defaults ot empty dictionary.
				-store_history (bool): determins whether to store the errors
					during the itertions.
		"""
		if kwargs.get('verbose', False):
			logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
		else:
			logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

		model = lmfit.Model(getattr(Fitting_Models, _lmfitModelFunctions[fitType]),\
							independent_variables = ['q'])

		#Initial fit over the wider qRange
		#qRange = kwargs.get('qRange',np.array([0,np.inf]) )
		mask = np.logical_and(self.q>min(qRange), self.q<max(qRange))
		tempq = self.q[mask]
		tempI = self.Iabs[mask]
		tempI = np.nan_to_num(tempI)

		#convFact = 1
		#if "Schultz" in fitType:
			#logging.info(" Due to the form of the Schultz distribution the data has to be shifted to m^(-1).\
 		#Therefore the intensity fit might be off. Not sure, will have to check.")
			#if self.qUnits is None:
				#logging.info('The wavevecto has no units, nanometers will be used')
				#convFact = _UnitsConv['m']/_UnitsConv['nm']
			#else:
				#convFact = _UnitsConv['m']/_UnitsConv[self.qUnits]
		#convFact = 1
		logging.debug('SAS.lmfit_fit() tempq size: ', tempq.shape)
		tempq = tempq.astype(np.float64)
		tempI = tempI.astype(np.float64)

		#np.array([II.astype(np.float64) for II in tempI])
		#fit_kws = kwargs.get('fit_kws',{})
		if 'nan_policy' not in fit_kws:
			fit_kws['nan_policy'] = 'omit'
		if kwargs.get('verbose', False):
			print 'Params before fit:'
			for k in params:
				print '{}: {}'.format(k, params[k].value)


		if store_history:
			error_history = []
			def storeError(pars, iter, resid, *args, **kwargs):
				error_history.append(resid)
			result = model.fit(tempI, params = params, q=tempq, weights = 1./tempI, fit_kws = fit_kws, iter_cb = storeError)
			logging.debug('_____________FIRST FIT____________\n The fit succeded: {}\n Lmfit message: {}\n The scipy return code is: {}\n The scipy message is: {}\n _______________________'.format(result.success,\
																																									   result.message,\
																																									   result.ier,\
																																							   result.lmdif_message))
			#self.fitResults[fitType]['error_history']=error_history
			

		else:
			result = model.fit(tempI, params = params, q=tempq, weights = 1./tempI, fit_kws = fit_kws)
			logging.debug('_____________FIRST FIT____________\n The fit succeded: {}\n Lmfit message: {}\n The scipy return code is: {}\n The scipy message is: {}\n _______________________'.format(result.success,\
																																									   result.message,\
																																									   result.ier,\
																																									    result.lmdif_message))


		if kwargs.get('verbose',False):
			print result.fit_report()
			print 'Params after first fit:'
			for k in result.params:
				print '{}: {}'.format(k, result.params[k].value)
		#Second fit done over the reduced range in order to correctly fit the
		#intensity
		if RedqRange is None:
			RedqRange = np.array([0,0.1*max(self.q)])
		mask = np.logical_and(self.q>min(RedqRange), self.q<max(RedqRange))
		tempq = self.q[mask]
		tempI = self.Iabs[mask]

		#if "Schultz" in fitType:
			#convFact = _UnitsConv['m']/_UnitsConv[self.qUnits]
		#tempq = tempq*convFact
		#All parameters except I0 are set to constant
		for p in result.params:
			if p != 'I0':
				result.params[p].set(vary = False)
		tempI = np.nan_to_num(tempI)
		tempq = tempq.astype(np.float64)
		tempI = tempI.astype(np.float64)
		#np.array([II.astype(np.float64) for II in tempI])

		resultI = model.fit(tempI, result.params, q = tempq, fit_kws = fit_kws)
		if kwargs.get('verbose', False):
			print 'Params from Second fit:'
			for k in resultI.params:
				print '{}: {}'.format(k, resultI.params[k].value)

		logging.debug('_____________SECOND FIT____________\n The fit succeded: {}\n Lmfit message: {}\n The scipy return code is: {}\n The scipy message is: {}\n _______________________'.format(result.success,\
																																									   result.message,\
																																									   result.ier,\
																																									   result.lmdif_message))
		#The results of the fit are stored in the correct dictionary
		#self.fitResults[fitType] = {p:result.params[p].value for p in result.params}
		self.fitResults[fitType] = result.params.valuesdict()
		self.fitResults[fitType]['I0'] = resultI.params['I0'].value
		self.fitResults[fitType].update({'{}_stderr'.format(p):result.params[p].stderr for p in result.params})
		self.fitResults[fitType]['I0_stderr'] = resultI.params['I0'].stderr
		self.fitResults[fitType]['qRange'] = qRange
		self.fitResults[fitType]['redchi'] = result.redchi
		self.fitResults[fitType]['chi'] = result.chisqr
		self.fitResults[fitType]['residual'] = result.residual
		self.fitResults[fitType]['Units'] = self.qUnits
		if store_history:
			self.fitResults[fitType]['error_history']=error_history

		#print 'The data for {} fit was updated to:\n{}'.format(fitType,self.fitResults[fitType])

		if kwargs.get('plotFit',False):
			self.plot_fit(fitType)

	def em_fit(self, qRange = np.array([0,np.inf]),Rmin = 0.1, Rmax = 0, numbElements = 100,
	 	eps = 1e-2, k= 3, maxIter = 10000, **kwargs):
		'''em_fit performs a recursive fit of the the SAS data using the algorithm
		 taken from DOI: 10.1137/15M1024354.
		 	Args:
				-qRange (np.array): The extremes between which to fit the data.
					Defaults to [0,np.inf]
				-Rmin (float): The smallest radius to consider for the fit. The
					units are those used for q.	Defaults to 0.1
				-Rmax (float): The largest radius to use for the fit. If set to 0
					the radius will be calculated based on the smallest q-value
					available. Defaults to 0
				-numbElements (int): the number of elements in the redius vector.
					a larger number will result in a higher precision, btu also
					in longer calculation times and possible artifacts arising.
					Defaults to 100
				-eps (float): one of the two parameters defined in the paper to
					avoid overfitting the data. Defaults to 1e-2
				-k (int): second parameter defined in the paper to avoid
					overfitting the data. Defaults to 3
				-maxIter (int): maximum number of iterations to do before stopping
					the algorithm even if the convergence was not achieved.

		'''
		if kwargs.get('verbose', False):
			logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
		else:
			logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)


		convFact = 1
		if self.qUnits is not None:
			convFact = _UnitsConv[kwargs.get('units','nm')]/_UnitsConv[self.qUnits]
		logging.debug('orignial number of points (q-values): {}'.format(len(self.q)))
		#qRange = kwargs.get('qRange',np.array([0,np.inf]) )
		logging.debug("the given qRange is: {}, while the extremes of the q vector are: {}, {}".format(qRange,self.q[0],self.q[-1]))
		mask = np.logical_and(self.q>min(qRange), self.q<max(qRange))
		tempq = self.q[mask]
		tempI = self.Iabs[mask]
		logging.debug('After applying qRange the number of points is : '.format(len(tempq)))
		xk, Rvec, H = EM.expected_maximization_method(tempq, tempI, Rmin = Rmin,\
							Rmax = Rmax, eps = eps, k = k, numbElements =numbElements,\
							maxIter = maxIter, verbose = kwargs.get('verbose',False))

		xk = xk.flatten()
		Volume = np.array(4./3.*Rvec**3*np.pi)
		totVolume = np.sum(xk*Volume)
		numbDist = np.array(xk/Volume)
		volDist = np.array(Volume*xk/totVolume)

		self.fitResults['EM'] = {'Rvec': Rvec, 'H': H, 'xk': xk, \
		'units': kwargs.get('units','nm'), 'qRange': qRange, 'volDist': volDist,\
		'numbDist': numbDist}

	def lmfit_from_em(self,forced_model = None, qRange = np.array([0,np.inf]),\
	fromVolume = False, **kwargs):
		'''Uses the results from the EM fit as a starting point for the lmfit method.
		looks for at most 2 peaks and uses gaussian distributions. Implemetation for the Schultz
		distribution still have to be made
			Args:
				-forced_model (String): model imposed (single or double, Gauss or
					Schultz). If None the best model based on the number of peaks
					founf. Defaults to None
				-qRange (np.ndarray): Range over which the fit should be done.
					Defaults to [0,np.inf]
				-fromVolume (bool): If True the fitted distribution will be the
					volume distribution, if False the number distribution will
					be used. Default is False.
				kwargs (dict):
					order (int): parameter to find the maxima defined as:
						'How many points on each side to use for the comparison to consider'
						Defaults to 3
					max_mod(String): tells argrelmax whether the data extremes should be considered
						as maxima ('wrap') or not ('clip') defaults to 'wrap'
		'''
		if kwargs.get('verbose', False):
			logging.basicConfig(format='%(levelname)s:%(message)s', level=10)
		else:
			logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
		#qRange = kwargs.get('qRange',np.array([0,np.inf]) )

		if self.fitResults['EM'] is None:
			logging.info('SAS.lmfit_from_em(): Expected maximization has not been performed on the sample. It will be performed using the default arguments')
			self.em_fit(qRange = qRange)

		Rvec = self.fitResults['EM']['Rvec']
		if fromVolume:
			xk = self.fitResults['EM']['volDist']
		else:
			xk = self.fitResults['EM']['numbDist'].flatten()
		paramSugg, fitType, redchi = Fitting_Models.guess_from_dist(Rvec,xk,fitType = forced_model, verbose = kwargs.get('verbose', False), goodness = True)
		self.fit_data(fitType = fitType, qRange = qRange, paramSugg=paramSugg,verbose = kwargs.get('verbose',False),fit_kws = kwargs.get('fit_kws',{}))

		return 0

	def plot_fit(self, fitType, ax = None, fig = None, yshift = 0, plotDistribution = True,\
	latexCode = False,  **kwargs):
		"""plot_fit: used to plot the results of the fit. Creates 2 axis, on the left the
		data plus the theoretical scattering curve (in solid line in the fitted range and as
		dashed line outside of the range). On the right plots the particle distribution as
		 a histogram. If the two axis are provided it used the given axes.
			Args:
				fitType (string): the type of fit to plot
				ax (list of :obj: plt.axes): the two axis used to plot the curves
					and the distribution. If not provided they will be created.
					Defaults to None
				yshift (int): used to shift the data vertically, useful when
					plotting multiple fits on the same graph
				plotDistribution: sets whether the distribution whould be plotted
					along iwht the scattering curve.
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
			#print 'logging set to debug'
			logging.basicConfig(format='%(levelname)s:%(message)s', level = logging.DEBUG)
		else:
			#print 'logging set to info'
			logging.basicConfig(format='%(levelname)s:%(message)s', level = logging.INFO)

		qConvFact = 1
		IConvFact = 1
		figTitle = '{} of {}'.format(_AvlbSASFitDicInv[fitType], self.sampleName)

		if (kwargs.get('qUnits',None) in _AvlbUnits) and (self.qUnits in _AvlbUnits):
			qConvCoeff = _UnitsConv[self.qUnits]/_UnitsConv[kwargs.get('qUnits',None)]
		else:
			qUnits = self.qUnits

		if (kwargs.get('IUnits',None) in _AvlbUnits) and (self.IUnits in _AvlbUnits):
			IConvCoeff = _UnitsConv[self.IUnits]/_UnitsConv[kwargs.get('IUnits',None)]
		else:
			IUnits = self.IUnits


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

			if self.fitResults[fitType] is None:
				logging.error("No fit using {} has been performed yet!".format(fitType))
				return 0
			else:
				params = self.create_params(fitType)
				for p in params:
					params[p].set(value = self.fitResults[fitType][p])
					#singleGauss[p])
				RedqRange = self.fitResults[fitType]['qRange']

		if fitType == "EM":
			if self.fitResults['EM'] is None:
				logging.error("No expected maximization has been performed yet!")
				return 0
			else:
				RedqRange = self.fitResults['EM']['qRange']
				distModel = 'EM'

		#The axes are set-up, the first one to plot the scattering curve+fit, the second
		#to plot the distribution
		if plotDistribution and (ax is None or len(ax) != 2):
		#if ax is None or len(ax) != 2:
			fig = Figure(figsize=kwargs.get('figsize',(22,10)))
			gs = matplotlib.gridspec.GridSpec(1,3)
			ax = [fig.add_subplot(gs[0,0:2]), fig.add_subplot(gs[0,2])]

		elif (not plotDistribution) and (ax is None):
			fig = Figure(figsize=kwargs.get('figsize',(10,10)))
			#gs = matplotlib.gridspec.GridSpec(1,3)
			ax = [fig.add_subplot(111)]

		ax[0].tick_params(labelsize = kwargs.get('labelSize',18),size = kwargs.get('tickSize',6),
		 				width = kwargs.get('tickWidth',3) )


		#ax[0].ticklabel_format(style='sci',scilimits=(-2,2),axis='both')
		yFrmt = ScalarFormatter()
		yFrmt.set_powerlimits((-2,2))
		yFrmt.set_scientific(True)
		ax[1].yaxis.set_major_formatter(yFrmt)
		if plotDistribution:
			ax[1].set_ylabel('Number Distribution', fontsize = kwargs.get('lableSize',20))
			if qUnits in _AvlbUnits:
				ax[1].set_xlabel('Particle Size ({})'.format(_UnitsSymbols[qUnits]), fontsize = kwargs.get('lableSize',20))
			else:
				ax[1].set_xlabel('Particle Size', fontsize = kwargs.get('lableSize',20))
			ax[1].tick_params(labelsize = kwargs.get('labelSize',18),size = kwargs.get('tickSize',6),
			 				width = kwargs.get('tickWidth',3) )

		if qUnits in _AvlbUnits:
			ax[0].set_xlabel('Wavevector ({}$^{{-1}}$)'.format(_UnitsSymbols[qUnits]), fontsize = kwargs.get('lableSize',20))
		else:
			ax[0].set_xlabel('Wavevector', fontsize = kwargs.get('lableSize',20))
		if IUnits in _AvlbUnits:
			ax[0].set_ylabel('Intensity ({}$^{{-1}}$)'.format(_UnitsSymbols[qUnits]), fontsize = kwargs.get('lableSize',20))
		else:
			ax[0].set_ylabel('Intensity', fontsize = kwargs.get('lableSize',20))
		#If the error is available the experimental data is plotted with error
		#bars, if not as normal scatter points

		if self.Ierr is not None and kwargs.get('plot_error',False):
			ax[0].errorbar(self.q*qConvFact, self.Iabs*IConvFact*10**yshift, yerr=self.Ierr*IConvFact,\
							color = self.fit_plot_dict['color'], linestyle = self.fit_plot_dict['linestyle'],\
							marker = self.fit_plot_dict['marker'], ms = self.fit_plot_dict['ms'],\
							markeredgecolor = self.fit_plot_dict['mec'])
			#color = kwargs.get('expColor','k'), linestyle = kwargs.get('expLinestyle','None'),\
			#marker = kwargs.get('marker','o'), ms = kwargs.get('ms',2), \
			#markeredgecolor = kwargs.get('expColor','k'))
			ax[0].set_xscale("log", nonposx = "clip")
			ax[0].set_yscale("log", nonposx = "clip")
		else:
			ax[0].loglog(self.q*qConvFact, self.Iabs*IConvFact*10**yshift,\
						color = self.fit_plot_dict['color'], linestyle = self.fit_plot_dict['linestyle'],linewidth = self.fit_plot_dict['linewidth'],\
						marker = self.fit_plot_dict['marker'], ms = self.fit_plot_dict['ms'],\
						markeredgecolor = self.fit_plot_dict['mec'])

		#print 'THe reduced q range {}'.format(RedqRange)
		mask = np.logical_and(self.q>min(RedqRange), self.q<max(RedqRange))
		tempq = self.q[mask]

		if plotDistribution:
			#If the radii range in which to plot the distribution is not given it has to be calculated.
			#This is done using the _find_distribution_limits function.
			rLimits = kwargs.get('radRange',None)
			if rLimits is None:
				if fitType in _lmfitModels:
					rLimits = self._find_distribution_limits(fitType, distModel,params,\
					 limitRatio = kwargs.get('limitRatio',1000.))
				else:
					rLimits = self._find_distribution_limits(fitType, distModel)
			logging.debug('Calculated limits for radii: {}'.format(rLimits))
			radii = np.arange(min(rLimits), max(rLimits), (max(rLimits)-min(rLimits))/2000)

		if fitType in _lmfitModels:
			logging.debug('Plotting {}'.format(fitType))
			ax[0].loglog(self.q*qConvFact, 10**yshift*model.eval(params = params, q=self.q),\
							color = self.fit_plot_dict['fit_color'], linestyle = '--', linewidth = self.fit_plot_dict['fit_linewidth'],\
							marker = self.fit_plot_dict['fit_marker'], ms = self.fit_plot_dict['fit_ms'],\
							markeredgecolor = self.fit_plot_dict['fit_mec'])
			ax[0].loglog(tempq*qConvFact, 10**yshift*model.eval(params = params, q = tempq),\
							color = self.fit_plot_dict['fit_color'], linestyle = self.fit_plot_dict['fit_linestyle'], linewidth = self.fit_plot_dict['fit_linewidth'],\
							marker = self.fit_plot_dict['fit_marker'], ms = self.fit_plot_dict['fit_ms'],\
							markeredgecolor = self.fit_plot_dict['fit_mec'])
			if plotDistribution:
				ax[1].plot(radii/qConvFact, distModel.eval(params = params, x = radii),\
							color = self.fit_plot_dict['fit_color'], linestyle = self.fit_plot_dict['fit_linestyle'], linewidth = self.fit_plot_dict['fit_linewidth'],\
							marker = self.fit_plot_dict['fit_marker'], ms = self.fit_plot_dict['fit_ms'],\
							markeredgecolor = self.fit_plot_dict['fit_mec'])
		elif fitType == 'EM':
			#self._plot_em_result(ax, RedqRange, **kwargs)
			#print 'tempq shape', tempq.shape
			#expI = np.dot(self.EMFit['H'], self.EMFit['xk']).squeeze()
			#print 'H shape', expI.shape
			ax[0].loglog(tempq*qConvFact, 10**yshift*np.dot(self.fitResults['EM']['H'], self.fitResults['EM']['xk'])*IConvFact,\
						color = self.fit_plot_dict['fit_color'], linestyle = self.fit_plot_dict['fit_linestyle'], linewidth = self.fit_plot_dict['fit_linewidth'],\
						marker = self.fit_plot_dict['fit_marker'], ms = self.fit_plot_dict['fit_ms'],\
						markeredgecolor = self.fit_plot_dict['fit_mec'])
			#The H matrix has to be recreated to calculate the whole q range
			Vvec = self.fitResults['EM']['Rvec']**3*4/3*np.pi
			FF = EM.sphere_form_factor(self.q,self.fitResults['EM']['Rvec'])
			H = FF*np.matlib.repmat(Vvec,len(self.q),1)

			ax[0].loglog(self.q*qConvFact, 10**yshift*np.dot(H, self.fitResults['EM']['xk'])*IConvFact,\
						color = self.fit_plot_dict['fit_color'], linestyle = '--', linewidth = self.fit_plot_dict['fit_linewidth'],\
						marker = self.fit_plot_dict['fit_marker'], ms = self.fit_plot_dict['fit_ms'],\
						markeredgecolor = self.fit_plot_dict['fit_mec'])

			if plotDistribution:
				if kwargs.get('volumeDist', False):
					ax[1].plot(self.fitResults['EM']['Rvec']/qConvFact, self.fitResults['EM']['volDist'],\
							color = self.fit_plot_dict['fit_color'], linestyle = self.fit_plot_dict['fit_linestyle'], linewidth = self.fit_plot_dict['fit_linewidth'],\
							marker = self.fit_plot_dict['fit_marker'], ms = self.fit_plot_dict['fit_ms'],\
							markeredgecolor = self.fit_plot_dict['fit_mec'])
					ax[1].set_ylabel('Volume Distribution', fontsize = kwargs.get('lableSize',20))
				else:
					ax[1].plot(self.fitResults['EM']['Rvec']/qConvFact, self.fitResults['EM']['numbDist'],\
							color = self.fit_plot_dict['fit_color'], linestyle = self.fit_plot_dict['fit_linestyle'], linewidth = self.fit_plot_dict['fit_linewidth'],\
							marker = self.fit_plot_dict['fit_marker'], ms = self.fit_plot_dict['fit_ms'],\
							markeredgecolor = self.fit_plot_dict['fit_mec'])

		else:
			logging.info("Plotting of other fitting methods has not been implemented yet")
			return 0
		#plt.show()
		fig.suptitle(figTitle, fontsize = kwargs.get('lableSize',20)+2 )
		if latexCode:
			lCode = fit_result_latex(fitType = fitType)
			return (ax,fig, lCode)
		else:
			return (ax,fig)


	def fit_result_latex(self,fitType = ['Sing_Gauss'], fileOutput = None, preamble = False, units = None, close = True):
		'''fit_result_latex: exports the fits for the samples in latex format.
		The list can be then written to a file specified in fileOutput. The user
		can set which fits to export, whether the latex preamble should be included
		anche which units to use for the size of the particles.
			Args:
				-fitType(list): The fits which should be exported in the table.
					Defaults to ['Sing_Gauss'].
				-fileOutput(string): If this is parameter is specified the list
					is then written to the specified file. Defaults to None
				-preamble(bool): Specifies whether the latex preamble containing
					the packages to use should be inserted. If fileOutput is
					specified it is automatically set to true. Defaults to False.
				-units(string): Specifies which units should be used for the
					output. If set to None the qUnits will be used.
					Defaults to None
				-close(bool): Specifies whether \end{document} should be appended
					at the end. If fileOutput is specified this is automatically
					set to True
			-Returns:
				List containing the latex code needed ot create the tables.

		'''
		convCoeff = 1
		if units is not None:
			if (units in _AvlbUnits) and (self.qUnits in _AvlbUnits):
				convCoeff = _UnitsConv[units]/_UnitsConv[self.qUnits]
			else:
				units = self.qUnits
		else:
			units = self.qUnits

		if fileOutput is not None:
			preamble = True
			close = True

		if preamble:
			lCode = [r'\documentclass[10pt,a4paper]{article}',\
			r'\usepackage[a4paper, left = 1.5cm,right = 1.5cm,heightrounded]{geometry}',\
			r'\usepackage[utf8]{inputenc}',\
			r'\usepackage{amsmath}',\
			r'\usepackage{amsfonts}',\
			r'\usepackage{amssymb}',\
			r'\usepackage{subcaption}',\
			r'\usepackage[section]{placeins}',\
			r'\usepackage{graphicx}',\
			r'\usepackage[output-decimal-marker={.},alsoload=synchem]{siunitx}',\
			r'\DeclareSIUnit\Molar{\textsc{m}}',\
			r'\usepackage{booktabs}',\
			r'\renewcommand{\arraystretch}{1.8}',\
			r'\begin{document}']
		else:
			lCode = []
		lCode.append(r'\begin{tabular}{l p{2cm} p{2cm} p{2cm} p{2cm} p{2cm} p{2cm}}')
		header = r'\textbf{{Sample Name}}& $\mathbf{{\bar{{R}}}}$ ({0})& $\mathbf{{\sigma}}$ ({0})& \textbf{{Polydisp. (\%)}}& \textbf{{Fit Type}}& $\mathbf{{Red\chi}}$& \textbf{{Comments}}\\'.format(units)
		lCode.append(header)
		lCode.append(r'\toprule')
		if ('Sing_Gauss' in fitType) and (self.fitResults['Sing_Gauss'] is not None):
			R = self.fitResults['Sing_Gauss']['R_av']*convCoeff
			Rstderr = self.fitResults['Sing_Gauss']['R_av_stderr']*convCoeff
			sigma = self.fitResults['Sing_Gauss']['sigma']
			sigmaStderr = self.fitResults['Sing_Gauss']['sigma_stderr']
			redChi = self.fitResults['Sing_Gauss']['redchi']
			tableLine = r'{}& {:.3g}$\pm${:.3g}& {:.3g}$\pm${:.3g}& {:.3g}& {}&{:.3e}&\\'.format(self.sampleName, R, Rstderr, sigma, sigmaStderr, float(sigma)/R, 'Single Gaussian Dist.', redChi)
			lCode.append(tableLine)
		if ('Double_Gauss' in fitType) and (self.fitResults['Double_Gauss'] is not None):
			R1 = self.fitResults['Double_Gauss']['R1_av']*convCoeff
			R1stderr = self.fitResults['Double_Gauss']['R1_av_stderr']*convCoeff
			sigma1 = self.fitResults['Double_Gauss']['sigma1']
			sigma1Stderr = self.fitResults['Double_Gauss']['sigma1_stderr']
			R2 = self.fitResults['Double_Gauss']['R2_av']*convCoeff
			R2stderr = self.fitResults['Double_Gauss']['R2_av_stderr']*convCoeff
			sigma2 = self.fitResults['Double_Gauss']['sigma2']
			sigma2Stderr = self.fitResults['Double_Gauss']['sigma2_stderr']
			ratio = self.fitResults['Double_Gauss']['ratio']
			redChi = self.fitResults['Double_Gauss']['redchi']
			tableLine = r'{}& {:.3g}$\pm${:.3g}\newline{:.3g}$\pm${:.3g}& {:.3g}$\pm${:.3g}\newline{:.3g}$\pm${:.3g}& {:.3g}\newline{:.3g}& {}& {:.3e}&{:.3g}\%\newline{:.3g}\%\\'.format(self.sampleName, R1, R1stderr, R2, R2stderr, sigma1, sigma1Stderr, sigma2, sigma2Stderr, float(sigma1)/R1, float(sigma2)/R2, 'Bimodal Gaussian Dist.',redChi, ratio*100,(1-ratio)*100)
			lCode.append(tableLine)
		if ('Sing_Schultz' in fitType) and (self.fitResults['Sing_Schultz'] is not None):
			R = self.fitResults['Sing_Schultz']['R_av']*convCoeff
			Rstderr = self.fitResults['Sing_Schultz']['R_av_stderr']*convCoeff
			sigma = R/np.sqrt(1.+self.fitResults['Sing_Schultz']['Z'])
			sigmaStderr = sigma*self.fitResults['Sing_Schultz']['Z_stderr']/self.fitResults['Sing_Schultz']['Z']
			logging.debug('The stderr on the sigma for a Schultz distribution is calculated using the \% stderr on Z.')
			redChi = self.fitResults['Sing_Schultz']['redchi']
			tableLine = r'{}& {:.3g}$\pm${:.3g}& {:.3g}$\pm${:.3g}& {:.3g}& {}& {:.3e}&\\'.format(self.sampleName, R, Rstderr, sigma, sigmaStderr, float(sigma)/R, 'Single Schultz Dist.', redChi)
			lCode.append(tableLine)
		if ('Double_Schultz' in fitType) and (self.fitResults['Double_Schultz'] is not None):
			R1 = self.fitResults['Double_Schultz']['R1_av']*convCoeff
			R1stderr = self.fitResults['Double_Schultz']['R1_av_stderr']*convCoeff
			sigma1 = R/np.sqrt(1.+self.fitResults['Double_Schultz']['Z1'])
			sigma1Stderr = sigma*self.fitResults['Double_Schultz']['Z1_stderr']/self.fitResults['Double_Schultz']['Z1']
			R2 = self.fitResults['Double_Schultz']['R2_av']*convCoeff
			R2stderr = self.fitResults['Double_Schultz']['R2_av_stderr']*convCoeff
			sigma2 = R/np.sqrt(1.+self.fitResults['Double_Schultz']['Z2'])
			sigma2Stderr = sigma*self.fitResults['Double_Schultz']['Z2_stderr']/self.fitResults['Double_Schultz']['Z2']
			ratio = self.fitResults['Double_Schultz']['ratio']
			logging.debug('The stderr on the sigma for a Schultz distribution is calculated using the \% stderr on Z.')
			redChi = self.fitResults['Double_Schultz']['redchi']
			tableLine = r'{}& {:.3g}$\pm${:.3g}\newline{:.3g}$\pm${:.3g}& {:.3g}$\pm${:.3g}\newline{:.3g}$\pm${:.3g}& {:.3g}\newline{:.3g}& {}& {:.3e}&{:.3g}\%\newline{:.3g}\%\\'.format(self.sampleName, R1, R1stderr, R2, R2stderr, sigma1, sigma1Stderr, sigma2, sigma2Stderr, float(sigma1)/R1, float(sigma2)/R2, 'Bimodal Schultz Dist.', redChi,ratio*100,(1-ratio)*100)
			lCode.append(tableLine)
		lCode.append(r'\end{tabular}')
		lCode.append(r'\newline')
		if close:
			lCode.append(r'\end{document}')
		if fileOutput is not None:
			with open(fileOutput,'w') as fn:
				for line in lCode:
					fn.write("%s\n" % line)
		return lCode


	def _find_distribution_limits(self, fitType, model, params = None, limitRatio = 1000., verbose = False):
		"""_find_distribution_limits: determins the limits of the distribution plot
		if they are not given. This is done by taking the peak (for bimodal dist.
		the left peak is used for the left limit while the right peak is used for
		the right limit), and calculating the value of the distribution.
		from 0 to the peak. The left boundary will be taken as the value at which the
		distribution reaches 1/n the value at the peak (n=100?), or 0. The right boundary
		will be taken by calculating the distribution between the rightmost peak and a
		position simmetrically opposit the left peak. The limit will be selected as the
		value at which the distribution reaches 1/n the value at the peak (n=1000) or the
		symmetric position.
			Args:
				-fitType (string): the type of fit to plot
				-model (lmfit.Model): the lmfit model used to calculate the distribution
				-limitRatio (int): The curve will be calculated up to where the value is
					equal to Imax/limitRatio. Default: 1000
				-params (lmfit.Parameters): the lmfit parameters object which contains
					the parameters of the distribution
			Returns:
				-Numpy.array with the two values in increasing order
		"""
		if verbose:
			#print 'logging set to debug'
			logging.basicConfig(format='%(levelname)s:%(message)s', level = logging.DEBUG)
		else:
			#print 'logging set to info'
			logging.basicConfig(format='%(levelname)s:%(message)s', level = logging.INFO)

		if fitType in _lmfitModels:
			if params is None:
				logging.error("To calculate the radii limits for an lmfit model you need to provide the prameters")
				return 0

			if "Sing" in fitType:
				R_av = params['R_av'].value
			elif "Double" in fitType:
				R_av = min(params['R1_av'].value, params['R2_av'].value)
				logging.debug('The smalles peak is at: {}'.format(R_av))
			xs = np.arange(0,R_av,R_av/100.)
			ys = model.eval(params = params, x = xs)

			lowerThan = xs[ys<ys[-1]/limitRatio]
			if lowerThan.size == 0:
				lLimit = 0
			else:
				lLimit = lowerThan[-1]
			if "Double" in fitType:
				R_av = max(params['R1_av'].value, params['R2_av'].value)

			step = 1
			cont = True
			while cont:
				xs = np.arange(R_av,2**step*R_av, R_av/(100.*step))
				ys = model.eval(params = params, x = xs)
				lowerThan = xs[ys<ys[0]/limitRatio]

				if lowerThan.size == 0:
					step += 1
					if step > 7:
						step = 7
						cont = False
						rLimit = 2**step*R_av
				else:
					rLimit = lowerThan[0]
					cont = False
		elif fitType == 'EM':
			#Should be changed so that the maxima are found, the first and second
			lLimit = self.fitResults['EM']['Rvec'][0]
			rLimit = self.fitResults['EM']['Rvec'][-1]


		return np.array([lLimit,rLimit])

	def fit_guinier_regime(self, error_lim = 0.01, min_q = 3,skip_initial = 0, plot = False, verbose = False):
		"""guinier_regime: fit the data in the guinier regime using the linear relationship
		   ln(I(q)) = -q^2*Rg/3
		   Rg = -m*3
			Args:
				error_lim (long): acceptable relative error between the linear fit and the
					data. the value has to be between 0 and 1, 0 excluded. Defaults to 0.01
				min_q (int): minimum number of points needed to validate the Guinier regime
					fit. Has to be greater than 3. Defaults to 3
				skip_initial (int): number of initial points to ignore
				plot_guinier (bool): whether to plot the guinier fit or not. Defaults to
					False
		"""
		if verbose:
			logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
		else:
			logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
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
		logging.debug('Starting to fit the guinier regime from the {} point'.format(curr_q))
		#guinierRegime = True
		advance = True
		while advance:
			params = linMod.guess(I_log[:curr_q], x = q_2[:curr_q])
			result = linMod.fit(I_log[:curr_q],params = params, x = q_2[:curr_q])
			logging.debug('Completed a fit with {} points.\n'.format(curr_q))
			if result.params['slope'].value>0:
				logging.error('Guinier approzimation cannot be calculated if the slope has a positive value. Try skipping more initial points.')
				return
			Rg = np.sqrt(-result.params['slope'].value*3)
			error = np.mean(result.residual/I_log[:curr_q])
			logging.debug('Completed a fit with {} points.\nRg:{}\nerror:{}'.format(curr_q,Rg,error))
			guinierRegime = Rg*self.q[(skip_initial+curr_q)]<1
			logging.debug('The Rg*q is: {}'.format(Rg*self.q[(skip_initial+curr_q)]))
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
		if plot:
			fig = Figure(figsize = (6,6))
			ax = fig.add_subplot(111)
			ax.plot(q_2[:(2*curr_q)], I_log[:(2*curr_q)], 'k', marker = 'o', linestyle = 'None')
			ax.plot(q_2[:curr_q], result.eval(x = q_2[:curr_q]), 'r', linewidth = 2)
			ax.plot(q_2[:(2*curr_q)], result.eval(x = q_2[:(2*curr_q)]), 'r', linestyle = '--',linewidth = 2)

		self.guinier_regime = {'fit_Result': result,'Numb_qs': curr_q,'q_lim': self.q[curr_q],'Guinier_valid': guinierRegime, 'Rg': Rg}
		return {'best_fit': result,'last_q_pos': skip_initial+curr_q, \
				'last_q': self.q[(skip_initial+curr_q)],'within_Guinier': guinierRegime,\
				'initial_points_skipped': skip_initial, 'Rg': Rg}

	def find_porod_invariant(self, drho, units = 'm', qRange = [0,np.inf], precision = 0.1, postExtrap = None,\
						preExtrap = None, plot = False):
		"""porod_invariant: calculate the porod invariant for the given data. The porod
		invariant is calculated as the integral of I*q^2. The data is integrated, and
		then extrapolated to infinity by giving a constant slope to the signal of the
		type I(q) = q^(slope). The complete equation is:
		Q = pi^2*phi*(1-phi)*(drho)^212
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

		if plot:
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

		if fitType not in _AvlbSASFit:
			logging.error('{} is not a valid fit type for SAS data'.format(fitType))
		if self.fitResults[fitType] is None:
			return None
		if param is None:
			return self.fitResults[fitType]
		else:
			if isinstance(param, list):
				OrderedDict([(i, i * i) for i in range(5)])
				return OrderedDict([(p, self.fitResults[fitType].get(p, None)) for p in param])
			elif isinstance(param, str):
				return {param: self.fitResults[fitType].get(param, None)}
			else:
				logging.error('The format in which the parameter was passed is not recognizable: \n {} \n Please provide either a string or a list of strings'.format(param))

	#def fit_avlb(self, fitType):
		#if fitType in _AvlbSASFit:
