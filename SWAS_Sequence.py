from config import _AvlbUnits, _UnitsConv, _AvlbSASFit, _AvlbWASFit, _lmfitModels, _lmfitModelFunctions, _lmfitDistFunctions 
from Scattering_Object import ScatteringObject
from SAS import SmallAngleScattering
from WAS import WideAngleScattering
from SWAS import SWAS
import os
import os.path
import pickle
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
	
class SWAS_Sequence(object):
	"""Wrapper class used to collect the data from a time resolved scattering experiment.
	It consists in a list of SWAS objects. It enables to do fits on all samples consecutively
	as well as plot the results as a function of time
	"""
	
	def __init__(self,SAS_fname_list = None, WAS_fname_list = None, SA_dict = None, WA_dict = None, **kwargs):
		#self.SAS_objects = None
		#self.WAS_objects = None
		self.SWAS_objects = None
		self.size = 0
		self.time = None
		self.Temp = None
		self.avlbCurves = None
		self.sampleName = kwargs.get('SampleName','Sample')
		
		if SAS_fname_list is not None or WAS_fname_list is not None:
			self.create_from_file_list(SAS_fname_list, WAS_fname_list, SA_dict, WA_dict)

	
	def create_from_file_list(self,SAS_fname_list = None, WAS_fname_list = None, SA_dict = None, WA_dict = None):
		"""create_from_file_list: function to create the SAS/WAS objects and place them in the
		appropriate list. It also initializes the time and temperature lists.
		"""
		if SAS_fname_list is not None and WAS_fname_list is not None:
			if len(SAS_fname_list) != len(WAS_fname_list):
				print 'In order to correctly create the object, the number of SAS and WAS curves\
				must be the same, or one of the two must not be given. SAS: {}, WAS: {}'.format(len(SAS_fname_list),\
																								len(WAS_fname_list))
				return 0
			self.avlbCurves = 'Both'
			self.SWAS_objects = []
			i = 0
			for sasFile, wasFile in zip(SAS_fname_list,WAS_fname_list):

				self.SWAS_objects.append(SWAS(Sample_Name = self.sampleName+'_{}'.format(i),\
										 SA_Fname = sasFile,\
										 WA_Fname = wasFile,\
										 SA_dict = SA_dict,\
										 WA_dict = WA_dict))
				i += 1
		elif SAS_fname_list is not None and WAS_fname_list is None:
			self.avlbCurves = 'SAS'
			self.SWAS_objects = []
			i = 0
			for sasFile in SAS_fname_list:
				self.SWAS_objects.append(SWAS(Sample_Name = self.sampleName+'_{}'.format(i),\
										 SA_Fname = sasFile,\
										 WA_Fname = None,\
										 SA_dict = SA_dict,\
										 WA_dict = None))
				i += 1
		elif SAS_fname_list is None and WAS_fname_list is not None:
			self.avlbCurves = 'WAS'
			self.SWAS_objects = []
			i = 0
			for sasFile in SAS_fname_list:
				self.SWAS_objects.append(Sample_Name = self.sampleName+'_{}'.format(i),\
										 SA_Fname = None,\
										 WA_Fname = wasFile,\
										 SA_dict = None,\
										 WA_dict = WA_dict)
				i += 1
		self.time = np.arange(len(self.SWAS_objects))
		self.Temp = np.array([23]*len(self.SWAS_objects))
		for i,obj in enumerate(self.SWAS_objects):
			if obj.get_Dt() is not None:
				self.time[i] = obj.get_Dt[i]
			if obj.get_Temp() is not None:
				self.Temp[i] = obj.get_Temp[i]
		self.size = len(self.SWAS_objects)
	
	def find_bckg_scale_SAS(self, qRange = [0,np.inf], applyNorm = False, plot=False, ax = None):
		for i in range(self.size):
			self[i].SAS.find_bckg_scale(qRange, applyNorm, plot, ax)
	
	def find_bckg_scale_WAS(self, qRange = [0,np.inf], applyNorm = False, plot=False, ax = None):
		for i in range(self.size):
			self[i].WAS.find_bckg_scale(qRange, applyNorm, plot, ax)
	
	def set_abs_scale_SAS(self, absCoeff, thickness = None, bckgCoeff = None):
		for i in range(self.size):
			print i
			self[i].SAS.set_abs_scale(absCoeff,thickness,bckgCoeff)
			
	def set_abs_scale_WAS(self, absCoeff, thickness = None, bckgCoeff = None):
		for i in range(self.size):
			self[i].WAS.set_abs_scale(absCoeff,thickness,bckgCoeff)
		
	def set_temperature(self, Temp):
		"""set_Temperature: sets the temperature vector, as well as copying the values in the idividual
		scattering objects:
			Args:
				Temp (int or list): ether an int (if the temperature is constant), or a list.
					If a list is provided it must have the same lengths a the number of objects
		"""
		if isinstance(Temp, int):
			self.Temp = Temp
			for i in xrange(len(self.SWAS_objects)):
				self.SWAS_objects[i].set_Temp(Temp)
		elif isinstance(Temp,list):
			if len(Temp) != len(self.SWAS_objects):
				print 'If  list of temperatures is given it must have the same size as the number\
				of curves. Numb. Curves: {}, Numb. Temp: {}'.format(len(self.SWAS_objects),len(Temp))
				return 0
			
			self.Temp = np.array(Temp)
			for i in xrange(len(self.SWAS_objects)):
				self.SWAS_objects[i].set_Temp(Temp[i])
	
	def set_time(self, time):
		"""set_time: sets the time vector, as well as copying the values in the individual
		scattering objects:
			Args:
				Temp (int or list): ether an int (if the temperature is constant), or a list.
					If a list is provided it must have the same lengths a the number of objects
		"""
		
		if len(time) != len(self.SWAS_objects):
			print 'If  list of time is given it must have the same size as the number\
			of curves. Numb. Curves: {}, Numb. Temp: {}'.format(len(self.SWAS_objects),len(time))
			return 0
			
		self.time = np.array(time)
		for i in xrange(len(self.SWAS_objects)):
				self.SWAS_objects[i].set_Dt(time[i])
	
	def plot_SAS_Sequence(self, qRange=[0,np.inf], frequency=10, yShift = 10, ax = None , colormap = matplotlib.cm.jet, **kwargs):
		if ax is None:
			fig = Figure(figsize = figSize)
			ax = fig.add_subplot(111)
		for i,idx in enumerate(range(0,self.size, frequency)):
			self[idx].SAS.plot_data(qRange, yShift = 10**(yShift*i), ax = ax, color = colormap(float(idx)/float(self.size)))
			
	def plot_WAS_Sequence(self, qRange=[0,np.inf], frequency=10, yShift = 10, ax = None , colormap = matplotlib.cm.jet, **kwargs):
		if ax is None:
			fig = Figure(figsize = figSize)
			ax = fig.add_subplot(111)
		for i,idx in enumerate(range(0,self.size, frequency)):
			self[idx].WAS.plot_data(qRange, yShift = yShift*i, ax = ax, color = colormap(float(idx)/float(self.size)))
		
	def fit_data(self, fitType, fitList = -1, memory = False, **kwargs):
		"""fit_data is used to fit the data in the SWAS_Sequence file
			Args:
				fitType (String): a string containg the fit type to use for the fit. It is
					passed to the individual fitting functions.
				fitList [list of int]: a list containing the position of objects to fit. If
					-1 the all the data objects are fitted. Defaults to -1
				memory (bool): tells the function whether it should use the previous fit result
					as the starting point for the current fit. Defaults to False
				**kwargs (dict): the various kwargs which can be passed to the fitting function
				"""
		if kwargs.get('verbose', False):
			logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
		else:
			logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
			
			
		if fitList == -1:
			fitList = range(self.size)
		if fitType in _AvlbSASFit:
			if memory:
				self.SWAS_objects[fitList[0]].fit_SAS(fitType, **kwargs)
				params = self.SWAS_objects[fitList[0]].SAS.get_fit_params(fitType)
				for i in fitList[1:]:
					logging.debug('Fitting curve {}'.format(i))
					try:
						self.SWAS_objects[i].fit_SAS(fitType, paramSugg = params)
						params = self.SWAS_objects[i].get_fit_params(fitType)
					except:
						logging.error('At object {}, lmfit failed to do the fit'.format(i))
						
			else:
				for i in fitList:
					logging.debug('Fitting curve {}'.format(i))
					try:
						self.SWAS_objects[i].fit_SAS(fitType, **kwargs)
					except:
						print 'At object {}, lmfit failed to do the fit'.format(i)
						
		elif fitType in _AvlbWASFit:
			if memory:
				self.SWAS_objects[fitList[0]].fit_WAS(fitType, **kwargs)
				params = self.SWAS_objects[fitList[0]].get_fit_params(fitType)
				for i in fitList[1:]:
					logging.debug('Fitting curve {}'.format(i))
					self.SWAS_objects[i].fit_WAS(fitType, **params)
					params = self.SWAS_objects[i].get_fit_params(fitType)
			else:
				for i in fitList:
					logging.debug('Fitting curve {}'.format(i))
					try :
						self[i].fit_WAS(**kwargs)
					except ZeroDivisionError:
						print 'At object {}, lmfit failed to do the fit'.format(i)
					
		else:
			print 'Fit type {} is not recognized'.format(fitType)
			
	def plot_param(self, fitType, param, ax = None, **kwargs):
		'''plot_param: function used to plot the evolution of a parameter
		over the sequence.
			Args:
			fitType (String): the fit for which the parameter should be
				plotted. Can be any one of the fits available
			param (String): the parameter which should be plotted. In case
				of WAS peak fitting this can either be the parameter
				of a single peak or the parameter of all the peaks
			ax (list: mpl.axes): a list of axis on which the parameters
				should be plotted.
			**kwargs (dict): the various parameters for the plotting of the data
		'''
		
		if fitType in _AvlbSASFit:
			if fitType in _lmfitModels:
				plotData = [obj.SAS.get_fit_params(fitType, param = param) for obj in self.SWAS_objects]
				numbPlots = len(plotData[0])
				if ax is None or len(ax)<numbPlots:
					ax = []
					fig = Figure(figsize = kwargs.get('figSize',(6,(6*numbPlots))))
					for i in xrange(numbPlots):
						ax.append(fig.add_subplot(numbPlots,1,(i+1)))
				self._plot_lmfit_param(plotData,ax)
			elif fitType == 'EM':
				plotData = [obj.SAS.EMFit for obj in self.SWAS_objects]
				self._plot_em_param(ax, plotData, **kwargs)
				

		if fitType in _AvlbWASFit:
			if fitType == 'PVM':
				plotData = [obj.WAS.get_fit_params(fitType, param = param) for obj in self.SWAS_objects]
				numbPlots = len(plotData[0])
				if ax is None or len(ax)<numbPlots:
					ax = []
					fig = Figure(figsize = kwargs.get('figSize',(6,(6*numbPlots))))
					for i in xrange(numbPlots):
						ax.append(fig.add_subplot(numbPlots,1,(i+1)))
				self._plot_PVM_param(plotData,ax)
	
	def _plot_PVM_param(self, paramList, ax, **kwargs):
		'''_plot_PVM_param: specific function used to plot the evolution
		of the parameters from a Pseudo-Voigt fit.
		'''
		numPlots = len(paramList[0])
		keys = paramList[0].keys()
		for i,k in enumerate(keys):
			tempPlot = [params.get(k,None) if params is not None else None for params in paramList]
			ax[i].scatter(self.time,tempPlot, c = self.time, cmap = kwargs.get('cmap', matplotlib.pyplot.cm.jet), edgecolor = 'face')
	
	def _plot_lmfit_param(self, paramList,ax, **kwargs):
		'''_plot_lmfit_param: specific function used to plot the evolution
		of the parameters from a lmfit fit
		'''
		numPlots = len(paramList[0])
		keys = paramList[0].keys()
		for i,k in enumerate(keys):
			tempPlot = [params.get(k,None) if params is not None else None for params in paramList]
			ax[i].scatter(self.time,tempPlot, c = self.time, cmap = kwargs.get('cmap', matplotlib.pyplot.cm.jet), edgecolor = 'face')
	
	def _plot_em_param(self, ax, paramList, **kwargs):

		if kwargs.get('verbose', False):
			print 'verbose is true'
			logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
		else:
			print ' verbose is false'
			logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
			

		#logging.info('Plotting the evolution of the EM curve is still not implemented')
		#if not isinstance(ax, Axes3D):
			#fig = kwargs.get('fig',mpl.figure())
			#ax = Axes3D(fig)
		Dist = np.empty((self.size, paramList[0]['Rvec'].size))
		print 'There are {} data sets'.format(self.size)
		print 'The radius vector is divided in {} points'.format(paramList[0]['Rvec'].size)
		for i in xrange(self.size):
			
			Dist[i] = paramList[i]['xk'].reshape((-1))
		ax.imshow(Dist.T)
		#ax.set_ylim(paramList[0]['Rvec'][0],paramList[0]['Rvec'][-1])
		#ax.set_xlim(self.time[0],self.time[-1])
	def save_to_file(self, directory = None, fileName = None):
		"""save_to_file: saves the class to a file.
			Args:
				directory (string): the directory in which to save
					the file. If None then the current working directory
					is used. Defaults to None
				fileName (string): the name of the file in which to save
					the class' instance. If None the sample name is used,
					If there is no sample name then saves under a number.
					Defaults to None
		"""
		if directory is None:
			directory = os.getcwd()
		if fileName is None:
			if self.sampleName is not None:
				fileName = self.sampleName
			else:
				n = 0
				exists = True
				while exists:
					if os.path.isfile(os.path.join(directory,str(n)+'.p')):
						n += 1
					else:
						fileName = str(n) + '.p'
						exists = False
		
		with open(os.path.join(directory, fileName),'wb') as f:
			pickle.dump(self,f)
		
	@staticmethod
	def load_from_file(directory, filename):
		if os.path.isfile(os.path.join(directory,filename)):
			with open(os.path.join(directory, filename),'rb') as f:
				return pickle.load(f)
		else:
			print os.path.join(directory, fileName), ' is not a file'
			return None
		
	def __getitem__(self,n):
		return self.SWAS_objects[n]