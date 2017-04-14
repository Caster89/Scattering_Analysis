from config import _AvlbUnits, _UnitsConv, _AvlbSASFit, _AvlbWASFit, _lmfitModels, _lmfitModelFunctions, _lmfitDistFunctions 
from Scattering_Object import ScatteringObject
from SAS import SmallAngleScattering
from WAS import WideAngleScattering
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
	
class SWAS(object):
	
	def __init__(self, Sample_Name = 'Sample', SA_Fname = None, WA_Fname = None, SA_dict = None, WA_dict = None, **kwarg):
		"""Initialized the object. This is basically a wrapper function to initialize
		either the SAS object or the WAS object or both.
			Args:
				SA_Fname (string): the name of the file containing the SAS data. If None
					the SAS is not
					
				WA_Fname (string): the name of the file containing the WAS
		"""
		self.SAS = None
		self.WAS = None
		self.sampleName = Sample_Name
		self.avlbCurves = None
		#Copying the dictionary is necessary if not when setting the sample name
		#the dictionary would be modified (python always "passes by reference").
		#This caouses problems with the SWAS Sequence class
		

		
		if SA_Fname is not None:
			tempSA_dict = SA_dict.copy()
			if SA_dict.get('SampleName', None) is None:
				tempSA_dict['SampleName'] = 'SA_'+ self.sampleName
			self.import_SA_from_file(SA_Fname, **tempSA_dict)
			self.avlbCurves = 'SAS'
		if WA_Fname is not None:
			tempWA_dict = WA_dict.copy()
			if WA_dict.get('SampleName', None) is None:
				tempWA_dict['SampleName'] = 'WA_'+ self.sampleName
			self.import_WA_from_file(WA_Fname, **tempWA_dict)
			if self.avlbCurves == 'SAS':
				self.avlbCurves = 'Both'
			else:
				self.avlbCurves = WAS
	def __getitem__(self,arg):
		"""__getitem__: funtion to overide the square braket operator
		here it is used to obtain the various vectors from both the
		small angles and the wide angles.
		"""
		if arg == 'qs':
			return self.SAS.q
		elif arg == 'Is':
			return self.SAS.IAbs
		elif arg == 'IsBckg':
			return self.SAS.Ibckg
		elif arg == 'qw':
			return self.WAS.q
		elif arg == 'Iw':
			return self.WAS.IAbs
		elif arg == 'IwBckg':
			return self.WAS.Ibckg
	
	def import_SA_from_file(self, fname, **kwargs):
		self.SAS = SmallAngleScattering(fname, **kwargs)
	
	def import_WA_from_file(self, fname, **kwargs):
		self.WAS = WideAngleScattering(fname, **kwargs)
		
	def plot_data(self, qsRange=[0,np.inf],qwRange = [0, np.inf], ysShift=1, ywShift = 1, axs = None, fig = None,\
				  xUnits = 'nm', yUnits ='m', **kwargs):
		"""plot_data: plots both the SAS and WAS data in the same image. Takes the same parameters as
		the individual plotting functions, and passes them over.
			Args:
				qsRange (list): the q range for the small angle plot. Defaults to [0, np.inf]
				qwRange (list): the q range for the wide angle plot. Defaults to [0,np.inf]
				ysShift (int): the vertical shift for the small angle curve. Defaults to 0
				ywShift (int): the vertical shift for the wide angle curve. Defaults to 0
				ax (list): list with the 2 axis for plotting the SAS and WAS curves. The
					list must have 2 axes. If None it create the axes. Defaults to None
				fig (matplotlig.figure): the figure onto which the data is plotted. If None
					creates the figure using the 'figSize' keyword. Defaults to None
				xUnits (string): units for the x-axis. Defaults to 'nm'
				yUnits (string): unist for the y-axis. Defaults to 'm'
			Returns:
				the axes list and the figure
		"""
		#print 
		#axs = kwargs.pop('axs', None)
		#print 'the axs array has length :', len(axs)
		#print 'and they are: ', axs
		
		
		#fig = kwargs.pop('fig', None)
		if fig is None:
			#print 'had to create the figure'
			fig = Figure(figsize = kwargs.pop('figSize', (12,6)))
		if axs is None or len(axs) != 2:
			#print 'had to create the axis'
			axs = []
			axs.append(fig.add_subplot(121))
			axs.append(fig.add_subplot(122))
		
		self.SAS.plot_data(qRange = qsRange, yShift = ysShift, ax = axs[0], **kwargs)
		self.WAS.plot_data(qRange = qwRange, yShift = ywShift, ax = axs[1], **kwargs)
		
		return [axs, fig]
	
	def fit_guinier_regime(self, error_lim = 0.01, min_q = 3,skip_initial = 0, plot_guinier = False):
		return self.SAS.fit_guinier_regime( error_lim, min_q,skip_initial, plot_guinier)
	
	def find_porod_invariant(self,drho, units = 'm', qRange = [0,np.inf], precision = 0.1, postExtrap = None,\
						preExtrap = None, plot_porod = False):
		return self.SAS.find_porod_invariant(drho, units, qRange, precision, postExtrap,\
						preExtrap, plot_porod)
	
	def fit_data(self, fitType, **kwargs):
		if fitType in _AvlbSASFit:
			self.fit_SAS(fitType, **kwargs)
		elif fitTYpe in _AvlbWASFit:
			self.fit_WAS(fitType, **kwargs)
		else:
			print 'Fit type {} is not recognized'.format(fitType)
	
	def fit_SAS(self, fitType = "Sing_Gauss",**kwargs):
		self.SAS.fit_data(fitType,**kwargs)
	
	def fit_WAS(self, fitType = 'PVM', **kwargs):
		self.WAS.fit_peaks(fitType, **kwargs)
	
	def plot_fit(self, fitType, ax = None, **kwargs):
		if fitType in _AvlbSASFit:
			plot_SAS_fit(fitType, ax, **kwargs)
		if fitType in _AvlbWASFit:
			plot_WAS_fit(fitType, ax, **kwargs)
	
	def plot_SAS_fit(self, fitType,ax = None, **kwargs):
		return self.SAS.plot_fit(fitType,ax, **kwargs)
	
	def plot_WAS_fit(self, fitType,ax = None, **kwargs):
		return self.WAS.plot_fit(fitType,ax, **kwargs)
	
	def get_Temp(self):
		if self.avlbCurves is 'Both':
			Temps = self.SAS.setup['Temp']
			Tempw = self.WAS.setup['Temp']
			if Temps == Tempw:
				return Temps
			else:
				print 'The two curves have a different temperature, cannot return a \
				unique value'
				return None
		elif self.avlbCurves is 'SAS':
			Temps = self.SAS.setup['Temp']
			return Temps
		elif self.avlbCurves is 'WAS':
			Tempw = self.WAS.setup['Temp']
			return Tempw
		
	
	def set_Temp(self, Temp):
		if self.SAS is not None:
			self.SAS.setup['Temp'] = Temp
		if self.WAS is not None:
			self.WAS.setup['Temp'] = Temp
	
	def get_Dt(self):
		if self.avlbCurves is 'Both':
			Dts = self.SAS.setup['Dt']
			Dtw = self.WAS.setup['Dt']
			if Dts == Dtw:
				return Dts
			else:
				print 'The two curves have a different times, cannot return a \
				unique value'
				return None
		elif self.avlbCurves is 'SAS':
			return self.SAS.setup['Dt']
		elif self.avlbCurves is 'WAS':
			return 	AS.setup['Dt']
	
	def set_Dt(self, Dt):
		if self.SAS is not None:
			self.SAS.setup['Dt'] = Dt
		if self.WAS is not None:
			self.WAS.setup['Dt'] = Dt
	
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
			with open(os.path.join(directory, fileName),'rb') as f:
				return pickle.load(f)
		else:
			print os.path.join(directory, fileName), ' is not a file'
			return None

