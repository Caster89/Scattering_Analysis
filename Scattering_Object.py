from config import _AvlbUnits, _UnitsConv, _AvlbSASFit, _AvlbWASFit, _lmfitModels, _lmfitModelFunctions, _lmfitDistFunctions 
import numpy as np
#import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Qt5Agg")
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import pickle
import os
import os.path
import logging
from mpl_toolkits.mplot3d import Axes3D
try :
	import lmfit
except ImportError:
	print "Could not import the lmfit package. Some fittings will be disabled."
	lmfitAvlb = False
else:
	lmfitAvlb = True


class ScatteringObject(object):
	"""Base Class used as parent for all the other scattering\n
	classes 
	"""
	
	def __init__(self, fname = None, **kwargs):
		"""Initialize all the necessary variables
		the child classes will add extra variables
		where they are needed.
		fname: The complete name of the file (directory+filename)
		kwarg:
			- qcol: the column of the q-values
			- Icol: the column of the intensities
			- skip_header: the number of header lines to skip
			- delimiter: the delimiter used in the csv
			- units: the units used in the file
		"""
		self.q = None
		self.Iraw = None
		self.Ibckg = None
		self.coeffBckg = None
		self.coeffAbs = None
		self.Inorm = None
		self.Iabs = None
		self.isAbs = False
		self.Ierr = None
		self.IbckgErr = None
		self.qUnits = kwargs.get('qUnits', None)
		self.IUnits = kwargs.get('IUnits', None)
		self.setup = {'DetDist' : kwargs.get('DetDist', None),\
					  'Temp' : kwargs.get('Temp', None), \
					  'Dt' : kwargs.get('Dt', None),\
					  'Thickness' : kwargs.get('Thickness', None),\
					  'Wavelength': kwargs.get('Wavelength',None)}
		self.plot_dict = {'color' : 'k', 'linestyle' : 'None', 'linewidth' : 2, 'marker' : 'o', 'ms' : 4,\
						  'mec' : 'None'}
		self.fit_plot_dict = {'color' : 'k', 'linestyle' : 'None', 'linewidth' : 2, 'marker' : 'o', 'ms' : 4,\
						  'mec' : 'None', 'fit_color' : 'r', 'fit_linestyle' : '-', 'fit_linewidth' : 2, \
						  'fit_marker' : 'None', 'fit_ms' : 2, 'fit_mec' : 'None'}
		self.sampleName = kwargs.get('SampleName', None)
		if fname is not None:
			#Create list of keywords required in create_from_file and use them to 
			#populate a dictionary to pass as argument
			createKeywords = ['qcol','Icol', 'skip_header', 'delimiter', 'bckgFname', 'units', 'Errcol']
			passArguments = {kw : kwargs[kw] for kw in createKeywords if kw in kwargs}
			self.create_from_file(fname, **passArguments)
		
	def create_from_file(self, fname, **kwargs):
		"""Reads the data from a file into the object.
		
		Args:
			fname (string): the complete path to the file which contains the scattering data
			qcol (int): the column which contains the q-wavevector. Defaults to 0
			Icol (int): the column which contains the intensity vector. Defaults to 1
			Errcol (int): the column which contains the estimated error. If
				-1 no error is imported. Defaults to -1
			skip_header (int): the number of lines to skip at the beginning of the
				file. defaults to 2
			delimiter 'string': the character used as a separator between columns.
				Defaults to ';'
			qunits (string): the units in which the q-vector is given. They have to be
				contained withing the _AvlbUnits list in the config.py file. If None
					Defaults to 'nm'
			Iunits (string): the units in which the intensity is given. They have to be
				contained within the _AvlbUnits list in the config.py file. If None
					Defaults to 'nm'
			bkcgFname (string): the filename containing the background intensity. If None
				no background is imported. Defaults to None
			ReadHeader (dict): dictionary containing the information ot read the header in
				order to retrieve the information to store in Setup dictionary:
				
				Elements:
					linePos (int): indicates which line to use for the header. linePos has to
						be <skip_header. Defaults to 1
					delimiter (string): the delimiter used in the header. Defaults to the delimiter
						used in the file
					field (dict): dictionary containing the name of the field and its position
		"""
		#The values used to import the data are either read, or if not passed
		#are set to their default value
		
		qcol = kwargs.get('qcol', 0)
		Icol = kwargs.get('Icol', 1)
		Errcol = kwargs.get('Errcol', -1)
		skip_header = kwargs.get('skip_header',2)
		delimiter = kwargs.get('delimiter', ';')
		bckgFname = kwargs.get('bckgFname', None)
		qunits = kwargs.get('qunits', None)
		Iunits = kwargs.get('Iunits', None)
		
		temp = np.genfromtxt(fname,skip_header=skip_header, delimiter = delimiter)
		self.q = temp[:, qcol]
		self.Iraw = temp[:, Icol]
		#print 'Importing raw I:', self.Iraw
		self.Inorm = temp[:,Icol]
		self.Iabs = temp[:,Icol]
		if Errcol != -1:
			self.Ierr = temp[:,Errcol]
			
		self.qUnits = 'nm'
		self.IUnits = 'm'
		if qunits is not None:
			if qunits in _AvlbUnits:
				self.qUnits = qunits
			else:
				print "{} is not a recognized unit. The availabel units are:\n".format(qunits)
				for au in _AvlbUnits:
					print "\t - {}".format(au)
				print r"q-Units are set to nm^(-1)"
				

			
		if Iunits is not None:
			if Iunits in _AvlbUnits:
				self.IUnits = Iunits
			else:
				print "{} is not a recognized unit. The availabel units are:\n".format(Iunits)
				for au in _AvlbUnits:
					print "\t - {}".format(au)
				print r"q-Units are set to m^(-1)"
		if bckgFname is not None:
			temp = np.genfromtxt(bckgFname,skip_header=skip_header, delimiter = delimiter)
			#temp = temp[:,Icol]
			if len(temp) == len(self.q):
				self.Ibckg = temp[:,Icol]
				#print 'Importing bck: ', self.Ibckg
				if Errcol != -1:
					self.IbckgErr = temp[:,Errcol]
			else:
				print r"Size mismatch: the Bckg vector must have the same size as the I \n \
				 vector, the background was not imported."

		if "ReadHeader" in kwargs:
			print "Header interpretation not implemented yet"

	def create_from_data(self, q, I, Ibckg = None, qunits = 'nm', Iunits = 'm'):
		"""Sets data by passing both the q values and the intesity
			Args:
				q (list): the q vector.
				I (list): the intensity vector
				Ibckg (list): the background intensity vector. Defaults to None
				qunits (string): the units of the q wavevector. Has to be in _AvlbUnits
					in the config.py file. Defaults to 'nm'
				Iunits (string): the units of the intensity vector. Has to be in _AvlbUnits
					in the config.py file. Defaults to 'm'
		"""
		if len(q) == len(I):
			self.q = np.array(q)
			self.Iraw = np.array(I)
		else:
			print "Size mismatch: the q vector must have the same size as the I vector"
			return 0
		if bckgI is not None:
			if len(Ibckg) == len(q):
				self.Ibckg = np.array(Ibackg)
			else:
				print r"Size mismatch: the Bckg vector must have the same size as the I \n \
				 vector, the background was not imported."
		if qunits is not None:
			if qunits in _AvlbUnits:
				self.qUnits = qunits
			else:
				print r"{} is not a recognized unit. The available units are:\n".format(qunits)
				for unt in AvlbUnits:
					print r"\t - ", unt
		if Iunits is not None:
			if Iunits in _AvlbUnits:
				self.IUnits = Iunits
			else:
				print r"{} is not a recognized unit. The available units are:\n".format(Iunits)
				for unt in AvlbUnits:
					print r"\t - ", unt
				
	def find_bckg_scale(self, qRange = [0,np.inf], applyNorm = False, plot=False, ax = None):
		"""Compares the raw intensity to the background within the 
		qRange specified. It finds the best coefficient by which to multiply the 
		background so as to minimize the difference between the two. If the lmfit module
		is present it does it using the lmfit functions, if not by using a simple average.
		if apply is set to True then the coefficient will be stored and the Inorm will
		be updated. If plot is selected then it will plot the raw data and the corrected
		background.
			Args:
				qRange (list): contains the q range limits used to normalize the background.
					Defaults to [0, np.inf].
				applyNorm (bool): tells whether the normalization should be applied. This sets
					the CoeffBckg to the value found, and the Inorm vector to Iraw-CoeffBckg*Ibckg.
					Defaults to False
				plot (bool): tells whether the beckground normaliztion should be plotted. Defaults
					to False
		"""
		tempq = self.q[np.logical_and(self.q>min(qRange), self.q<max(qRange))]
		tempI = self.Iraw[np.logical_and(self.q>min(qRange), self.q<max(qRange))]
		tempb = self.Ibckg[np.logical_and(self.q>min(qRange), self.q<max(qRange))]
		
		if lmfitAvlb:
			p = lmfit.Parameters()
			p.add('a',value = 1, min = 0.0)
			fit = lmfit.minimize(lambda pars,int,bckg: int-pars['a'].value*bckg ,p ,args = (tempI,tempb))
			coeff = fit.params['a'].value
		else:
			coeff = np.mean(tempI/tempb)
		
		print "The normalizing coefficient for the background is: {}".format(coeff)
		if applyNorm:
			self.coeffBckg = coeff
			self.Inorm = self.Iraw-self.coeffBckg*self.Ibckg
			if self.IbckgErr is not None:
				self.IbckgErr = self.IbckgErr*coeff
		
		if plot:
			#print 'plotting'
			if ax is None:
				fig = Figure(figsize = (6,6))
				ax = fig.add_subplot(111)
			ax.loglog(self.q, self.Iraw, 'k', marker = 'o', linestyle = 'None')
			ax.loglog(self.q, self.Ibckg, 'b')
			ax.loglog(self.q, coeff*self.Ibckg, 'r')
			ax.legend(['Exp. Data', 'Exp. Bckg', 'Norm. Backg'], frameon= False)
			#plt.show()
			
	def set_abs_scale(self, absCoeff, thickness = None, bckgCoeff = None):
		"""set_abs_scale: applies the corrections necessary to set the data in absolute
		scale. Requires that the coeffBckg and thickness be set.
			Args:
				absCoeff (int): the integer by which the signal has to be multiplied in
					order to obtain the data in absolute scale
				thickness (int): the thickness of the capillary. It can either be set or
					changed here. Defaults to None
				bckgCoeff (int): the coefficient by which to multiply the background
					for a correct subtraction. It can be set or changed from here.
					Defgaults to None
		"""
		if self.Ibckg is None:
			print "Cannot calculate the signal at absolute scale without a background"
			return 0
		if thickness is not None:
			print "Setting the thickness to {} {}. Make sure the units are correct".format(thickness,self.IUnits)
			self.setup['Thickness'] = thickness
		if bckgCoeff is not None:
			self.coeffBckg = bckgCoeff
			self.Inorm = self.Iraw-self.coeffBckg*self.Ibckg
			self.IbckgErr = self.IbckgErr*self.coeffBckg
			
		if self.setup['Thickness'] is None:
			print "Cannot calculate the signal at absolute scale without the sample's thickness"
			return 0
		if self.coeffBckg is None:
			print "Cannot calculate the signal at absolute scale without the background"
			return 0
		self.coeffAbs = absCoeff
		self.Iabs = (self.Iraw-self.coeffBckg*self.Ibckg)*self.setup['Thickness']*self.coeffAbs
		if self.Ierr is not None:
			self.Ierr =self.Ierr*self.setup['Thickness']*self.coeffAbs
		if self.IbckgErr is not None:
			self.IbckgErr =self.IbckgErr*self.setup['Thickness']*self.coeffAbs
	
	def plot_data(self, **kwargs):
		print "Method not implemented in parent classe"
	
	def fit_data(self, ):
		print "Method not implemented in parent classe"
	
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
	def load_from_file(directory, fileName):
		if os.path.isfile(os.path.join(directory,fileName)):
			with open(os.path.join(directory, fileName),'rb') as f:
				return pickle.load(f)
		else:
			print os.path.join(directory, fileName), ' is not a file'
			return None