from config import _AvlbUnits, _UnitsConv, _UnitsSymbols, _AvlbSASFit, _AvlbWASFit, _lmfitModels, _lmfitModelFunctions, _lmfitDistFunctions
import numpy as np
#import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Qt5Agg")
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import pickle
import collections
from collections import OrderedDict
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
		fname: The complete name of the file (directory+filename). This can also
			be a list of filenames, in case multiple images were taken for the
			sample.
		kwarg:
			- qcol: the column of the q-values
			- Icol: the column of the intensities
			- skip_header: the number of header lines to skip
			- delimiter: the delimiter used in the csv
			- qUnits: the units for the wavevector used in the file
			- IUnits: the units used for the intensity in the file
			- DetDistance: the distance between the sample and detector
			- Temp: the temperature of the sample
			- Dt: the time elapsed since the beginning of the measurement
			- Wavelength: the wavelength of the incident beam
			- readHeader (dict): dictionary containing the information to read
				the header in order to retrieve the information to store in
				setup dictionary:
				Elements:
					- linePos (int): indicates which line to use for the header.
						linePos has to be < skip_header. Defaults to 0
					- delimiter (string): the delimiter used in the header.
						Defaults to the delimiter used in the file
					- field (dict): dictionary containing the name of the field
						and its position
		"""

		"""
		The standard values for all scattering objects are:
			-q: The values of the wavevector
			-Iframes: The raw data from the radial average of each individual frame
			-Iraw: the raw data from the radial average of the scattering pattern
			-Ibckg: the background intensity such as empty capillary or solvent
			-Idark: The dark of the ccd camera, where applicable
			-coeffBckg: The coefficient used to correct the intensity of the background
			-coeffAbs: The coefficient used to bring the data's intensity to absolute units
			-Inorm: The raw data with the background subtracted
			-Iabs: The normalized date corrected for absolute values
			-Ierr: The error on the scattering intensity
			-qerr: The error on the q-values
			-IbckgErr: The error on the background intensity
			-qUnits: The units used for the wavevector
			-IUnits: The units used for the scattered intensity
			-setup: The data on the instumental setup:
				-DetDistance: The sample detector distance
				-Temp: The temperature during the measurement
				-Dt: The time at which the measurement was done (time resolved exp)
				-Thickness: The capillary's thickness
				-Wavelength: The wavelength of the x-ray
			-plot_dict: All parameters for plotting the scattered data, must be
						keywords used in matplotlib
			-fit_plot_dict: All parameterd for plotting the fitted curve, must
						be keywords used in matplotlib
			-sampleName: The name of the sample

		"""
		self.q = None
		self.Iframes = None
		self.Iraw = None
		self.IbckgFrames = None
		self.Ibckg = None
		self.Idark = None
		self._coeffBckg = 1
		self._coeffAbs = 1
		self.Inorm = None
		self.Iabs = None
		self.isAbs = False
		self.Ierr = None
		self.qerr = None
		self.IbckgErr = None
		self._qUnits = kwargs.get('qUnits', None)
		self._IUnits = kwargs.get('IUnits', None)
		self._setup = {'DetDist' : kwargs.get('DetDist', None),\
					  'Temp' : kwargs.get('Temp', None), \
					  'Dt' : kwargs.get('Dt', None),\
					  'Thickness' : kwargs.get('Thickness', 1),\
					  'Wavelength': kwargs.get('Wavelength',None)}
		self._plot_dict = {'color' : 'k', 'linestyle' : 'None', 'linewidth' : 2, 'marker' : 'o', 'ms' : 4,\
						  'mec' : 'None'}
		self._fit_plot_dict = {'color' : 'k', 'linestyle' : 'None', 'linewidth' : 2, 'marker' : 'o', 'ms' : 4,\
						  'mec' : 'None', 'fit_color' : 'r', 'fit_linestyle' : '-', 'fit_linewidth' : 2, \
						  'fit_marker' : 'None', 'fit_ms' : 2, 'fit_mec' : 'None'}
		print self.fit_plot_dict
		print self._fit_plot_dict
		self._sampleName = kwargs.get('SampleName', None)
		if fname is not None:
			#Create list of keywords required in create_from_file and use them to
			#populate a dictionary to pass as argument
			createKeywords = ['qcol','Icol', 'skip_header', 'delimiter', 'bckgFname',\
			'qUnits', 'IUnits', 'Errcol', 'readHeader', 'verbose']
			passArguments = {kw : kwargs[kw] for kw in createKeywords if kw in kwargs}
			self.create_from_file(fname, **passArguments)

	def create_from_file(self, fname, **kwargs):
		"""Reads the data from a file into the object.
		Args:
			fname (string/stirng list): the complete path to the file which
				contains the scattering data. The data can be either a string
				or a list of strings.
			qcol (int): the column which contains the q-wavevector. Defaults to 0
			Icol (int): the column which contains the intensity vector. Defaults to 1
			Errcol (int): the column which contains the estimated error. If
				-1 no error is imported. Defaults to -1
			skip_header (int): the number of lines to skip at the beginning of the
				file. defaults to 2
			delimiter 'string': the character used as a separator between columns.
				Defaults to ';'
			qUnits (string): the units in which the q-vector is given. They have to be
				contained withing the _AvlbUnits list in the config.py file. If None
					Defaults to 'nm'
			IUnits (string): the units in which the intensity is given. They have to be
				contained within the _AvlbUnits list in the config.py file. If None
					Defaults to 'nm'
			bkcgFname (string): the filename containing the background intensity. If None
				no background is imported. Defaults to None
			readHeader (dict): dictionary containing the information to read the header in
				order to retrieve the information to store in Setup dictionary:

				Elements:
					linePos (int): indicates which line to use for the header. linePos has to
						be < skip_header. Defaults to 0
					delimiter (string): the delimiter used in the header. Defaults to the delimiter
						used in the file
					field (dict): dictionary containing the name of the field and its position
		"""
		if kwargs.get('verbose', False):
			logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
		else:
			logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
		#The values used to import the data are either read, or if not passed
		#are set to their default value
		qcol = kwargs.get('qcol', 0)
		Icol = kwargs.get('Icol', 1)
		Errcol = kwargs.get('Errcol', -1)
		skip_header = kwargs.get('skip_header',2)
		delimiter = kwargs.get('delimiter', None)
		bckgFname = kwargs.get('bckgFname', None)
		readHeader = kwargs.get('readHeader',{})
		if ('linePos' not in readHeader) or ('delimiter' not in readHeader):
			readHeader = None
		#qunits = kwargs.get('qunits', None)
		#Iunits = kwargs.get('Iunits', None)
		#If fname is a string then only one file is read
		if isinstance(fname, str):
			temp = np.genfromtxt(fname,skip_header=skip_header, delimiter = delimiter)
			self.q = temp[:, qcol]
			self.Iraw = temp[:, Icol]
			self.IFrames = temp[:,Icol]
			if Errcol != -1:
				self.Ierr = temp[:,Errcol]

		elif isinstance(fname, collections.Iterable):
			'''
			If fname is an iterable then it is supposed that each element is a
			file and is loaded in Iframes.
			'''
			temp = np.genfromtxt(fname[0],skip_header=skip_header, delimiter = delimiter)
			self.q = temp[:, qcol]
			self.Iframes = np.empty( (len(self.q), len(fname)) )
			for cc, fn in enumerate(fname):
				self.Iframes[:,cc] = np.genfromtxt(fn,skip_header=skip_header, delimiter = delimiter, usecols = Icol)
			self.Iraw = self.Iframes.mean(axis = 1)
			##TODO: Add the error in case of multi frames. In this case it
			#should be either and "average" of the errors per frame or a std
			#deviation between the different frames
		'''
		Inorm and Iabs are initiated as equal to the raw data in order for the
		class to be usable in case the user does not intend to do any
		normalisation on the data.
		'''
		self.Inorm = temp[:,Icol]
		self.Iabs = temp[:,Icol]

		#This block should have been substituted by setting qUnits and IUnits
		#as properties
		'''
		if kwargs.get('qUnits', None) is not None:
			if kwargs.get('qUnits') in _AvlbUnits:
				self.qUnits = kwargs.get('qUnits')
			else:
				print "{} is not a recognized unit. The availabel units are:\n".format(kwargs.get('qUnits'))
				for au in _AvlbUnits:
					print "\t - {}".format(au)
				print r"q-Units are set to None"

		if kwargs.get('IUnits', None) is not None:
			if kwargs.get('IUnits') in _AvlbUnits:
				self.IUnits = kwargs.get('IUnits')
			else:
				#self.IUnits = Iunits
			#else:
				print "{} is not a recognized unit. The availabel units are:\n".format(kwargs.get('IUnits'))
				for au in _AvlbUnits:
					print "\t - {}".format(au)
				print r"I-Units are set to None"
		'''
		#Import the background if a file is provided
		if bckgFname is not None:
			if isinstance(bckgFname, str):
				temp = np.genfromtxt(bckgFname,skip_header=skip_header, delimiter = delimiter)
				if len(temp[:,0]) == len(self.q):
					self.Ibckg = temp[:,Icol]
					#print 'Importing bck: ', self.Ibckg
					if Errcol != -1:
						self.IbckgErr = temp[:,Errcol]
				else:
					logging.error("Size mismatch: the Bckg vector must have the same size as the I vector, the background was not imported.")

			elif isinstance(bckgFname, collections.Iterable):
				temp = np.genfromtxt(bckgFname[0],skip_header=skip_header, delimiter = delimiter)
				if len(temp[:,0]) == len(self.q):
					self.IbckgFrames = np.empty( (len(self.q), len(bckgFname)) )
					for cc, fn in enumerate(bckgFname):
						self.IbckgFrames[:,cc] = np.genfromtxt(fn,skip_header=skip_header, delimiter = delimiter, usecols = Icol)
					self.Ibckg = self.IbckgFrames.mean(axis = 1)
				else:
					logging.error("Size mismatch: the Bckg vector must have the same size as the I vector, the background was not imported.")
					self.Ibckg = np.zeros(self.Iraw.shape)
					self.coeffBckg = 1

		else:
			self.Ibckg = np.zeros(self.Iraw.shape)
			self.coeffBckg = 1

		if readHeader is not None:
			if isinstance(fname, str):
				temp = np.genfromtxt(fname, use_rows = readHeader.pop('linePos'), delimiter = readHeader.pop('delimiter'))
			elif isinstance(fname, collections.Iterable):
				temp = np.genfromtxt(fname[0], use_rows = readHeader.pop('linePos'), delimiter = readHeader.pop('delimiter'))
			for kword in readHeader:
				self._setup[kword] = readHeader[kword]

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
		else:
			self.Ibckg = np.zeros(self.Iraw.shape)
			self.coeffBckg = 1
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

	def find_bckg_scale(self, qRange = [0,np.inf], applyNorm = False, plot=False,\
	 					ax = None, verbose = False):
		"""Compares the raw intensity to the background within the
		qRange specified. It finds the best coefficient by which to multiply the
		background so as to minimize the difference between the two. If the lmfit module
		is present it does it using the lmfit functions, if not by using a simple average.
		If apply is set to True then the coefficient will be stored and the Inorm will
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
		if verbose is False:
			logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
		else:
			logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

		mask = np.logical_and(self.q>min(qRange), self.q<max(qRange))
		tempq = self.q[mask]
		tempI = self.Iraw[mask]
		tempb = self.Ibckg[mask]

		if lmfitAvlb:
			p = lmfit.Parameters()
			p.add('a',value = 1)
			fit = lmfit.minimize(lambda pars,intensity,bckg: intensity-pars['a'].value*bckg ,p ,args = (tempI,tempb), nan_policy = 'omit')
			coeff = fit.params['a'].value
		else:
			coeff = np.mean(tempI/tempb)

		print "The normalizing coefficient for the background is: {}".format(coeff)

		if applyNorm:
			logging.info('The value of the background coeff was automatically updated')
			self.coeffBckg = coeff
			self.Inorm = self.Iraw-self.coeffBckg*self.Ibckg
			if self.IbckgErr is not None:
				self.IbckgErr = self.IbckgErr*coeff

		if plot:
			#print 'plotting'
			if ax is None:
				fig = Figure(figsize = (6,6))
				ax = fig.add_subplot(111)
			ax.loglog(self.q, self.Iraw, **self.plot_dict)
			ax.loglog(self.q, self.Ibckg, 'b')
			ax.loglog(self.q, coeff*self.Ibckg, 'r')
			ax.legend(['Exp. Data', 'Exp. Bckg', 'Norm. Backg'], frameon= False)
			#plt.show()

	def set_abs_scale(self, absCoeff=1, thickness = None, bckgCoeff = None):
		"""set_abs_scale: applies the corrections necessary to set the data in absolute
		scale. Requires that the coeffBckg and thickness be set.
			Args:
				absCoeff (int): the integer by which the signal has to be multiplied in
					order to obtain the data in absolute scale
				thickness (int): the thickness of the capillary. It can either be set or
					changed here. Defaults to None
				bckgCoeff (int): the coefficient by which to multiply the background
					for a correct subtraction. It can be set or changed from here.
					Defaults to None
		"""
		#if self.Ibckg is None:
			#print "Cannot calculate the signal at absolute scale without a background"
			#return 0
		if thickness is not None:
			print "Setting the thickness to {} {}. Make sure the units are correct".format(thickness,self.IUnits)
			self.setup['Thickness'] = thickness
		if bckgCoeff is not None:
			self.coeffBckg = bckgCoeff
			#self.Inorm = self.Iraw-self.coeffBckg*self.Ibckg
			#if self.IbckgErr is not None:
				#self.IbckgErr = self.IbckgErr*self.coeffBckg

		if self.setup['Thickness'] is None:
			print "Cannot calculate the signal at absolute scale without the sample's thickness"
			return 0
		#if self.coeffBckg is None:
			#print "Cannot calculate the signal at absolute scale without the background"
			#return 0
		self.coeffAbs = absCoeff
		#self.Iabs = (self.Iraw-self.coeffBckg*self.Ibckg)*self.setup['Thickness']*self.coeffAbs
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

	'''
	DEFINITION OF ALL THE PROPERTIES OF THE CLASS
	'''

	###UNITS USED FOR THE INTENSITY###
	@property
	def IUnits(self):
		"""Stores the units in which the intensity is expressed"""
		return self._IUnits

	@IUnits.setter
	def IUnits(self, value):
		if value in _AvlbUnits:
			self._IUnits = value
		else:
			print '{} is not an available unit. Please select between the following units:'.format(value)
			for un in _AvlbUnits:
				print '-{}: {}'.format(un, _UnitsSymbols[un])

	###UNITS USED FOR THE WAVEVECTOR###
	@property
	def qUnits(self):
		"""Stores the units in which the wavevector is expressed"""
		return self._qUnits

	@qUnits.setter
	def qUnits(self, value):
		if value in _AvlbUnits:
			self._qUnits = value
		else:
			print '{} is not an available unit. Please select between the following units:'.format(value)
			for un in _AvlbUnits:
				print '-{}: {}'.format(un, _UnitsSymbols[un])

	###NAME USED FOR THE SAMPLE###
	@property
	def sampleName(self):
		"""Stores the name of the sample"""
		return self._sampleName

	@sampleName.setter
	def sampleName(self, value):
		self._sampleName = value

	###DETAILS REGARDING THE EXPERIMENTAL SETUP###
	@property
	def setup(self):
		"""Stores the details of the experimental setup (eg. Detector distance, wavelength, etc.)"""
		return self._setup

	@setup.setter
	def setup(self, value):
		#Should check whether this is a dictionary (probably done by update)
		self._setup.update(value)

	###DETAILS FOR PLOTTING THE DATA###
	@property
	def plot_dict(self):
		"""Stores the options for plotting the data using matplotlib"""
		return self._plot_dict

	@plot_dict.setter
	def plot_dict(self, value):
		self._plot_dict.update(value)

	###DETAILS FOR PLOTTING THE FIT###
	@property
	def fit_plot_dict(self):
		"""Stores the options for plotting the fit using matplotlib"""
		return self._fit_plot_dict

	@fit_plot_dict.setter
	def fit_plot_dict(self, value):
		self._fit_plot_dict.update(value)

	###COEFFICIENT FOR THE NORMALISATION OF THE BACKGROUND###
	@property
	def coeffBckg(self):
		"""Stores coefficient used to subtract the background"""
		return self._coeffBckg

	@coeffBckg.setter
	def coeffBckg(self, value):
		self._coeffBckg = value
		self.Inorm = self.Iraw - self.coeffBckg*self.Ibckg
		self.Iabs = self.coeffAbs*self.Inorm
		if self.IbckgErr is not None:
			self.IbckgErr = self.IbckgErr*self.coeffBckg

	###COEFFICIENT FOR THE ABSOLUTE VALUE OF THE SIGNAL###
	@property
	def coeffAbs(self):
		"""Stores coefficient used to obtain the intensity in absolute value"""
		return self._coeffAbs

	@coeffAbs.setter
	def coeffAbs(self, value):
		self._coeffAbs = value
		self.Iabs = (self.Iraw-self.coeffBckg*self.Ibckg)*self.setup['Thickness']*self.coeffAbs
		if self.Ierr is not None:
			self.Ierr =self.Ierr*self.setup['Thickness']*self.coeffAbs
		if self.IbckgErr is not None:
			self.IbckgErr =self.IbckgErr*self.setup['Thickness']*self.coeffAbs

	###THICKNESS OF THE SAMPLE###
	@property
	def thickness(self):
		"""Stores the thickness ofthe sample"""
		return self._setup['Thickness']

	@thickness.setter
	def thickness(self, value):
		self._setup['Thickness'] = value
		self.Iabs = (self.Iraw-self._coeffBckg*self.Ibckg)*self._setup['Thickness']*self._coeffAbs
		if self.Ierr is not None:
			self.Ierr =self.Ierr*self.setup['Thickness']*self.coeffAbs
		if self.IbckgErr is not None:
			self.IbckgErr =self.IbckgErr*self.setup['Thickness']*self.coeffAbs
