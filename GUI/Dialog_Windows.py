import numpy as np
from PyQt5 import QtCore, QtGui, QtWidgets, Qt
from ..SWAS_Sequence import SWAS_Sequence
from ..SAS import SmallAngleScattering
from ..WAS import WideAngleScattering
from ..SWAS import SWAS
from ..SWAS_Sequence import SWAS_Sequence
from ..config import _AvlbUnits, _UnitsConv, _AvlbSASFit,_AvlbWASFit,_AvlbSASFitDic,_AvlbWASFitDic, _lmfitModels, _lmfitModelFunctions, _lmfitDistFunctions 
import PlottingWidgets as PW

SAS_Fitting = {'Single Gaussian Distribution' : 'Sing_Gauss',\
                   'Double Gaussian Distribution' : 'Double_Gauss',\
                   'Single Schultz Distribution' : 'Sing_Schultz',\
                   'Double Schultz Distribution' : 'Double_Schultz',\
                   'Expected Minimization' : 'EM'}
SAS_Fitting_Names = ['Single Gaussian Distribution',\
                   'Double Gaussian Distribution',\
                   'Single Schultz Distribution',\
                   'Double Schultz Distribution',\
                   'Expected Minimization']
class Fitting_Dialog(QtWidgets.QDialog):
    '''Class used to create the fitting dialog needed to fit the data
    '''
   
    def __init__(self, parent, scatteringObj, angleRange = 'SAS', **kwargs):
        super(Fitting_Dialog, self).__init__(parent, **kwargs)
        self.scatteringObj = scatteringObj
        self.angleRange = angleRange
        #Temporary fix
        self.scatteringObj = self.scatteringObj['Object']
        self.layout = QtWidgets.QVBoxLayout()
        
        self.fitType = QtWidgets.QComboBox()
        if angleRange == 'SAS':
            for f in SAS_Fitting_Names:
                self.fitType.addItem(f, SAS_Fitting[f])  
        else:
            print 'WAS fitting not implemented yet'
        self.fitType.currentIndexChanged.connect(self.fit_type_changed)
        
        self.topPanel = QtWidgets.QStackedWidget()
        self.bottomPanel = QtWidgets.QGroupBox()
        self.buttonPanel = QtWidgets.QGroupBox()
        self.buttonLayout = QtWidgets.QHBoxLayout()
        
        self.createTopPanel()
        self.createBottomPanel()
        
        self.fitButton = QtWidgets.QPushButton('Fit')
        self.cancelButton = QtWidgets.QPushButton('Cancel')
        self.cancelButton.clicked.connect(self.reject)
        self.fitButton.clicked.connect(self.accept)
        self.buttonLayout.addStretch(1)
        self.buttonLayout.addWidget(self.cancelButton)
        self.buttonLayout.addWidget(self.fitButton)
        self.buttonLayout.addStretch(1)
        self.buttonPanel.setLayout(self.buttonLayout)
        
        
        self.layout.addWidget(self.fitType)
        self.layout.addWidget(self.topPanel)
        self.layout.addWidget(self.bottomPanel)
        self.layout.addWidget(self.buttonPanel)
        
        self.setLayout(self.layout)
        
        self.min_limit = 0
        self.min_fit_limit = 0
        self.min_line = None
        self.max_limit = np.inf
        self.max_fit_limit = np.inf
        self.max_line = None
        
        #self.setGeometry(100,200,300,600)
    #def setResult(self,value):
        #if value == 0:
            #super(Fitting_Dialog,self).setResult(value)
        #else
    def result(self):
        if super(Fitting_Dialog,self).result() == 0:
            return super(Fitting_Dialog,self).result() == 0
        else:
            fit_params = self.topPanel.currentWidget().return_params()
            print fit_params
            print float(self.min_ledit.text())
            print float(self.max_ledit.text())
            fit_params['qRange'] = np.array([float(self.min_ledit.text()),float(self.max_ledit.text())])
            fit_type = str(self.fitType.currentData())
            #SAS_Fitting[SAS_Fitting_Names[self.fitType.currentIndex()]]
            return [fit_type, fit_params]
    def fit_type_changed(self):
        i = self.fitType.currentIndex()
        self.topPanel.setCurrentIndex(i)
       
    def createBottomPanel(self):
        self.bottomLayout = QtWidgets.QVBoxLayout(self)
        self.ctrlLayout = QtWidgets.QHBoxLayout(self)
        
        if isinstance(self.scatteringObj,SWAS_Sequence):
            plotData = {'Object' : self.scatteringObj, 'plot_type' : self.angleRange}
            if self.angleRange == 'SAS':
                plotData['SAS'] = range(self.scatteringObj.size)
            else:
                plotData['WAS'] = range(self.scatteringObj.size)
            self.plotWidget = PW.MultiplePlot(self.bottomPanel,plotData = plotData)
        elif isinstance(self.scatteringObj, SWAS):
            if self.angleRange == 'SAS':
                self.plotWidget = PW.SinglePlot(self.bottomPanel,plotData = self.scatteringObj.SAS)
            else:
                self.plotWidget = PW.SinglePlot(self.bottomPanel,plotData = self.scatteringObj.WAS)
        else:
            self.plotWidget = PW.SinglePlot(self.bottomPanel,plotData = self.scatteringObj)
        
        extremes = self.plotWidget.ax_x_lim()
        self.min_limit = min(extremes)
        self.max_limit = max(extremes)
        self.bottomLayout.addWidget(self.plotWidget)
        
        self.min_label = QtWidgets.QLabel('Min: ',self.bottomPanel)
        self.min_ledit = QtWidgets.QLineEdit(self.bottomPanel)
        self.min_ledit.setText('0')
        #self.min_ledit.set
        self.min_ledit.editingFinished.connect(self.min_limit_changed)
        
        self.max_label = QtWidgets.QLabel('Max: ',self.bottomPanel)
        self.max_ledit = QtWidgets.QLineEdit(self.bottomPanel)
        self.max_ledit.setText('inf')
        self.max_ledit.editingFinished.connect(self.max_limit_changed)
        
        self.ctrlLayout.addStretch(2)
        self.ctrlLayout.addWidget(self.min_label)
        self.ctrlLayout.addWidget(self.min_ledit)
        self.ctrlLayout.addStretch(1)
        self.ctrlLayout.addWidget(self.max_label)
        self.ctrlLayout.addWidget(self.max_ledit)
        self.ctrlLayout.addStretch(2)
        self.bottomLayout.addLayout(self.ctrlLayout)
        self.bottomPanel.setLayout(self.bottomLayout)
        #self.add_line()
        
    def add_line(self):
        self.plotWidget.plotWidget.plot_canvas.axes.axvline(0.1, color ='r')
        self.plotWidget.plotWidget.plot_canvas.draw()
        
    def createTopPanel(self):
        if self.angleRange == 'SAS':
            self.singGaussPanel = Sing_Gauss_Params(self.topPanel)
            self.doubleGaussPanel = Double_Gauss_Params(self.topPanel)
            self.singSchultzPanel = Sing_Schultz_Params(self.topPanel)
            self.doubleSchultzPanel = Double_Schultz_Params(self.topPanel)
            self.expMinimPanel = Exp_Minim_Params(self.topPanel)
            i = self.topPanel.addWidget(self.singGaussPanel)
            self.topPanel.addWidget(self.doubleGaussPanel)
            self.topPanel.addWidget(self.singSchultzPanel)
            self.topPanel.addWidget(self.doubleSchultzPanel)
            self.topPanel.addWidget(self.expMinimPanel)
            
        self.topPanel.setCurrentIndex(i)
            #self.doubleGaussPanel = QtWidgets.Panel()
        
    def min_limit_changed(self):

        
        new_min = float(self.min_ledit.text())

        if new_min < self.max_fit_limit:
            self.min_fit_limit = new_min
            
            if self.min_line is not None:
                
                self.plotWidget.remove_line(self.min_line)
                self.min_line = None
            if self.min_fit_limit > self.min_limit:
                
                self.min_line = self.plotWidget.plotWidget.plot_canvas.axes.axvline(self.min_fit_limit, color ='r')
                self.plotWidget.plotWidget.plot_canvas.draw()
                #self.min_line = self.plotWidget.plotWidget.plot_canvas.axes.axvline(x=0.1)  
                #self.plotWidget.plotWidget.redraw()
             
    def max_limit_changed(self):
        new_max = float(self.max_ledit.text())
        
        if new_max > self.min_fit_limit:
            self.max_fit_limit = new_max
            if self.max_line is not None:
                self.plotWidget.remove_line(self.max_line)
                self.max_line = None
            if self.max_fit_limit < self.max_limit:
                self.max_line = self.plotWidget.plotWidget.plot_canvas.axes.axvline(self.max_fit_limit, color ='r')
                self.plotWidget.plotWidget.plot_canvas.draw()
       
class Sing_Gauss_Params(QtWidgets.QScrollArea):
    
    def __init__(self, parent):
        super(Sing_Gauss_Params, self).__init__(parent)
        self.scrollAreaContents = QtWidgets.QWidget()
        self.gridLayout = QtWidgets.QGridLayout(self)
        self.leftLayout = QtWidgets.QFormLayout()
        
        self.rightLayout = QtWidgets.QFormLayout()
        
        self.avgRadius_lbl = QtWidgets.QLabel('Avg. Radius (nm)')
        self.avgRadius_ledit = MinMax_LineEdit(self)
        
        self.Sigma_lbl = QtWidgets.QLabel('Sigma (nm)')
        self.Sigma_ledit = MinMax_LineEdit(self)
        
        self.Intensity_lbl = QtWidgets.QLabel('Intensity')
        self.Intensity_ledit = MinMax_LineEdit(self)
        
        self.varyBckg_lbl = QtWidgets.QLabel('Vary Background')
        self.varyBckg_ckbx = QtWidgets.QCheckBox()
        self.varyBckg_ckbx.stateChanged.connect(self.bckg_state_changed)
        
        self.Bckg_lbl = QtWidgets.QLabel('Background')
        self.Bckg_lbl.setEnabled(False)
        self.Bckg_ledit = MinMax_LineEdit(self)
        self.Bckg_ledit.setEnabled(False)
        
        self.separator = QtWidgets.QFrame()
        self.separator.setGeometry(QtCore.QRect(0,0,2,10))
        self.separator.setFrameShape(QtWidgets.QFrame.VLine)
        self.separator.setFrameShadow(QtWidgets.QFrame.Sunken)
        
        self.gridLayout.addWidget(self.avgRadius_lbl,0,0,QtCore.Qt.AlignBottom)
        self.gridLayout.addWidget(self.avgRadius_ledit,0,1,QtCore.Qt.AlignBottom)
        
        self.gridLayout.addWidget(self.Sigma_lbl,1,0,QtCore.Qt.AlignBottom)
        self.gridLayout.addWidget(self.Sigma_ledit,1,1,QtCore.Qt.AlignBottom)
        
        self.gridLayout.addWidget(self.separator,0,2,-1,1, QtCore.Qt.AlignHCenter)
        
        self.gridLayout.addWidget(self.Intensity_lbl,0,3,QtCore.Qt.AlignBottom)
        self.gridLayout.addWidget(self.Intensity_ledit,0,4,QtCore.Qt.AlignBottom)
        
        self.gridLayout.addWidget(self.varyBckg_lbl,1,3,QtCore.Qt.AlignBottom)
        self.gridLayout.addWidget(self.varyBckg_ckbx,1,4,QtCore.Qt.AlignBottom)
        
        self.gridLayout.addWidget(self.Bckg_lbl,2,3,QtCore.Qt.AlignBottom)
        self.gridLayout.addWidget(self.Bckg_ledit,2,4,QtCore.Qt.AlignBottom)
        
        self.gridLayout.setSizeConstraint(QtWidgets.QLayout.SetMinAndMaxSize)
        self.scrollAreaContents.setLayout(self.gridLayout)
        self.setWidget(self.scrollAreaContents)
        
    def bckg_state_changed(self):
        
        if self.varyBckg_ckbx.isChecked():
            self.Bckg_lbl.setEnabled(True)
            self.Bckg_ledit.setEnabled(True)
        else:
            self.Bckg_lbl.setEnabled(False)
            self.Bckg_ledit.setEnabled(False)

    def return_params(self):
        params = {}
        minRadius, Radius, maxRadius = self.avgRadius_ledit.getValues()
        if minRadius is not None:
            params['R_min'] = minRadius
        if maxRadius is not None:
            params['R_max'] = maxRadius
        if Radius is not None:
            params['R_av'] = Radius
            
        minSigma, Sigma, maxSigma = self.Sigma_ledit.getValues()
        if minSigma is not None:
            params['sigma_min'] = minSigma
        if maxSigma is not None:
            params['sigma_max'] = maxSigma
        if Sigma is not None:
            params['sigma'] = Sigma
            
        minInt, Int, maxInt = self.Intensity_ledit.getValues()
        if minInt is not None:
            params['I0_min'] = minInt
        if maxInt is not None:
            params['I0_max'] = maxInt
        if Int is not None:
            params['I0'] = Int
            
        if  self.varyBckg_ckbx.isChecked():
            params['bckg_vary'] = True
            minBckg, Bckg, maxBckg = self.Bckg_ledit.getValues()
            if minBckg is not None:
                params['bckg_min'] = minBckg
            if maxBckg is not None:
                params['bckg_max'] = maxBckg
            if Bckg is not None:
                params['bckg'] = Bckg
        return params
        print 'This function should return the parameters for the Sing_Gauss_Params panel'

class Double_Gauss_Params(QtWidgets.QScrollArea):
    
    def __init__(self, parent):
        super(Double_Gauss_Params, self).__init__(parent)
        self.scrollAreaContents = QtWidgets.QWidget()
        self.gridLayout = QtWidgets.QGridLayout(self)
        self.leftLayout = QtWidgets.QFormLayout()
        
        self.rightLayout = QtWidgets.QFormLayout()
        
        self.avgRadius_lbl = QtWidgets.QLabel('1st Avg. Radius (nm)')
        self.avgRadius_ledit = MinMax_LineEdit(self)
        
        self.avgRadius2_lbl = QtWidgets.QLabel('2nd Avg. Radius (nm)')
        self.avgRadius2_ledit = MinMax_LineEdit(self)
        
        self.Sigma_lbl = QtWidgets.QLabel('1st Sigma (nm)')
        self.Sigma_ledit = MinMax_LineEdit(self)
        
        self.Sigma2_lbl = QtWidgets.QLabel('2nd Sigma (nm)')
        self.Sigma2_ledit = MinMax_LineEdit(self)
        
        self.Intensity_lbl = QtWidgets.QLabel('Intensity')
        self.Intensity_ledit = MinMax_LineEdit(self)
        
        self.varyBckg_lbl = QtWidgets.QLabel('Vary Background')
        self.varyBckg_ckbx = QtWidgets.QCheckBox()
        self.varyBckg_ckbx.stateChanged.connect(self.bckg_state_changed)
        
        self.Bckg_lbl = QtWidgets.QLabel('Background')
        self.Bckg_lbl.setEnabled(False)
        self.Bckg_ledit = MinMax_LineEdit(self)
        self.Bckg_ledit.setEnabled(False)
        
        self.separator = QtWidgets.QFrame()
        self.separator.setGeometry(QtCore.QRect(0,0,2,10))
        self.separator.setFrameShape(QtWidgets.QFrame.VLine)
        self.separator.setFrameShadow(QtWidgets.QFrame.Sunken)
        
        self.gridLayout.addWidget(self.avgRadius_lbl,0,0,QtCore.Qt.AlignBottom)
        self.gridLayout.addWidget(self.avgRadius_ledit,0,1,QtCore.Qt.AlignBottom)
        
        self.gridLayout.addWidget(self.avgRadius2_lbl,1,0,QtCore.Qt.AlignBottom)
        self.gridLayout.addWidget(self.avgRadius2_ledit,1,1,QtCore.Qt.AlignBottom)
        
        self.gridLayout.addWidget(self.Sigma_lbl,2,0,QtCore.Qt.AlignBottom)
        self.gridLayout.addWidget(self.Sigma_ledit,2,1,QtCore.Qt.AlignBottom)
        
        self.gridLayout.addWidget(self.Sigma2_lbl,3,0,QtCore.Qt.AlignBottom)
        self.gridLayout.addWidget(self.Sigma2_ledit,3,1,QtCore.Qt.AlignBottom)
        
        self.gridLayout.addWidget(self.separator,0,2,-1,1, QtCore.Qt.AlignHCenter)
        
        self.gridLayout.addWidget(self.Intensity_lbl,0,3,QtCore.Qt.AlignBottom)
        self.gridLayout.addWidget(self.Intensity_ledit,0,4,QtCore.Qt.AlignBottom)
        
        self.gridLayout.addWidget(self.varyBckg_lbl,1,3,QtCore.Qt.AlignBottom)
        self.gridLayout.addWidget(self.varyBckg_ckbx,1,4,QtCore.Qt.AlignBottom)
        
        self.gridLayout.addWidget(self.Bckg_lbl,2,3,QtCore.Qt.AlignBottom)
        self.gridLayout.addWidget(self.Bckg_ledit,2,4,QtCore.Qt.AlignBottom)
        
        self.gridLayout.setSizeConstraint(QtWidgets.QLayout.SetMinAndMaxSize)
        self.scrollAreaContents.setLayout(self.gridLayout)
        self.setWidget(self.scrollAreaContents)
        
    def bckg_state_changed(self):
        
        if self.varyBckg_ckbx.isChecked():
            self.Bckg_lbl.setEnabled(True)
            self.Bckg_ledit.setEnabled(True)
        else:
            self.Bckg_lbl.setEnabled(False)
            self.Bckg_ledit.setEnabled(False)

    def return_params(self):
        params = {}
        minRadius, Radius, maxRadius = self.avgRadius_ledit.getValues()
        if minRadius is not None:
            params['R1_min'] = minRadius
        if maxRadius is not None:
            params['R1_max'] = maxRadius
        if Radius is not None:
            params['R1_av'] = Radius
            
        minRadius, Radius, maxRadius = self.avgRadius2_ledit.getValues()
        if minRadius is not None:
            params['R2_min'] = minRadius
        if maxRadius is not None:
            params['R2_max'] = maxRadius
        if Radius is not None:
            params['R2_av'] = Radius
            
        minSigma, Sigma, maxSigma = self.Sigma_ledit.getValues()
        if minSigma is not None:
            params['sigma1_min'] = minSigma
        if maxSigma is not None:
            params['sigma1_max'] = maxSigma
        if Sigma is not None:
            params['sigma1'] = Sigma
        
        minSigma, Sigma, maxSigma = self.Sigma2_ledit.getValues()
        if minSigma is not None:
            params['sigma2_min'] = minSigma
        if maxSigma is not None:
            params['sigma2_max'] = maxSigma
        if Sigma is not None:
            params['sigma2'] = Sigma
            
        minInt, Int, maxInt = self.Intensity_ledit.getValues()
        if minInt is not None:
            params['I0_min'] = minInt
        if maxInt is not None:
            params['I0_max'] = maxInt
        if Int is not None:
            params['I0'] = Int
            
        if  self.varyBckg_ckbx.isChecked():
            params['bckg_vary'] = True
            minBckg, Bckg, maxBckg = self.Bckg_ledit.getValues()
            if minBckg is not None:
                params['bckg_min'] = minBckg
            if maxBckg is not None:
                params['bckg_max'] = maxBckg
            if Bckg is not None:
                params['bckg'] = Bckg
        return params
        print 'This function should return the parameters for the Sing_Gauss_Params panel'

class Sing_Schultz_Params(QtWidgets.QScrollArea):
    
    def __init__(self, parent):
        super(Sing_Schultz_Params, self).__init__(parent)
        self.scrollAreaContents = QtWidgets.QWidget()
        self.gridLayout = QtWidgets.QGridLayout(self)
        self.leftLayout = QtWidgets.QFormLayout()
        
        self.rightLayout = QtWidgets.QFormLayout()
        
        self.avgRadius_lbl = QtWidgets.QLabel('Avg. Radius (nm)')
        self.avgRadius_ledit = MinMax_LineEdit(self)
        
        self.Z_lbl = QtWidgets.QLabel('Z')
        self.Z_ledit = MinMax_LineEdit(self)
        
        self.Intensity_lbl = QtWidgets.QLabel('Intensity')
        self.Intensity_ledit = MinMax_LineEdit(self)
        
        self.varyBckg_lbl = QtWidgets.QLabel('Vary Background')
        self.varyBckg_ckbx = QtWidgets.QCheckBox()
        self.varyBckg_ckbx.stateChanged.connect(self.bckg_state_changed)
        
        self.Bckg_lbl = QtWidgets.QLabel('Background')
        self.Bckg_lbl.setEnabled(False)
        self.Bckg_ledit = MinMax_LineEdit(self)
        self.Bckg_ledit.setEnabled(False)
        
        self.separator = QtWidgets.QFrame()
        self.separator.setGeometry(QtCore.QRect(0,0,2,10))
        self.separator.setFrameShape(QtWidgets.QFrame.VLine)
        self.separator.setFrameShadow(QtWidgets.QFrame.Sunken)
        
        self.gridLayout.addWidget(self.avgRadius_lbl,0,0,QtCore.Qt.AlignBottom)
        self.gridLayout.addWidget(self.avgRadius_ledit,0,1,QtCore.Qt.AlignBottom)
        
        self.gridLayout.addWidget(self.Z_lbl,1,0,QtCore.Qt.AlignBottom)
        self.gridLayout.addWidget(self.Z_ledit,1,1,QtCore.Qt.AlignBottom)
        
        self.gridLayout.addWidget(self.separator,0,2,-1,1, QtCore.Qt.AlignHCenter)
        
        self.gridLayout.addWidget(self.Intensity_lbl,0,3,QtCore.Qt.AlignBottom)
        self.gridLayout.addWidget(self.Intensity_ledit,0,4,QtCore.Qt.AlignBottom)
        
        self.gridLayout.addWidget(self.varyBckg_lbl,1,3,QtCore.Qt.AlignBottom)
        self.gridLayout.addWidget(self.varyBckg_ckbx,1,4,QtCore.Qt.AlignBottom)
        
        self.gridLayout.addWidget(self.Bckg_lbl,2,3,QtCore.Qt.AlignBottom)
        self.gridLayout.addWidget(self.Bckg_ledit,2,4,QtCore.Qt.AlignBottom)
        
        self.gridLayout.setSizeConstraint(QtWidgets.QLayout.SetMinAndMaxSize)
        self.scrollAreaContents.setLayout(self.gridLayout)
        self.setWidget(self.scrollAreaContents)
        
    def bckg_state_changed(self):
        
        if self.varyBckg_ckbx.isChecked():
            self.Bckg_lbl.setEnabled(True)
            self.Bckg_ledit.setEnabled(True)
        else:
            self.Bckg_lbl.setEnabled(False)
            self.Bckg_ledit.setEnabled(False)

    def return_params(self):
        params = {}
        minRadius, Radius, maxRadius = self.avgRadius_ledit.getValues()
        if minRadius is not None:
            params['R_min'] = minRadius
        if maxRadius is not None:
            params['R_max'] = maxRadius
        if Radius is not None:
            params['R_av'] = Radius
            
        minZ, Z, maxZ = self.Z_ledit.getValues()
        if minZ is not None:
            params['Z_min'] = minZ
        if maxZ is not None:
            params['Z_max'] = maxZ
        if Z is not None:
            params['Z'] = Z
            
        minInt, Int, maxInt = self.Intensity_ledit.getValues()
        if minInt is not None:
            params['I0_min'] = minInt
        if maxInt is not None:
            params['I0_max'] = maxInt
        if Int is not None:
            params['I0'] = Int
            
        if  self.varyBckg_ckbx.isChecked():
            params['bckg_vary'] = True
            minBckg, Bckg, maxBckg = self.Bckg_ledit.getValues()
            if minBckg is not None:
                params['bckg_min'] = minBckg
            if maxBckg is not None:
                params['bckg_max'] = maxBckg
            if Bckg is not None:
                params['bckg'] = Bckg
        return params
        print 'This function should return the parameters for the Sing_Gauss_Params panel'

class Double_Schultz_Params(QtWidgets.QScrollArea):
    
    def __init__(self, parent):
        super(Double_Schultz_Params, self).__init__(parent)
        self.scrollAreaContents = QtWidgets.QWidget()
        self.gridLayout = QtWidgets.QGridLayout(self)
        self.leftLayout = QtWidgets.QFormLayout()
        
        self.rightLayout = QtWidgets.QFormLayout()
        
        self.avgRadius_lbl = QtWidgets.QLabel('1st Avg. Radius (nm)')
        self.avgRadius_ledit = MinMax_LineEdit(self)
        
        self.avgRadius2_lbl = QtWidgets.QLabel('2nd Avg. Radius (nm)')
        self.avgRadius2_ledit = MinMax_LineEdit(self)
        
        self.Z_lbl = QtWidgets.QLabel('1st Z (nm)')
        self.Z_ledit = MinMax_LineEdit(self)
        
        self.Z2_lbl = QtWidgets.QLabel('2nd Z (nm)')
        self.Z2_ledit = MinMax_LineEdit(self)
        
        self.Intensity_lbl = QtWidgets.QLabel('Intensity')
        self.Intensity_ledit = MinMax_LineEdit(self)
        
        self.varyBckg_lbl = QtWidgets.QLabel('Vary Background')
        self.varyBckg_ckbx = QtWidgets.QCheckBox()
        self.varyBckg_ckbx.stateChanged.connect(self.bckg_state_changed)
        
        self.Bckg_lbl = QtWidgets.QLabel('Background')
        self.Bckg_lbl.setEnabled(False)
        self.Bckg_ledit = MinMax_LineEdit(self)
        self.Bckg_ledit.setEnabled(False)
        
        self.separator = QtWidgets.QFrame()
        self.separator.setGeometry(QtCore.QRect(0,0,2,10))
        self.separator.setFrameShape(QtWidgets.QFrame.VLine)
        self.separator.setFrameShadow(QtWidgets.QFrame.Sunken)
        
        self.gridLayout.addWidget(self.avgRadius_lbl,0,0,QtCore.Qt.AlignBottom)
        self.gridLayout.addWidget(self.avgRadius_ledit,0,1,QtCore.Qt.AlignBottom)
        
        self.gridLayout.addWidget(self.avgRadius2_lbl,1,0,QtCore.Qt.AlignBottom)
        self.gridLayout.addWidget(self.avgRadius2_ledit,1,1,QtCore.Qt.AlignBottom)
        
        self.gridLayout.addWidget(self.Z_lbl,2,0,QtCore.Qt.AlignBottom)
        self.gridLayout.addWidget(self.Z_ledit,2,1,QtCore.Qt.AlignBottom)
        
        self.gridLayout.addWidget(self.Z2_lbl,3,0,QtCore.Qt.AlignBottom)
        self.gridLayout.addWidget(self.Z2_ledit,3,1,QtCore.Qt.AlignBottom)
        
        self.gridLayout.addWidget(self.separator,0,2,-1,1, QtCore.Qt.AlignHCenter)
        
        self.gridLayout.addWidget(self.Intensity_lbl,0,3,QtCore.Qt.AlignBottom)
        self.gridLayout.addWidget(self.Intensity_ledit,0,4,QtCore.Qt.AlignBottom)
        
        self.gridLayout.addWidget(self.varyBckg_lbl,1,3,QtCore.Qt.AlignBottom)
        self.gridLayout.addWidget(self.varyBckg_ckbx,1,4,QtCore.Qt.AlignBottom)
        
        self.gridLayout.addWidget(self.Bckg_lbl,2,3,QtCore.Qt.AlignBottom)
        self.gridLayout.addWidget(self.Bckg_ledit,2,4,QtCore.Qt.AlignBottom)
        
        self.gridLayout.setSizeConstraint(QtWidgets.QLayout.SetMinAndMaxSize)
        self.scrollAreaContents.setLayout(self.gridLayout)
        self.setWidget(self.scrollAreaContents)
        
    def bckg_state_changed(self):
        
        if self.varyBckg_ckbx.isChecked():
            self.Bckg_lbl.setEnabled(True)
            self.Bckg_ledit.setEnabled(True)
        else:
            self.Bckg_lbl.setEnabled(False)
            self.Bckg_ledit.setEnabled(False)

    def return_params(self):
        params = {}
        minRadius, Radius, maxRadius = self.avgRadius_ledit.getValues()
        if minRadius is not None:
            params['R1_min'] = minRadius
        if maxRadius is not None:
            params['R1_max'] = maxRadius
        if Radius is not None:
            params['R1_av'] = Radius
            
        minRadius, Radius, maxRadius = self.avgRadius2_ledit.getValues()
        if minRadius is not None:
            params['R2_min'] = minRadius
        if maxRadius is not None:
            params['R2_max'] = maxRadius
        if Radius is not None:
            params['R2_av'] = Radius
            
        minZ, Z, maxZ = self.Z_ledit.getValues()
        if minZ is not None:
            params['Z1_min'] = minZ
        if maxZ is not None:
            params['Z1_max'] = maxZ
        if Z is not None:
            params['Z1'] = Z
        
        minZ, Z, maxZ = self.Z2_ledit.getValues()
        if minZ is not None:
            params['Z2_min'] = minZ
        if maxZ is not None:
            params['Z2_max'] = maxZ
        if Z is not None:
            params['Z2'] = Z
            
        minInt, Int, maxInt = self.Intensity_ledit.getValues()
        if minInt is not None:
            params['I0_min'] = minInt
        if maxInt is not None:
            params['I0_max'] = maxInt
        if Int is not None:
            params['I0'] = Int
            
        if  self.varyBckg_ckbx.isChecked():
            params['bckg_vary'] = True
            minBckg, Bckg, maxBckg = self.Bckg_ledit.getValues()
            if minBckg is not None:
                params['bckg_min'] = minBckg
            if maxBckg is not None:
                params['bckg_max'] = maxBckg
            if Bckg is not None:
                params['bckg'] = Bckg
        return params
        print 'This function should return the parameters for the Sing_Gauss_Params panel'

class Exp_Minim_Params(QtWidgets.QScrollArea):
    
    def __init__(self, parent):
        super(Exp_Minim_Params, self).__init__(parent)
        self.scrollAreaContents = QtWidgets.QWidget()
        self.gridLayout = QtWidgets.QGridLayout(self)
        self.leftLayout = QtWidgets.QFormLayout()
        
        self.rightLayout = QtWidgets.QFormLayout()
        
        self.Eps_lbl = QtWidgets.QLabel('Eps')
        self.Eps_ledit = QtWidgets.QLineEdit(self)
        
        self.K_lbl = QtWidgets.QLabel('K')
        self.K_ledit = QtWidgets.QLineEdit(self)
        
        self.numbEl_lbl = QtWidgets.QLabel('Number Elements')
        self.numbEl_ledit = QtWidgets.QLineEdit(self)
        
        self.Iter_lbl = QtWidgets.QLabel('Max. Iter')
        self.Iter_ledit = QtWidgets.QLineEdit(self)
        
        self.varyBckg_lbl = QtWidgets.QLabel('Vary Background')
        self.varyBckg_ckbx = QtWidgets.QCheckBox()
        self.varyBckg_ckbx.stateChanged.connect(self.bckg_state_changed)
        
        self.Bckg_lbl = QtWidgets.QLabel('Background')
        self.Bckg_lbl.setEnabled(False)
        self.Bckg_ledit = QtWidgets.QLineEdit(self)
        self.Bckg_ledit.setEnabled(False)
        
        self.separator = QtWidgets.QFrame()
        self.separator.setGeometry(QtCore.QRect(0,0,2,10))
        self.separator.setFrameShape(QtWidgets.QFrame.VLine)
        self.separator.setFrameShadow(QtWidgets.QFrame.Sunken)
        
        self.gridLayout.addWidget(self.Eps_lbl, 0, 0, QtCore.Qt.AlignBottom)
        self.gridLayout.addWidget(self.Eps_ledit, 0, 1, QtCore.Qt.AlignBottom)
        
        self.gridLayout.addWidget(self.K_lbl, 1, 0, QtCore.Qt.AlignBottom)
        self.gridLayout.addWidget(self.K_ledit, 1, 1, QtCore.Qt.AlignBottom)
        
        self.gridLayout.addWidget(self.numbEl_lbl, 2, 0, QtCore.Qt.AlignBottom)
        self.gridLayout.addWidget(self.numbEl_ledit, 2, 1, QtCore.Qt.AlignBottom)
        
        self.gridLayout.addWidget(self.separator,0,2,-1,1, QtCore.Qt.AlignHCenter)
        
        self.gridLayout.addWidget(self.Iter_lbl,0,3,QtCore.Qt.AlignBottom)
        self.gridLayout.addWidget(self.Iter_ledit,0,4,QtCore.Qt.AlignBottom)
        
        self.gridLayout.addWidget(self.varyBckg_lbl,1,3,QtCore.Qt.AlignBottom)
        self.gridLayout.addWidget(self.varyBckg_ckbx,1,4,QtCore.Qt.AlignBottom)
        
        self.gridLayout.addWidget(self.Bckg_lbl,2,3,QtCore.Qt.AlignBottom)
        self.gridLayout.addWidget(self.Bckg_ledit,2,4,QtCore.Qt.AlignBottom)
        
        self.gridLayout.setSizeConstraint(QtWidgets.QLayout.SetMinAndMaxSize)
        self.scrollAreaContents.setLayout(self.gridLayout)
        self.setWidget(self.scrollAreaContents)
        
    def bckg_state_changed(self):
        
        if self.varyBckg_ckbx.isChecked():
            self.Bckg_lbl.setEnabled(True)
            self.Bckg_ledit.setEnabled(True)
        else:
            self.Bckg_lbl.setEnabled(False)
            self.Bckg_ledit.setEnabled(False)

    def return_params(self):
        params = {}
        if self.Eps_ledit.text() != '':
            params['eps'] = float(self.Eps_ledit.text())
        if self.K_ledit.text() != '':
            params['k'] = float(self.K_ledit.text())
        if self.numbEl_ledit.text() != '':
            params['numbElements'] = float(self.numbEl_ledit.text())
        if self.Iter_ledit.text() != '':
            params['maxIter'] = float(self.Iter_ledit.text())
            
        return params
        print 'This function should return the parameters for the Sing_Gauss_Params panel'

class MinMax_LineEdit(QtWidgets.QWidget):
    
    def __init__(self,parent):
        super(MinMax_LineEdit, self).__init__(parent)
        self.layout = QtWidgets.QGridLayout()
        self.width = 50
        self.height = 20
        
        self.lblMin = QtWidgets.QLabel('Min:')
        self.lblMax = QtWidgets.QLabel('Max:')
        self.leditMin = QtWidgets.QLineEdit()
        self.leditMin.setFixedSize(self.width,self.height)
        self.leditMax = QtWidgets.QLineEdit()
        self.leditMax.setFixedSize(self.width,self.height)
        self.leditValue = QtWidgets.QLineEdit()
        self.leditValue.setFixedSize(self.width,self.height)
        
        self.layout.addWidget(self.lblMin,0,0)
        self.layout.addWidget(self.lblMax,0,2)
        self.layout.addWidget(self.leditMin,1,0)
        self.layout.addWidget(self.leditValue,1,1)
        self.layout.addWidget(self.leditMax,1,2)
        self.layout.setContentsMargins(0,0,0,0)
        self.layout.setVerticalSpacing(0)
        self.setLayout(self.layout)
        
    def setEnabled(self, newValue):
        self.lblMin.setEnabled(newValue)
        self.lblMax.setEnabled(newValue)
        self.leditMin.setEnabled(newValue)
        self.leditMax.setEnabled(newValue)
        self.leditValue.setEnabled(newValue)
        
    def getValues(self):
        minimum = None
        maximum = None
        avg = None
        if self.leditMin.text() != '':
            minimum = float(self.leditMin.text())
        if self.leditMax.text() != '':
            maximum = float(self.leditMax.text())
        if self.leditValue.text() != '':
            avg = float(self.leditMax.text())
        return [minimum, avg, maximum]
    
class Choose_Fit_Plot(QtWidgets.QDialog):
    
    def __init__(self, parent = None, avlbAngles = 'SWAS'):
        super(Choose_Fit_Plot,self).__init__(parent)
        self.layout = QtWidgets.QVBoxLayout()
        self.buttonLayout = QtWidgets.QHBoxLayout()
        
        self.angleComboBox = QtWidgets.QComboBox()
        self.angleComboBox.addItem('Select angle range...',None)
        if avlbAngles == 'SWAS':
            self.angleComboBox.addItem('Small-angles', 'SAS')
            self.angleComboBox.addItem('Wide-angles', 'WAS')
        elif avlbAngles == 'SAS':
            self.angleComboBox.addItem('Small-angles', 'SAS')
        elif avlbAngles == 'WAS':
            self.angleComboBox.addItem('Wide-angles', 'WAS')
        self.angleComboBox.currentIndexChanged.connect(self.angleSelected)
        self.fitTypeComboBox = QtWidgets.QComboBox()
        self.fitTypeComboBox.setEnabled(False)
        self.okButton = QtWidgets.QPushButton('Ok')
        self.okButton.setEnabled(False)
        self.okButton.clicked.connect(self.accept)
        self.cancelButton = QtWidgets.QPushButton('Cancel')
        self.cancelButton.clicked.connect(self.reject)
        
        self.layout.addWidget(self.angleComboBox)
        self.layout.addWidget(self.fitTypeComboBox)
        self.buttonLayout.addWidget(self.cancelButton)
        self.buttonLayout.addWidget(self.okButton)
        self.layout.addLayout(self.buttonLayout)
        self.setLayout(self.layout)
        
        
    def angleSelected(self):
        for i in xrange(self.fitTypeComboBox.count()):
                    self.fitTypeComboBox.removeItem(0)          
        if self.angleComboBox.currentData() is not None:
            self.fitTypeComboBox.setEnabled(True)
            self.okButton.setEnabled(True)
            
            currAngle = str(self.angleComboBox.currentData())
            if currAngle == 'SAS':
                for ft in _AvlbSASFitDic:
                    self.fitTypeComboBox.addItem(ft,_AvlbSASFitDic[ft])
            elif currAngle == 'WAS':
                for ft in _AvlbWASFitDic:
                    self.fitTypeComboBox.addItem(ft,_AvlbWASFitDic[ft])
        else:
            self.fitTypeComboBox.setEnabled(False)
            self.okButton.setEnabled(False)
            
        
        
