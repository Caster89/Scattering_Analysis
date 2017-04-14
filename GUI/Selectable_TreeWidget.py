from PyQt5 import QtCore, QtGui, QtWidgets, Qt
from ..SWAS_Sequence import SWAS_Sequence
from ..SAS import SmallAngleScattering
from ..WAS import WideAngleScattering
from ..SWAS import SWAS
from ..SWAS_Sequence import SWAS_Sequence
from Dialog_Windows import Fitting_Dialog, Choose_Fit_Plot
import logging
_SAS = 0
_WAS = 1
class Scattering_TreeWidgetItem(QtWidgets.QTreeWidgetItem):
    """Class used for the tree widget panel to store the scattering object.
    Based on whether the object is SAS, WAS, SWAS or a sequence, the object
    will then decide how to provide the data. The standard display has SWAS
    sequence object sorted by angle. THe top items will therefore have 2
    children: one with all the SAS data and one with the WAS data.
    """
        
    def __init__(self, parent, scatteringObject, **kwargs):
        '''Initialize the object. The first step is to initialize the object
        using the super initializer. The scattering object is stored in a
        variable. The sample name is set as the text in position 0. The
        initialization then splits based on the scattering object past.
        If the scattering object is either SAS or WAS then nothing is done.
        If the scattering object is a SWAS object then one/two children are
        initialized based on whether the SWAS object has SAS,WAS or both.
        Finally if the object is a sequence then one/two childred are created,
        and then children are added for each curve of the sequence.
        '''
        
        super(Scattering_TreeWidgetItem, self).__init__(parent, **kwargs)
        
        self.setText(0, scatteringObject.sampleName)
        self.setFlags(self.flags() | QtCore.Qt.ItemIsTristate | QtCore.Qt.ItemIsUserCheckable)
        self.setCheckState(0,0)
        #the default sorting is based on the scattering angle (e.g. SAS vs WAS)
        
        self.sorting = None
        self.scatteringObj = scatteringObject
        
        
        #First case in which the object is a SWAS object. The item will therefore
        #have 2 children, a SAS and a WAS
        if isinstance(scatteringObject,SWAS):
            #Need to add 2 Children, a SAS and a WAS child
            self.addChild(Scattering_TreeWidgetItem(self,scatteringObject.SAS))
            self.child(_SAS).setFlags(self.child(_SAS).flags() | QtCore.Qt.ItemIsUserCheckable)
            if scatteringObject.SAS is None:
                self.child(_SAS).setFlags(self.child(_SAS).flags() ^ QtCore.Qt.ItemIsUserCheckable)
            
            self.addChild(Scattering_TreeWidgetItem(self,scatteringObject.WAS))
            self.child(_WAS).setFlags(self.child(_WAS).flags() | QtCore.Qt.ItemIsUserCheckable)
            if scatteringObject.WAS is  None:
                self.child(_WAS).setFlags(self.child(_WAS).flags() ^ QtCore.Qt.ItemIsUserCheckable)
        
        if isinstance(scatteringObject, SWAS_Sequence):
            if kwargs.get('sorting', 'angle'):
                self.sort_by_angle()
                self.sorting = 'angle'
            else:
                self.sort_by_time()
                self.sorting = 'time'

    
    def sort_by_angle(self):
        if self.sorting == 'angle':
            return
        for n in xrange(self.childCount()):
            self.removeChild(self.child(0))
        #If the object is a SWAS_Sequence object. since the standard setup
        #is with sorting by angle, each SWAS_Sequence will have 2 children
        if isinstance(self.scatteringObj,SWAS_Sequence):
            
            self.addChild(QtWidgets.QTreeWidgetItem(self,['SAXS']))
            self.child(_SAS).setFlags(self.child(_SAS).flags()| QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsTristate)
            self.child(_SAS).setCheckState(0,QtCore.Qt.Unchecked)
            
            self.addChild(QtWidgets.QTreeWidgetItem(self,['WAXS']))
            self.child(_WAS).setFlags(self.child(_SAS).flags()| QtCore.Qt.ItemIsUserCheckable|QtCore.Qt.ItemIsTristate)
            self.child(_WAS).setCheckState(0,QtCore.Qt.Unchecked)
            
            #If the SWAS_Sequence has both small angle and wide angle data
            if self.scatteringObj.avlbCurves == 'Both':
                for i in xrange(self.scatteringObj.size):
                    
                    self.child(_SAS).addChild(Scattering_TreeWidgetItem(self.child(_SAS), self.scatteringObj[i].SAS))
                    numbChildren = self.child(_SAS).childCount()
                    
                    self.child(_SAS).child(numbChildren-1).setFlags(self.child(_SAS).child(numbChildren-1).flags()| QtCore.Qt.ItemIsUserCheckable)
                    self.child(_SAS).child(numbChildren-1).setCheckState(0,QtCore.Qt.Unchecked)

                    self.child(_WAS).addChild(Scattering_TreeWidgetItem(self.child(_WAS), self.scatteringObj[i].WAS))
                    
                    self.child(_WAS).child(numbChildren-1).setFlags(self.child(_SAS).child(numbChildren-1).flags()| QtCore.Qt.ItemIsUserCheckable)
                    self.child(_WAS).child(numbChildren-1).setCheckState(0,QtCore.Qt.Unchecked)
            #If the sequence only has small angle data
            elif scatteringObject.avlbCurves == 'SAS':
                #The WAS child cannot be checked anymore
                self.child(_WAS).setFlags(self.child(_WAS).flags() ^ QtCore.Qt.ItemIsUserCheckable)
                
                for i in xrange(scatteringObject.size):
                    self.child(_SAS).addChild(Scattering_TreeWidgetItem(self.child(_SAS), self.scatteringObj[i].SAS))
                    numbChildren = self.child(_SAS).childCount()
                    
                    self.child(_SAS).child(numbChildren-1).setFlags(self.child(_SAS).child(numbChildren-1).flags()| QtCore.Qt.ItemIsUserCheckable)
                    self.child(_SAS).child(numbChildren-1).setCheckState(0,QtCore.Qt.Unchecked)  
            #If the sequence only has wide angle data
            elif scatteringObject.avlbCurves == 'WAS':
                #The SAS child cannot be checked anymore
                self.child(_SAS).setFlags(self.child(_SAS).flags() ^ QtCore.Qt.ItemIsUserCheckable)
                for i in xrange(scatteringObject.size):
                    self.child(_WAS).addChild(Scattering_TreeWidgetItem(self.child(_SAS), self.scatteringObj[i].WAS))
                    numbChildren = self.child(_WAS).childCount()
                    
                    self.child(_WAS).child(numbChildren-1).setFlags(self.child(_SAS).child(numbChildren-1).flags()| QtCore.Qt.ItemIsUserCheckable)
                    self.child(_WAS).child(numbChildren-1).setCheckState(0,QtCore.Qt.Unchecked)
        self.sorting = 'angle'
        
    def sort_by_time(self):
        if self.sorting == 'time':
            return
        for n in xrange(self.childCount()):
            self.removeChild(self.child(0))
        #If the object is a SWAS_Sequence object. The item will have as many
        #children as the number of scattering objects in the sequence 
        if isinstance(self.scatteringObj,SWAS_Sequence):
            numbObj = self.scatteringObj.size
            for n in xrange(numbObj):
                self.addChild(Scattering_TreeWidgetItem(self, self.scatteringObj[n]))
        self.sorting = 'time'
    
    def checked_children(self):
        '''checked_children checks the children of the item and finds the ones which have
        a check state higher than 0 (fully or partially selected). It then returns a dict.
        with the position of the selected items.  
        '''
    
        checkedChild = {}
        if isinstance(self.scatteringObj, SWAS_Sequence):
            checkedChild['SAS'] = []
            checkedChild['WAS'] = []
            if self.sorting == 'angle':
                numbChild = self.child(_SAS).childCount()
                for n in xrange(numbChild):
                    if self.child(_SAS).child(n).checkState(0) > 0:
                        checkedChild['SAS'].append(n)
                    if self.child(_WAS).child(n).checkState(0) > 0:
                        checkedChild['WAS'].append(n)
            if self.sorting =='time':
                numbChild = self.childCount()
                for n in xrange(numbChild):
                    if self.child(n).child(_SAS).checkState(0) > 0:
                        checkedChild['SAS'].append(n)
                for n in xrange(numbChild):
                    if self.child(n).child(_WAS).checkState(0) > 0:
                        checkedChild['WAS'].append(n)
        if isinstance(self.scatteringObj, SWAS):
            checkedChild['SAS'] = []
            checkedChild['WAS'] = []
            if self.child(_SAS).checkState(0) > 0:
                checkedChild['SAS'].append(0)
            if self.child(_WAS).checkState(0) > 0:
                checkedChild['WAS'].append(0)
        return checkedChild
    
                
    def set_sorting(self, sorting = 'angle'):
        if sorting == 'angle':
            self.sort_by_angle()
        elif sorting == 'time':
            self.sort_by_time()
        else:
            print '{} is an unrecognized sorting order'.format(sorting)
    
class TreeWidget_Header (QtWidgets.QHeaderView):
    
    def __init__(self, parent):
        super(TreeWidget_Header,self).__init__(QtCore.Qt.Horizontal, parent)
        self.layout = QtWidgets.QHBoxLayout(self)
        self.angleSort = QtWidgets.QPushButton('Angle')
        self.angleSort.clicked.connect(self.sort_by_angle)
        self.timeSort = QtWidgets.QPushButton('Time')
        self.timeSort.clicked.connect(self.sort_by_time)
        self.layout.addWidget(self.angleSort)
        self.layout.addWidget(self.timeSort)
        self.setLayout(self.layout)
        self.setSectionResizeMode(QtWidgets.QHeaderView.Stretch)
    
    def sizeHint(self):
        baseSize = super(TreeWidget_Header,self).sizeHint()
        baseSize.setHeight(30)
        return baseSize
    
    def sort_by_angle(self):
        self.parent().sort_by_angle()
    def sort_by_time(self):
        self.parent().sort_by_time()
        
class Selectable_TreeWidget(QtWidgets.QTreeWidget):
    
    def __init__(self, parent, **kwargs):
        if kwargs.get('verbose', False):
			self.logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
        else:
			self.logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
			
            
        super(Selectable_TreeWidget, self).__init__(parent)
        self.setHeader(TreeWidget_Header(self))
    
    def add_scattering_object(self, SO):
        self.addTopLevelItem(Scattering_TreeWidgetItem(self,SO))
        #for n in xrange(self.topLevelItemCount()):
            #self.itemAt(0,n).sort_by_time()
    def sort_by_angle(self):
        nTopLvl = self.topLevelItemCount()
        for n in xrange(nTopLvl):
            self.itemAt(0,n).sort_by_angle()
    
    def sort_by_time(self):
        nTopLvl = self.topLevelItemCount()
        for n in xrange(nTopLvl):
            self.itemAt(0,n).sort_by_time()
    
    def fit(self, plot_widget):
        '''fit is used to fit the currently selected scattering objects.
        If first displays a fitting dialog widget to allow the user to select
        the type of fit as well as impose limits on the parameters
        '''
        checkedItems = self.checked_items()
        checkedItem = checkedItems[0]
        
        #Ask the user which angles (small o wide)
        #should be fitted.
        #TODO this window should only appear if a SWAS or SWAS_Sequence object is selected
        if isinstance(checkedItem['Object'], SAS):
            angles = 'SAS'
        elif isinstance(checkItem['Object'], WAS):
            angles = 'WAS'
        else:
            chooseFitAngle = QtWidgets.QDialog(self)
            fitTypeCombo = QtWidgets.QComboBox()
            fitTypeCombo.addItem('Small Angles', 'SAS')
            fitTypeCombo.addItem('Wide Angles', 'WAS')
            okButton = QtWidgets.QPushButton('OK')
            cancelButton = QtWidgets.QPushButton('Cancel')
            dialogLayout = QtWidgets.QVBoxLayout()
            dialogLayout.addWidget(fitTypeCombo)
            buttonLayout = QtWidgets.QHBoxLayout()
            buttonLayout.addWidget(cancelButton)
            buttonLayout.addWidget(okButton)
            dialogLayout.addLayout(buttonLayout)
            okButton.clicked.connect(chooseFitAngle.accept)
            cancelButton.clicked.connect(chooseFitAngle.reject)
            chooseFitAngle.setLayout(dialogLayout)
            if chooseFitAngle.exec_():
                angles = fitTypeCombo.currentData
            self.logging.debug('Selected {} angles'.format(angles))
        
        diagResult = []
        dialog = Fitting_Dialog(plot_widget, checkedItem, angleRange = angles)
        #dialog.finished.connect(dialog.result)
        if dialog.exec_():
            fitType, params = dialog.result()
            self.logging.debug('dialog result:\n -Fit Type:\n {}\n -Params:\n {} '.format(fitType, params))

        
        if angles == 'SAS':
            self.logging.debug('selected SAS')
            if isinstance(checkedItem['Object'], SAS):
                self.logging.debug('fitting SAS')
                checkedItem['Object'].fit_data(fitType, **params)
            elif isinstance(checkedItem['Object'], SWAS):
                self.logging.debug('fitting SAS from SWAS')
                checkedItem['Object'].fit_SAS(fitType, **params)
        if angles == 'WAS':
            self.logging.debug('selected WAS')
            if isinstance(checkedItem['Object'], WAS):
                self.logging.debug('fitting WAS')
                checkedItem['Object'].fit_peaks(fitType, **params)
            elif isinstance(checkedItem['Object'], SWAS):
                self.logging.debug('fitting WAS from SWAS')
                checkedItem['Object'].fit_WAS(fitType, **params)
        
        
    def plot_data(self, plot_widget):
        '''plot is used to pass the selected data to the plot_widget to plot
        based on how the selection was made
        '''
        checkedItems = self.checked_items()
        for ci in checkedItems:
            plot_widget.add_plot(ci)
        
    
    def plot_fit(self,plot_widget):
        checkedItems = self.checked_items()
        if len(checkedItems) > 0:
            checkedItem = checkedItems[0]
        else:
            return 0
        chooseFitAngle = Choose_Fit_Plot(self)
        if chooseFitAngle.exec_():
            fitType = str(chooseFitAngle.fitTypeComboBox.currentData())
            plot_widget.add_fit(checkedItem, fitType)
        
        
        
        #dialog = Fitting_Dialog(plot_widget, checked_Item, diagResult)
        #dialog.finished.connect(dialog.result)
        #if dialog.exec_():
            #fit_type, params = dialog.result()
        
    def fitDialog(self, params):
        print params
    
    def checked_items(self):
        checkedItems = []
        nTopLvl = self.topLevelItemCount()
        selectedItems = []
        for i in range(nTopLvl):
            currItem = self.itemAt(0,i)
            currStatus = currItem.checkState(0)
            if currStatus > 0:
                selectedItems.append({'Object' : currItem.scatteringObj})
                selectedItems[-1].update(currItem.checked_children())
        return selectedItems
    