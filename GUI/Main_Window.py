from PyQt5 import QtCore, QtGui, QtWidgets
from Selectable_TreeWidget import Selectable_TreeWidget,Scattering_TreeWidgetItem
from PlottingWidgets import CentralWidget
import Dialog_Windows
class Main_Window(QtWidgets.QMainWindow):

 
    def __init__(self, **kwargs):
        super(Main_Window,self).__init__(**kwargs)

        
        self.splitter = QtWidgets.QSplitter(self)

        self.objectsList = Selectable_TreeWidget(self)
        self.objectsList.setColumnCount(1)
        self.objectsList.setHeaderLabels([''])
        #self.objectsList.itemClicked.connect(self.selection_changed)
        
        
        self.centralArea = CentralWidget(self)
        self.centralArea.setSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding,QtWidgets.QSizePolicy.MinimumExpanding)
        #QtWidgets.QMdiArea(self)
        self.optionsPanel = QtWidgets.QTreeView(self)
        
        self.splitter.addWidget(self.objectsList)
        self.splitter.addWidget(self.centralArea)
        self.splitter.addWidget(self.optionsPanel)
        self.splitter.setSizes([self.width()*1./5., self.width()+3./5., self.width()*1./5.])
        #layout = QtWidgets.QVBoxLayout(self)
        #layout.addWidget(self.splitter)
        self.setCentralWidget(self.splitter)
        self.statusBar()
        
        #Setup the menu bar
        self.menuBar()
        self.menuBar().setNativeMenuBar(False)
        
        #File menu actions
        loadFromDataAction = QtWidgets.QAction('&Load from Data', self) 
        loadFromDataAction.setShortcut('Ctrl+L')
        loadFromDataAction.setStatusTip('Loads a text file with the data')
        
        loadFromPickleAction = QtWidgets.QAction('&Load from Pickle', self) 
        loadFromDataAction.setStatusTip('Loads a pickle file with the data')
        
        saveDatasetAction = QtWidgets.QAction('&Save to Pickle', self)
        saveDatasetAction.setShortcut('Ctrl+S')
        saveDatasetAction.setStatusTip('Saves dataset to pickle file')
        
        savePlotAction = QtWidgets.QAction('&Save Plot', self) 
        savePlotAction.setStatusTip('Saves plot to image')
        
        exitAction = QtWidgets.QAction('&Exit', self)        
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(QtWidgets.qApp.quit)
        
        self.fileMenu = self.menuBar().addMenu('File')
        self.fileMenu.addAction(loadFromDataAction)
        self.fileMenu.addAction(loadFromPickleAction)
        self.fileMenu.addAction(saveDatasetAction)
        self.fileMenu.addAction(savePlotAction)
        self.fileMenu.addAction(exitAction)
        
        #Actions menu
        fitDataAction = QtWidgets.QAction('&Fit Data', self)
        fitDataAction.setStatusTip('Fit data from the selected object')
        fitDataAction.triggered.connect(self.fit_data)
        #fitDataAction.setEnabled(False)
        
        
        plotDataAction = QtWidgets.QAction('&Plot Data', self)
        plotDataAction.setStatusTip('Plot data from the selected object')
        plotDataAction.triggered.connect(self.plot_data)
        #plotDataAction.setEnabled(False)
        
        plotFitAction = QtWidgets.QAction('&Plot Fit', self)
        plotFitAction.setStatusTip('Plots the data fit from the selected object')
        plotFitAction.triggered.connect(self.plot_fit)
        
        
        
        
        self.actionsMenu = self.menuBar().addMenu('Actions')
        
        self.actionsMenu.addAction(fitDataAction)
        self.actionsMenu.addAction(plotDataAction)
        self.actionsMenu.addAction(plotFitAction)

        #self.fileMenu.addAction(exitAction)
        
    def fit_data(self):
        print 'initiate Fit plot'
        self.objectsList.fit(self.centralArea)

    def plot_fit(self):
        #dialog = Fitting_Dialog(window, data[0])
        #dialog.finished.connect()
        #dialog.exec_()
        self.objectsList.plot_fit(self.centralArea)
            #print self.objectsList.selectedItems()
    def plot_data(self):
        self.objectsList.plot_data(self.centralArea)
    

        
if __name__ == '__main__':

    import sys
    app = QtWidgets.QApplication(sys.argv)
    window = Main_Window()
    window.setGeometry(500, 300, 300, 300)
    window.show()
    sys.exit(app.exec_())