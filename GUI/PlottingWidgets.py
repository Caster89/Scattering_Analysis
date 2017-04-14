import matplotlib as mpl
mpl.use("Qt5Agg")
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib import lines
import logging
import numpy as np
from PyQt5 import QtCore, QtGui, QtWidgets, Qt
from ..SWAS_Sequence import SWAS_Sequence
from ..SAS import SmallAngleScattering
from ..WAS import WideAngleScattering
from ..SWAS import SWAS
from ..SWAS_Sequence import SWAS_Sequence

plot_types = ('SAS','WAS','SWAS','SAS_Fit','WAS_Fit','SWAS_Fit')
class CentralWidget(QtWidgets.QWidget):
    """Class used for the tree widget panel to store the scattering object.
    Based on whether the object is SAS WAS SWAS or a sequence the object
    will then decide how to provide the data.
    """
    def __init__(self, parent, scattering_object = None, **kwargs):
        if kwargs.get('verbose', False):
			self.logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
        else:
			self.logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
			
        super(CentralWidget, self).__init__(parent, **kwargs)
        self.layout = QtWidgets.QVBoxLayout(self)
        self.layout.setContentsMargins(0,0,0,0)
        self.stackedWidget = QtWidgets.QStackedWidget(self)
        self.toolBar = QtWidgets.QToolBar(self)
        self.selectPlot = QtWidgets.QComboBox(self.toolBar)
        self.selectPlot.activated.connect(self.selected_plot)
        self.toolBar.addWidget(self.selectPlot)
        #self.selectPlot.insertItem(0,'Test1')
        #self.selectPlot.insertItem(1,'Test2')
        self.layout.addWidget(self.stackedWidget)
        self.layout.addWidget(self.toolBar)
        self.setLayout(self.layout)
    
    def add_plot(self, plotData, plotLayout = 'individual', **kwargs):
        """add_plot is used to add a new widget to the stacked widget along
        with a new entry in the dropdown list. It decides which plotting widget
        to use based on the object to be plotted  and how it should be displayed.
        Args:
            plotData(dict): contains all the information relative to the object
                to plot:
                    'Object': one of the existing scattering objects
                    'SAS': the list of SAS objects which needs to be plotted. (Valid
                        only if the objects is a SWAS or SWAS_Sequence object)
                    'WAS': the list of WAS objects which needs to be plotted. (Valid
                        only if the objects is a SWAS or SWAS_Sequence object)
                    If WAS and SAS are both in the dictionary but are different then the
                    plot will not be done and ean error message will be printed
        """
        
       
        if isinstance(plotData['Object'], SWAS_Sequence):
            if plotData['SAS'] or plotData['WAS']:
                #Plot both the SAS and WAS patterns
                if plotData['SAS'] == plotData['WAS']:
                    if len(plotData['SAS']) == 1:
                        self.logging.debug('plotting single SWAS')
                        currPos = self.stackedWidget.addWidget(DoublePlot(self.stackedWidget, plotData = plotData['Object'][plotData['SAS'][0]],\
                                                                          plot_type = 'SWAS'))
                        self.selectPlot.insertItem(currPos, '{}_plot'.format(plotData['Object'][plotData['SAS'][0]].sampleName))
                    else:
                        self.logging.debug('plotting multiple SWAS')
                        currPos = self.stackedWidget.addWidget(MultiplePlot(self.stackedWidget, plotData = plotData, plot_type = 'SWAS_Sequence'))
                        self.selectPlot.insertItem(currPos, '{}_plot'.format(plotData['Object'].sampleName))
                #Plot only SAS data
                elif not plotData['WAS']:
                    if len(plotData['SAS']) == 1:
                        self.logging.debug('plotting single SAS')
                        currPos = self.stackedWidget.addWidget(SinglePlot(self.stackedWidget, plotData = plotData['Object'][plotData['SAS'][0]].SAS,\
                                                                          plot_type = 'SAS'))
                        self.selectPlot.insertItem(currPos, '{}_plot'.format(plotData['Object'][plotData['SAS'][0]].SAS.sampleName))
                    else:
                        self.logging.debug('plotting multiple SAS')
                        currPos = self.stackedWidget.addWidget(MultiplePlot(self.stackedWidget, plotData = plotData, plot_type = 'SAS'))
                        self.selectPlot.insertItem(currPos, '{}_plot'.format(plotData['Object'].sampleName))
                #Plot only WAS data
                elif not plotData['SAS']:
                    if len(plotData['WAS']) == 1:
                        self.logging.debug('plotting single WAS')
                        currPos = self.stackedWidget.addWidget(SinglePlot(self.stackedWidget, plotData = plotData['Object'][plotData['WAS'][0]].SAS,\
                                                                          plot_type = 'WAS'))
                        self.selectPlot.insertItem(currPos, '{}_plot'.format(plotData['Object'][plotData['WAS'][0]].WAS.sampleName))
                    else:
                        self.logging.debug('plotting multiple WAS')
                        currPos = self.stackedWidget.addWidget(MultiplePlot(self.stackedWidget, plotData = plotData, plot_type = 'WAS'))
                        self.selectPlot.insertItem(currPos, '{}_plot'.format(plotData['Object'].sampleName))
                else:
                    self.logging.error('Cannot understend how/what to plot with the given selection')
                    return
            elif isinstance(plotData['Object'], SWAS):
                if plotData['SAS'] and plotData['WAS']:
                    self.logging.debug('plotting SWAS')
                    currPos = self.stackedWidget.addWidget(DoublePlot(self.stackedWidget, plotData = plotData['Object'],\
                                                                        plot_type = 'SWAS'))
                    self.selectPlot.insertItem(currPos, '{}_plot'.format(plotData['Object'].sampleName))
                elif plotData['SAS']:
                    self.logging.debug('plotting single SAS')
                    currPos = self.stackedWidget.addWidget(SinglePlot(self.stackedWidget, plotData = plotData['Object'].SAS,\
                                                                      plot_type = 'SAS'))
                    self.selectPlot.insertItem(currPos, '{}_plot'.format(plotData['Object'].SAS.sampleName))
                elif plotData['WAS']:
                    self.logging.debug('plotting single WAS')
                    currPos = self.stackedWidget.addWidget(SinglePlot(self.stackedWidget, plotData = plotData['Object'].WAS,\
                                                                      plot_type = 'WAS'))
                    self.selectPlot.insertItem(currPos, '{}_plot'.format(plotData['Object'].WAS.sampleName))
            elif isinstance(plotData['Object'], SAS):
                self.logging.debug('plotting single SAS')
                currPos = self.stackedWidget.addWidget(SinglePlot(self.stackedWidget, plotData = plotData['Object'].SAS,\
                                                                      plot_type = 'SAS'))
            elif isinstance(plotData['Object'], WAS):
                self.logging.debug('plotting single WAS')
                currPos = self.stackedWidget.addWidget(SinglePlot(self.stackedWidget, plotData = plotData['Object'].WAS,\
                                                                      plot_type = 'WAS'))
            else:
                self.logging.info('No recognizable element selected')
                return
        
    def add_fit(self, plotData, fitType, **kwargs ):
        '''add_fit is used ot add a new plot containing the fitted data. In this case
        all SAS plots are double (the data and the distribution).
        Args:
            plotData(dict): contains all the information relative to the object for which
                the fitting was done.
            fitType(string): the name of the fitting method which shoudl be plotted
        '''

        if isinstance(plotData['Object'], SWAS_Sequence):
            if not plotData['WAS']:
                if len(plotData['SAS']) == 1:
                    self.logging.debug('plotting single SAS')
                    currPos = self.stackedWidget.addWidget(FitPlot(self.stackedWidget, plotData = plotData['Object'][plotData['SAS'][0]].SAS,\
                                                                      plot_type = 'SAS', fitType = fitType))
                    self.selectPlot.insertItem(currPos, '{}_fit'.format(plotData['Object'][plotData['SAS'][0]].SAS.sampleName))
                else:
                    self.logging.debug('plotting multiple SAS')
                    currPos = self.stackedWidget.addWidget(MultiplePlot(self.stackedWidget, plotData = plotData, plot_type = 'SAS', fitType = fitType))
                    self.selectPlot.insertItem(currPos, '{}_fit'.format(plotData['Object'].sampleName))
                #Plot only WAS data
            elif not plotData['SAS']:
                if len(plotData['WAS']) == 1:
                    self.logging.debug('plotting single WAS')
                    currPos = self.stackedWidget.addWidget(FitPlot(self.stackedWidget, plotData = plotData['Object'][plotData['WAS'][0]].SAS,\
                                                                          plot_type = 'WAS'))
                    self.selectPlot.insertItem(currPos, '{}_plot'.format(plotData['Object'][plotData['WAS'][0]].WAS.sampleName))
                else:
                    self.logging.debug('plotting multiple WAS')
                    currPos = self.stackedWidget.addWidget(MultiplePlot(self.stackedWidget, plotData = plotData, plot_type = 'WAS'))
                    self.selectPlot.insertItem(currPos, '{}_plot'.format(plotData['Object'].sampleName))
            else:
                self.logging.error('Cannot understend how/what to plot with the given selection')
                return
        elif isinstance(plotData['Object'], SWAS):
            if not plotData['WAS']:
                currPos = self.stackedWidget.addWidget(FitPlot(self.stackedWidget, plotData = plotData['Object'].SAS,\
                                                                plot_type = 'SAS', fitType = fitType))
                self.selectPlot.insertItem(currPos, '{}_fit'.format(plotData['Object'].SAS.sampleName))
            elif not plotData['SAS']:
                currPos = self.stackedWidget.addWidget(FitPlot(self.stackedWidget, plotData = plotData['Object'].SAS,\
                                                                  plot_type = 'WAS'))
                self.selectPlot.insertItem(currPos, '{}_plot'.format(plotData['Object'].WAS.sampleName))
        elif isinstance(plotData['Object'], SAS):
            currPos = self.stackedWidget.addWidget(FitPlot(self.stackedWidget, plotData = plotData['Object'],\
                                                            plot_type = 'SAS', fitType = fitType))
            self.selectPlot.insertItem(currPos, '{}_fit'.format(plotData['Object'].sampleName))
        elif isinstance(plotData['Object'], WAS):
            currPos = self.stackedWidget.addWidget(FitPlot(self.stackedWidget, plotData = plotData['Object'],\
                                                            plot_type = 'WAS', fitType = fitType))
            self.selectPlot.insertItem(currPos, '{}_fit'.format(plotData['Object'].sampleName))

        else:
            self.logging.info('No recognizable element selected')
            return
            
    def selected_plot(self):
        
        selIndx = self.selectPlot.currentIndex()
        self.stackedWidget.setCurrentIndex(selIndx)
        self.logging.debug('changed index to {}'.format(selIndx))

class SingleCanvasPlot(FigureCanvas):
    '''Drawing Canvas for plotting one axes on a single figure. Can be used to plot
    SAS or WAS data 
    '''
    def __init__(self, parent=None, figsize = (5,4), **kwargs):
                 #firstAx = [[0,1],[0,1]], secondAx = [[0,1],[1,2]] ):
        '''Initiates the canvas.
            Args:
                parent (QtWidget): the parent widget to which the canvas is associated
                    Defauts to None
                figsize (list of int): the size of the figure in the form: (width,height).
                    This will be used ot create the figure. Defaults to (5,4)
                
                
        '''
        if kwargs.get('verbose', False):
			self.logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
        else:
			self.logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
			

        self.fig = Figure(figsize=figsize)
        self.axes = self.fig.add_subplot(111)
        
        super(SingleCanvasPlot, self).__init__(self.fig)
        #FigureCanvas.__init__(self, fig)
        self.setParent(parent)

        self.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                QtWidgets.QSizePolicy.Expanding)
        self.updateGeometry()
        
class SinglePlot(QtWidgets.QWidget):
    '''Wrapping class used to place the single canvas plot in a widget
    '''
    def __init__(self, parent=None,**kwargs):
        if kwargs.get('verbose', False):
			self.logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
        else:
			self.logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
			
        super(SinglePlot, self).__init__(parent)
        self.scatteringObject = kwargs.get('plotData', None)
        self.plot_layout = QtWidgets.QVBoxLayout(self)
        self.plot_canvas = SingleCanvasPlot(self, figsize = (5,4))

        self.navi_toolbar = NavigationToolbar(self.plot_canvas, self)
        self.plot_layout.addWidget(self.plot_canvas)  # the matplotlib canvas
        self.plot_layout.addWidget(self.navi_toolbar)
        self.setLayout(self.plot_layout)
        
        #If the widget is created with an Onject than it can directly be plotted
        if self.scatteringObject is not None:
            self.plot_data()
        
    def plot_data(self):
        self.scatteringObject.plot_data(ax = self.plot_canvas.axes)
    
    def remove_line(self, line):
        '''remove_line is used ot remove a particular object from the axes given
        it's handle. Useful to remove vertical lines drawn to provide visial aid
        for fitting limits
        '''
        if isinstance(line, lines.Line2D):
            if line in self.plot_canvas.axes.lines:
                self.plot_canvas.axes.lines.remove(line)
            else:
                self.logging.error('{} is not in the axes'.format(line))
    
    def axvline(self,x):
        '''axvline simply calls axvline on the widgets axes
        '''
        self.plot_canvas.axes.axvline(x, color = 'k', linewidth = 2)
                
    def cla(self):
        '''cla cleans the axes in the widget
        '''
        self.plot_canvas.axes.cla()
    
    def redraw(self):
        '''Redraws the canvas in case something was added or removed from it
        '''
        self.plot_canvas.draw()
    
    def ax_x_lim(self):
        '''returns the x limits of the canvas' axes
        '''
        return self.plot_canvas.axes.get_xlim()
    
    def ax_y_lim(self):
        '''returns the y limits of the canvas' axes
        '''
        return self.plot_canvas.axes.get_ylim()
        
class DoubleCanvasPlot(FigureCanvas):
    '''Drawing Canvas for potting two axis on the same figure. Can be used to plot
    SAS and WAS data at the same time or the scattering curve plus a plot of the fitting
    (e.g. a SAS curve and the distribution of sizes of the fit)
    '''
    def __init__(self, parent=None, figsize = (5,4), rows = 1, cols = 2, rowSpan = 1, colSpan = 1, **kwargs):
                 #firstAx = [[0,1],[0,1]], secondAx = [[0,1],[1,2]] ):
        '''Initiates the canvas.
            Args:
                parent (QtWidget): the parent widget to which the canvas is associated
                    Defauts to None
                figsize (list of int): the size of the figure in the form: (width,height).
                    This will be used ot create the figure. Defaults to (5,4)
                rows (int): the number of rows in which the two axis should be disposed
                    Defaults to 1
                cols (int): the number of columns over which the axis should be disposed.
                    Defaults to 2.
                rowSpan (int): number of rows spanned by the first axis. The remaining
                    rows are attributed to the second axis. Defaults to 1.
                colSpan (int): number of cols spanned by the first axis. The remaining
                    cols are attributed to the second axis. Defaults to 1
                
        '''
        if kwargs.get('verbose', False):
			self.logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
        else:
			self.logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
			
            
        firstAx = [[0,rowSpan],[0,colSpan]]
        if rowSpan == rows:
            secondAx = [[0,rowSpan]]
        else:
            secondAx = [[rowSpan,rows]]
        if colSpan == cols:
            secondAx.append([0,colSpan])
        else:
            secondAx.append([colSpan,cols])
        self.fig = Figure(figsize=figsize)
        gs = mpl.gridspec.GridSpec(rows,cols)
        self.axes1 = self.fig.add_subplot(gs[slice(*firstAx[0]),slice(*firstAx[1])])
        self.axes2 = self.fig.add_subplot(gs[slice(*secondAx[0]),slice(*secondAx[1])])
        
        super(DoubleCanvasPlot, self).__init__(self.fig)
        #FigureCanvas.__init__(self, fig)
        self.setParent(parent)

        self.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                QtWidgets.QSizePolicy.Expanding)
        self.updateGeometry()
        
class DoublePlot(QtWidgets.QWidget):
    '''Wrapper class to place the DoubleCanvasPlot in a widget
    '''
    def __init__(self, parent=None, **kwargs):
        if kwargs.get('verbose', False):
			self.logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
        else:
			self.logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
			
        super(DoublePlot, self).__init__(parent)
        self.scatteringObject = kwargs.get('plotData', None)
        
        self.plot_layout = QtWidgets.QVBoxLayout(self)
        self.plot_canvas = DoubleCanvasPlot(self, figsize = (5,4))

        self.navi_toolbar = NavigationToolbar(self.plot_canvas, self)
        self.plot_layout.addWidget(self.plot_canvas)  # the matplotlib canvas
        self.plot_layout.addWidget(self.navi_toolbar)
        self.setLayout(self.plot_layout)
        if self.scatteringObject is not None and fitType is not None:
            self.plot_fit()
        elif self.scatteringObject is not None:
            self.plot_data()
    
    def cla(self):
        '''Clears bothof the axis in the figure
        '''
        self.plot_canvas.axes1.cla()
        self.plot_canvas.axes2.cla()
    
    def remove_line(self, line):
        '''searches both axes to see if the given line is in either.
        If it is it is removed
        '''
        if isinstance(line, lines.Line2D):
            if line in self.plot_canvas.axes1.lines:
                self.plot_canvas.axes1.lines.remove(line)
            elif line in self.plot_canvas.axes2.lines:
                self.plot_canvas.axes2.lines.remove(line)
            else:
                self.logging.error('{} was not found in either axis'.format(line))
 
    def plot_data(self):
        '''Uses the plotting function of the scattering object to plot the data on the two
        available axis
        '''
        self.scatteringObject.plot_data(axs = [self.plot_canvas.axes1,self.plot_canvas.axes2])
    
    def plot_fit(self):
        '''Uses the fitting plot function of the scattering object to plot the data on the
        two available axis
        '''
        self.scatteringObject.plot_fit()
        
    def redraw(self):
        '''Redraws the canvas after lines have been added or removed from the figure
        '''
        self.plot_canvas.draw()
        
class MultiplePlot(QtWidgets.QWidget):
    '''Widget used to plot a sequence of scattering objects. It is composed of a plotting widget,
    two buttons to move the current position forward and backwards by one, and a slider to select
    any avalable fitting 
    '''
    def __init__(self, parent=None, **kwargs):
        if kwargs.get('verbose', False):
			self.logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
        else:
			self.logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
			
        super(MultiplePlot, self).__init__(parent)
        
        #Setup the data for the plotting
        self.plotData = kwargs.get('plotData', None)
        if self.plotData is not None:
            self.scatteringObject = self.plotData['Object']

        self.plotType = kwargs.get('plot_type', 'SAS')
        
        self.currPlot = 0
        
        if self.plotType in ('SAS','WAS'):
            self.plotWidget = SinglePlot(self)
        else:
            self.plotWidget = DoublePlot(self)

        
        
        #Create the widget in which the data is going to be plotted
        self.layout = QtWidgets.QVBoxLayout(self)
        self.layout.setContentsMargins(0,0,0,0)
        self.layout.setSpacing(0)
        
        
        self.scroll_toolBar = QtWidgets.QToolBar(self)
        self.scrollBar = QtWidgets.QScrollBar(QtCore.Qt.Horizontal,self.scroll_toolBar)
        self.scrollBar.setSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding,QtWidgets.QSizePolicy.MinimumExpanding)
        self.scrollBar.sliderMoved.connect(self.sliderMoving)
        self.scrollBar.sliderReleased.connect(self.sliderChanged)
        
        self.text_toolBar = QtWidgets.QToolBar(self)
        self.text_toolBar_layout = QtWidgets.QHBoxLayout(self.text_toolBar)
        
        spacerL = QtWidgets.QWidget()
        spacerR = QtWidgets.QWidget()
        spacerL.setSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding,QtWidgets.QSizePolicy.MinimumExpanding)
        spacerR.setSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding,QtWidgets.QSizePolicy.MinimumExpanding)
        self.prev_button = QtWidgets.QPushButton('<<', self.text_toolBar)
        self.prev_button.clicked.connect(self.prevPlot)
        self.next_button = QtWidgets.QPushButton('>>', self.text_toolBar)
        self.next_button.clicked.connect(self.nextPlot)
        self.lineEdit = QtWidgets.QLabel(self.text_toolBar)
        self.lineEdit.setText('LineEditText')
        

        self.scroll_toolBar.addWidget(self.scrollBar)

        self.text_toolBar.addWidget(spacerL)
        self.text_toolBar.addWidget(self.prev_button)
        self.text_toolBar.addWidget(self.lineEdit)
        self.text_toolBar.addWidget(self.next_button)
        self.text_toolBar.addWidget(spacerR)

        self.text_toolBar_layout.setAlignment(self.lineEdit, QtCore.Qt.AlignHCenter)
        
        self.layout.addWidget(self.scroll_toolBar)
        self.layout.addWidget(self.text_toolBar)
        self.layout.addWidget(self.plotWidget)
        
       
        self.numb_curves = 0
        if self.scatteringObject is not None:
            self.InitializeValues()
        
    def InitializeValues(self):
        '''InitializeValues sets all the variables needed to move between the available
        data sets and visialize them
        '''
        #self.numb_curves = self.scatteringObject.size
        if self.plotType == 'SAS':
            self.selectedPlots = self.plotData['SAS']
            self.numbCurves = len(self.selectedPlots)
        elif self.plotType == 'WAS':
            self.selectedPlots = self.plotData['WAS']
            self.numbCurves = len(self.selectedPlots)
        else:
            self.selectedPlots = self.plotData['SAS']
            self.numbCurves = len(self.selectedPlots)
            
        self.typeObject = self.scatteringObject.avlbCurves
        self.lineEdit.setText('1/{}'.format(self.numbCurves+1))
        self.scrollBar.setMinimum(1)
        self.scrollBar.setMaximum(self.numbCurves)
        self.scrollBar.setValue(1)
        if not isinstance(self.scatteringObject,SWAS_Sequence):
            self.text_toolBar.hide()
            self.scroll_toolBar.hide()
        if self.fitType is None:
            self.plot_data()
        else:
            self.plot_fit()
    
    def set_data(self, scatteringObj, selectedPlots, plotType):
        '''set_data is a quick setter function to set the scattering object,
        the data to plot and the type of plot and fit
        '''
        self.scatteringObject = scatteringObj
        self.selectedPlots = selectedPlots
        self.plotType = plotType
        
    
    def plot_data(self, **kwargs):
        '''plot_data used the plotting functions of the object to plot
        the data on the available axis after clearing them
        '''
        
        #print self.scatteringObject[self.currPlot].SAS.q
        self.plotWidget.cla()
        if self.plotType == 'SAS':
            if isinstance(self.scatteringObject,SmallAngleScattering):
                self.scatteringObject.plot_data(ax = self.plotWidget.plot_canvas.axes, **kwargs)
            if isinstance(self.scatteringObject,SWAS):
                self.scatteringObject.SAS.plot_data(ax = self.plotWidget.plot_canvas.axes, **kwargs)
            if isinstance(self.scatteringObject,SWAS_Sequence):
                self.scatteringObject[self.selectedPlots[self.currPlot]].SAS.plot_data(ax = self.plotWidget.plot_canvas.axes, **kwargs)
            else:
                self.logging.error('Cannot plot Small angles for the selected scattering data')
        if self.plotType == 'WAS':
            if isinstance(self.scatteringObject,WideAngleScattering):
                self.scatteringObject.plot_data(ax = self.plotWidget.plot_canvas.axes, **kwargs)
            if isinstance(self.scatteringObject,SWAS):
                self.scatteringObject.WAS.plot_data(ax = self.plotWidget.plot_canvas.axes, **kwargs)
            if isinstance(self.scatteringObject,SWAS_Sequence):
                self.scatteringObject[self.selectedPlots[self.currPlot]].WAS.plot_data(ax = self.plotWidget.plot_canvas.axes, **kwargs)
            else:
                self.logging('Cannot plot wide angles for the selected scattering data')
        
        if self.plotType == 'SWAS_Sequence':
            #print 'plotting ', self.currPlot, 'which is object ',self.selectedPlots[self.currPlot]
            #print 'of ',self.scatteringObject, ' ', self.scatteringObject[self.selectedPlots[self.currPlot]].SAS.q
            if isinstance(self.scatteringObject,SWAS_Sequence):
                self.scatteringObject[self.selectedPlots[self.currPlot]].plot_data(axs = [self.plotWidget.plot_canvas.axes1, self.plotWidget.plot_canvas.axes2],\
                                                                                       fig = self.plotWidget.plot_canvas.fig,**kwargs)

        self.plotWidget.redraw()
                
    def nextPlot(self):
        '''Sets the current plot to one after the current one. If the last plot is
        currently being shown it does nothing. 
        '''
        if (self.currPlot+1)<self.numbCurves:
            self.currPlot += 1
            self.scrollBar.setValue(self.currPlot+1)
            self.lineEdit.setText('{}/{}'.format(self.currPlot+1, self.numbCurves))
            self.plot_data()

            
    def prevPlot(self):
        '''Sets the current plot to one before the current one. If the first plot is
        currently being shown it does nothing. 
        '''
        if (self.currPlot-1) >= 0:
            self.currPlot -= 1
            self.scrollBar.setValue(self.currPlot+1)
            self.lineEdit.setText('{}/{}'.format(self.currPlot+1, self.numbCurves))
            self.plot_data()
    
    def sliderMoving(self):
        '''Updates in real time the number being shown on the text as the slider is moved
        '''
        self.lineEdit.setText('{}/{}'.format(self.scrollBar.value(), self.numbCurves))
    
    def sliderChanged(self):
        '''Updates the current plot based on the selection done with the slider
        '''
        self.currPlot = self.scrollBar.value()-1
        self.plot_data()
        
    def ax_x_lim(self):
        '''Returns the x limits of the current axes/axis
        '''
        if isinstance(self.plotWidget, DoublePlot):
            return [self.plotWidget.plot_canvas.axes1.get_xlim(),self.plotWidget.plot_canvas.axes1.get_xlim()]
        else:
            return self.plotWidget.plot_canvas.axes.get_xlim()
    
    def axvline(self,x):
        if isinstance(self.plotWidget, SinglePlot):
            self.plotWidget.axvline(x)
    
    def remove_line(self, line):
        if isinstance(line,lines.Line2D):
            self.plotWidget.remove_line(line)
        else:
            self.logging.error('{} is not a matplotlib 2D line'.format(line))
            
class FitPlot(QtWidgets.QWidget):
    
    def __init__(self, parent=None, **kwargs):
        super(FitPlot, self).__init__(parent)
        self.scatteringObject = kwargs.get('plotData', None)
        self.fitType = kwargs.get('fitType',None)
        self.plot_layout = QtWidgets.QVBoxLayout(self)
        self.plot_canvas = DoubleCanvasPlot(self, figsize = (5,4))

        self.navi_toolbar = NavigationToolbar(self.plot_canvas, self)
        self.plot_layout.addWidget(self.plot_canvas)  # the matplotlib canvas
        self.plot_layout.addWidget(self.navi_toolbar)
        self.setLayout(self.plot_layout)
        if self.scatteringObject is not None and self.fitType is not None:
            self.plot_fit()
        elif self.scatteringObject is not None:
            self.plot_data()
    
    def cla(self):
        self.plot_canvas.axes1.cla()
        self.plot_canvas.axes2.cla()
    
    def remove_line(self, line):
        if isinstance(line, lines.Line2D):
            if line in self.plot_canvas.axes1.lines:
                self.plot_canvas.axes1.lines.remove(line)
            elif line in self.plot_canvas.axes2.lines:
                self.plot_canvas.axes2.lines.remove(line)
            else:
                print line, ' is not in either axes'
    
    def plot_data(self):
        self.scatteringObject.plot_fit(self.fitType, axs = [self.plot_canvas.axes1,self.plot_canvas.axes2])
    
    def plot_fit(self):
        self.scatteringObject.plot_fit(self.fitType, axs = [self.plot_canvas.axes1,self.plot_canvas.axes2])
    def redraw(self):
        self.plot_canvas.draw()
        
class MultipleFitPlot(MultiplePlot):
    
    def __init__(self,  parent=None, **kwargs):
        super(MultipleFitPlot, self).__init__(parent, **kwargs)
        self.fitType = kwargs.get('fitType')
        self.plotWidget = DoublePlot(self)
        
    def set_data(self, scatteringObj, selectedPlots, plotType, fitType):
        '''set_data is a quick setter function to set the scattering object,
        the data to plot, the type of plot and fit
        '''
        self.scatteringObject = scatteringObj
        self.selectedPlots = selectedPlots
        self.plotType = plotType
        self.fitType = fitType
        
    
    def plot_data(self, **kwargs):
        '''plot_data used the plotting functions of the object to plot
        the data on the available axis after clearing them
        '''
        
        #print self.scatteringObject[self.currPlot].SAS.q
        self.plotWidget.cla()
        if self.plotType == 'SAS':
            if isinstance(self.scatteringObject,SmallAngleScattering):
                self.scatteringObject.plot_fit(self.fitType, ax = [self.plotWidget.plot_canvas.axes1, self.plotWidget.plot_canvas.axes2], **kwargs)
            if isinstance(self.scatteringObject,SWAS):
                self.scatteringObject.SAS.plot_fit(self.fitType, ax = [self.plotWidget.plot_canvas.axes1, self.plotWidget.plot_canvas.axes2], **kwargs)
            if isinstance(self.scatteringObject,SWAS_Sequence):
                self.scatteringObject[self.selectedPlots[self.currPlot]].SAS.plot_fit(self.fitType, ax = [self.plotWidget.plot_canvas.axes1, self.plotWidget.plot_canvas.axes2],\
                                                                                       **kwargs)
            else:
                self.logging.error('Cannot plot Small angles for the selected scattering data')
        if self.plotType == 'WAS':
            if isinstance(self.scatteringObject,WideAngleScattering):
                self.scatteringObject.plot_fit(self.fitType, ax = self.plotWidget.plot_canvas.axes, **kwargs)
            if isinstance(self.scatteringObject,SWAS):
                self.scatteringObject.WAS.plot_fit(self.fitType, ax = [self.plotWidget.plot_canvas.axes1, self.plotWidget.plot_canvas.axes2], **kwargs)
            if isinstance(self.scatteringObject,SWAS_Sequence):
                self.scatteringObject[self.selectedPlots[self.currPlot]].WAS.plot_fit(self.fitType, ax = [self.plotWidget.plot_canvas.axes1, self.plotWidget.plot_canvas.axes2], **kwargs)
            else:
                self.logging('Cannot plot wide angles for the selected scattering data')
        
        if self.plotType == 'SWAS_Sequence':
            #print 'plotting ', self.currPlot, 'which is object ',self.selectedPlots[self.currPlot]
            #print 'of ',self.scatteringObject, ' ', self.scatteringObject[self.selectedPlots[self.currPlot]].SAS.q
            if isinstance(self.scatteringObject,SWAS_Sequence):
                self.scatteringObject[self.selectedPlots[self.currPlot]].plot_fit(self.fitType, axs = [self.plotWidget.plot_canvas.axes1, self.plotWidget.plot_canvas.axes2],\
                                                                                       fig = self.plotWidget.plot_canvas.fig,**kwargs)

        self.plotWidget.redraw()
        