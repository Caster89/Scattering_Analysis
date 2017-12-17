import unittest
import os, sys
import numpy as np
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(CURRENT_DIR))

from Scattering_Object import ScatteringObject
from SAS import SmallAngleScattering
import Fitting_Models

class TestSAS(unittest.TestCase):
    def setUP(self):
        pass

    def test_loadFromFile(self):
        fileName = r'../Data/NC_B_91.dat'
        data = np.genfromtxt(fileName, skip_header=2)
        sasData = SmallAngleScattering(fileName, qUnits='nm', IUnits='m',
            DetDist=5, Temp = 25, skip_header=0)
        self.assertTrue((sasData.q==data[:,0]).all())
        self.assertTrue((sasData.Iraw == data[:,1]).all())
        self.assertTrue((sasData.Inorm == data[:,1]).all())
        self.assertTrue((sasData.Iabs == data[:,1]).all())

        self.assertEqual(sasData.qUnits,'nm')
        self.assertEqual(sasData.IUnits,'m')

    def test_loadFromData(self):
        fileName = r'../Data/NC_B_91.dat'
        data = np.genfromtxt(fileName)
        sasData = SmallAngleScattering()
        sasData.create_from_data(q=data[:,0], I=data[:,1], qUnits='nm', IUnits='m',
            DetDist=5, Temp = 25)
        self.assertTrue((sasData.q==data[:,0]).all())
        self.assertTrue((sasData.Iraw == data[:,1]).all())
        self.assertTrue((sasData.Inorm == data[:,1]).all())
        self.assertTrue((sasData.Iabs == data[:,1]).all())

        self.assertEqual(sasData.qUnits,'nm')
        self.assertEqual(sasData.IUnits,'m')

    def test_plotData(self):
        q= np.logspace(-2,1,1000)
        R_av=10.0
        sigma=1.0
        I0=1.0
        I = Fitting_Models.single_gauss_spheres(q,R_av = R_av,sigma = sigma,I0 = I0,bckg=0)
        Inoise = I*np.random.normal(0,0.001,I.size)
        sasData = SmallAngleScattering()
        sasData.create_from_data(q=q, I=Inoise, qUnits='nm', IUnits='m',
            DetDist=5, Temp = 25)
        sasData.plot_data(verbose=True)
"""
    def test_SingleGauss(self):
        q= np.logspace(-2,1,1000)
        R_av=10.0
        sigma=3.0
        I0=1.0
        bckg = 0
        I = Fitting_Models.single_gauss_spheres(q,R_av = R_av,sigma = sigma,I0 = I0,bckg=0)
        Inoise = I*np.random.normal(0,0.01,I.size)
        Inoise[Inoise<=0.]=1e-10
        sasData = SmallAngleScattering()
        sasData.create_from_data(q=q, I=Inoise, qUnits='nm', IUnits='m',
            DetDist=5, Temp = 25)
        paramSugg = {'R_av':{'value': 10, 'min':9, 'max':11},
                     'sigma':{'value': 3.},
                     'I0':{'value': 1., 'max': 5, 'min': 0},}
        sasData.fit_data('Sing_Gauss', paramSugg=paramSugg, verbose = True)
        result = sasData.fitResults['Sing_Gauss']
        self.assertAlmostEqual(R_av, result['R_av'], places=3)
        self.assertAlmostEqual(sigma, result['sigma'], places=3)
        self.assertAlmostEqual(I0, result['I0'], places=3)
        self.assertAlmostEqual(bckg, result['bckg'], places=3)
"""

if __name__ == '__main__':
    unittest.main()
