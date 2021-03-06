import numpy as np


_AvlbUnits = ['A','nm','um','mm','cm','m']

_UnitsSymbols = {'A': '$\AA$','nm':'nm','um':'$\mu$m','mm':'mm','cm':'cm','m':'m'}

_UnitsConv = {'A':1e-10,'nm':1e-9,'um':1e-6,'mm':1e-3,'cm':1e-2,'m':1e0}

_AvlbSASFit = ['Sing_Gauss',
               'Double_Gauss',
               'Sing_Schultz',
               'Double_Schultz',
               'Mono_Cube',
               'Gauss_Cube',
               'EM',
               'Dist_from_EM',]

_AvlbSASFitDic = {'Single Gaussian Distribution': 'Sing_Gauss',
                  'Double Gaussian Distribution': 'Double_Gauss',
                  'Single Schultz Distribution': 'Sing_Schultz',
                  'Double Schultz Distribution': 'Double_Schultz',
                  'Monodisperse Cubes': 'Mono_Cube',
                  'Single Gaussian Cubes': 'Gauss_Cube',
                  'Expected Maximization': 'EM',
                  'Distribution From Expected Maximization': 'Dist_from_EM'}

_AvlbSASFitDicInv = {'Sing_Gauss': 'Single Gaussian Distribution',
                     'Double_Gauss': 'Double Gaussian Distribution',
                     'Sing_Schultz': 'Single Schultz Distribution',
                     'Double_Schultz': 'Double Schultz Distribution',
                     'Mono_Cube': 'Monodisperse Cube',
                     'Gauss Cube': 'Single Gaussian Cubes',
                     'EM': 'Expected Maximization',
                     'Dist_from_EM': 'Distribution From Expected Maximization'}

_AvlbWASFitDic = {'Pseudo-Voigt Model': 'PVM'}

_AvlbWASFit = ['PVM']

_lmfitModels = ['Sing_Gauss',
                'Double_Gauss',
                'Sing_Schultz',
                'Double_Schultz',
                'Mono_Cube',
                'Gauss_Cube']

_lmfitModelsMapping = {'Sing_Gauss': 'SingleGaussianSphereModel',
                        'Double_Gauss': 'DoubleGaussianSphereModel',
                        'Sing_Schultz': 'SingleSchultzSphereModel',
                        'Double_Schultz': 'DoubleSchultzSphereModel',
                        'Mono_Cube': 'MonodisperseCubeModel',
                        'Gauss Cube': 'SingleGaussianCubeModel',}


_lmfitModelFunctions = {'Sing_Gauss': 'single_gauss_spheres',
                        'Double_Gauss': 'double_gauss_spheres',
                        'Sing_Schultz': 'single_schultz_spheres',
                        'Double_Schultz': 'double_schultz_spheres',
                        'Mono_Cube': 'monodisperse_cube',
                        'Gauss Cube': 'single_gaussian_cube',}

_lmfitModelParams = {'Sing_Gauss': {'R_av': {'value': 1, 'min': 0, 'max': np.inf, 'vary': True},
                                    'sigma': {'value': 1, 'min': 0, 'max': np.inf, 'vary': True},
                                    'I0': {'value': 1, 'min': -np.inf, 'max': np.inf, 'vary': True},
                                    'bckg': {'value': 0, 'min': -np.inf, 'max': np.inf, 'vary': False}, },
                     'Double_Gauss': {'R1_av': {'value': 1, 'min': 0, 'max': np.inf, 'vary': True},
                                      'sigma1': {'value': 1, 'min': 0, 'max': np.inf, 'vary': True},
                                      'R2_av': {'value': 1, 'min': 0, 'max': np.inf, 'vary': True},
                                      'sigma2': {'value': 1, 'min': 0, 'max': np.inf, 'vary': True},
                                      'I0': {'value': 1, 'min': -np.inf, 'max': np.inf, 'vary': True},
                                      'ratio': {'value': 0.5,'min': 0, 'max':1., 'vary': True},
                                      'bckg': {'value': 0, 'min': -np.inf, 'max': np.inf, 'vary': False}, },
                     'Sing_Schultz': {'R_av': {'value': 1, 'min': 0, 'max': np.inf, 'vary': True},
                                      'Z': {'value': 1, 'min': 0, 'max': np.inf, 'vary': True},
                                      'I0': {'value': 1, 'min': -np.inf, 'max': np.inf, 'vary': True},
                                      'bckg': {'value': 0, 'min': -np.inf, 'max': np.inf, 'vary': False}, },
                     'Double_Schultz': {'R1_av': {'value': 1, 'min': 0, 'max': np.inf, 'vary': True},
                                        'Z1': {'value': 1, 'min': 0, 'max': np.inf, 'vary': True},
                                        'R2_av': {'value': 1, 'min': 0, 'max': np.inf, 'vary': True},
                                        'Z2': {'value': 1, 'min': 0, 'max': np.inf, 'vary': True},
                                        'I0': {'value': 1, 'min': -np.inf, 'max': np.inf, 'vary': True},
                                        'ratio': {'value': 0.5,'min': 0, 'max':1., 'vary': True},
                                        'bckg': {'value': 0, 'min': -np.inf, 'max': np.inf, 'vary': False}, },
                     'Mono_Cube': {'L': {'value': 1, 'min': 0, 'max': np.inf, 'vary': True},
                                   'I0': {'value': 1, 'min': -np.inf, 'max': np.inf, 'vary': True},
                                   'bckg': {'value': 0, 'min': -np.inf, 'max': np.inf, 'vary': False}, },
                     'Gauss Cube': {'L': {'value': 1, 'min': 0, 'max': np.inf, 'vary': True},
                                    'sigma':{'value': 1, 'min': 0, 'max': np.inf, 'vary': True},
                                    'I0': {'value': 1, 'min': -np.inf, 'max': np.inf, 'vary': True},
                                    'bckg': {'value': 0, 'min': -np.inf, 'max': np.inf, 'vary': False}, }
                     }

_lmfitParamsOrder = ['R_av', 'R1_av', 'L', 'sigma', 'sigma1', 'Z', 'Z1', 'R2_av', 'sigma2', 'Z2', 'I0', 'ratio', 'bckg']

_lmfitDistFunctions = {'Sing_Gauss': 'single_gauss_distribution',
                       'Double_Gauss': 'double_gauss_distribution',
                       'Sing_Schultz': 'single_schultz_distribution',
                       'Double_Schultz': 'double_schultz_distribution'}
