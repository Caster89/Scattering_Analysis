_AvlbUnits = ['A','nm','mm','cm','m']
_UnitsConv = {'A':1e-10,'nm':1e-9,'mm':1e-3,'cm':1e-2,'m':1e0}
_AvlbSASFit = ['Sing_Gauss', 'Double_Gauss', 'Sing_Schultz', 'Double_Schultz', 'EM']
_AvlbSASFitDic = {'Single Gaussian Distribution' : 'Sing_Gauss',\
                   'Double Gaussian Distribution' : 'Double_Gauss',\
                   'Single Schultz Distribution' : 'Sing_Schultz',\
                   'Double Schultz Distribution' : 'Double_Schultz',\
                   'Expected Minimization' : 'EM'}
_AvlbWASFitDic = {'Pseudo-Voigt Model': 'PVM'}
_AvlbWASFit = ['PVM']
_lmfitModels = ['Sing_Gauss', 'Double_Gauss', 'Sing_Schultz', 'Double_Schultz']
_lmfitModelFunctions = {'Sing_Gauss': 'single_gauss_spheres', 'Double_Gauss': 'double_gauss_spheres',\
						'Sing_Schultz': 'single_schultz_spheres', 'Double_Schultz': 'double_schultz_spheres'}
_lmfitDistFunctions = {'Sing_Gauss': 'single_gauss_distribution', 'Double_Gauss': 'double_gauss_distribution',\
                       'Sing_Schultz': 'single_schultz_distribution', 'Double_Schultz': 'double_schultz_distribution'}