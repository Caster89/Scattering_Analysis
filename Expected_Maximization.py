from config import _AvlbUnits, _UnitsConv, _AvlbSASFit, _AvlbWASFit, _lmfitModels, _lmfitModelFunctions, _lmfitDistFunctions
import numpy as np
import numpy.matlib
import logging

def sphere_form_factor(q, R):
   """sphere_form_factor: Calculates the form factor for a sphere of given
   radius over the q-values specified.
      Args:
         q (numpy.array): vector containing the q values for which the form
            factor has to be calculated
         R (numpy.array): the radii of the spheres
      Returns:
         The form factor array with the same size as qxR
   """
   RR, qq = np.meshgrid(R,q)
   F = (3*(np.sin(qq*RR)-qq*RR*np.cos(qq*RR))/(qq*RR)**3)**2
   return F

def expected_maximization_method(q, I, Rmin, Rmax=0, Rvec = [], eps = 1e-12, k=3, numbElements=100, maxIter = 10000, verbose = False):
   """espected_maximisation_method: method for finding the ditribution of particles
   given the scattering curve. The method is taken from:
   "A Robust Inversion Method According to A New Notion of Regularization For Poisson
   Data With an Application to Nanoparticle volume Determination"
   Benvenuto, F. Haddar, H., Lantz, B.
      Args:
         q (numpy.array): the vector containing the values of the wave-vector
         I (numpy.array): is the intensity of the signal from the scattering
            after the radial averaging. Has to have th same size as q
         Rmin (int): the minimum size of the particle
         Rmax (int): the maximum size of the particle. If Rmax is set to 0 then
            its value is calculated from the resolution of the q vector. Defaults
            to 0
         eps (int): one of the two stopping criteria. Defaults to 1e-12
         k (int): second parameter used for hte stopping criteria. Defaults to 3
         numbElements (int): How many bins ot divide the radii range in.
            Defaults to 100
         maxIter (int): maximum number of iterations before stopping.
            Defaults to 10000
         verobse (bool): determins whether the algorithm outputs messages during
            the iterations. Defaults to False
      Returns:
         xk1 (numpy.array): the probability distribution array
         Rvec (numpy.array): the vector with the radii
         H (numpy.array): the matrix witht he form factors, needed to plot the
            scattering curve.
   """
   if verbose:
		logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
   else:
		logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
   ## Comments with a double '#' refer to parts of the code which can be used
   ## to insert an extra radius with value 0. This radius has a constant form
   ## factor of 1, and might be used to add/remove a constant background which
   ## was not correctly subtracted.

   #I has to be a vertical vector
   I.shape=(len(I),1)


   #The maximum radius detectable is defined by 1/(4dq) where dq is the step of the availabel q vectors
   #in order to get the biggest radius possible, the smallest delta q is selected
   if Rmax==0:
      logging.debug('Setting Rmax based on dq')
      deltaQ = q[1:]-q[:-1]
      Rmax = 1/(4*min(deltaQ))
   #Creation of the radius vector, setting the step based on
   #the number of elements desired
   if not Rvec:
      Rvec=np.arange(Rmin,Rmax,float(Rmax-Rmin)/numbElements)
   else:
      Rvec = Rvec
   ##Rvec2 = np.insert(Rvec,0,0)

   #Creation of the Volume vector from the radius vector
   Vvec = Rvec**3*4/3*np.pi
   ##Vvec2 = Rvec2**3*4/3*np.pi
   ##Vvec2[0] = 1
   ##Vvec2 = np.insert(Vvec,0,1)

   #the form factor matrix is created. q increases descending the column, while the radius increases
   #moving left along the rows
   FF = sphere_form_factor(q,Rvec)

   ##bckVec = np.ones((FF.shape[0],1))
   ##FF2 = np.hstack((bckVec,FF))

   #The base H matrix is the product of the elements of the form factor matrix multiplied by
   #the respective volume
   H = FF*np.matlib.repmat(Vvec,len(q),1)
   ##H2= FF2*np.matlib.repmat(Vvec2,len(q),1)

   #The basic components of the algorithm are created
   HT = H.transpose()
   ##HT2 = H2.transpose()
   one = np.ones((HT.shape[1],1))
   ##one2 = np.ones((HT2.shape[1],1))
   Hvec = np.dot(HT,one)
   ##Hvec2 = np.dot(HT2,one2)
   xk = Vvec
   ##xk2 = Vvec2
   xk.shape=(len(xk),1)
   ##xk2.shape=(len(xk2),1)

   cont = True
   i = 0
   while cont:
      #The k+1 step is obtained from the k step following the formula:
      #x_k+1 = x_k/(H^T*[1])*H^T*I/(H*x_k)

      Hxk=np.dot(H,xk)

      #starting from the right
      xk1 = I/Hxk
      xk1=np.dot(HT,xk1)
      xk1 = xk1/Hvec
      xk1=xk*xk1

      #the stopping criteria depends on two functions A and B
      #A = (H^T*(I/(H*x_k)-1))^2
      #B = (H^T)^2*[1]/(H*x_k)
      A = I/np.dot(H,xk1)
      A = A-np.ones(A.shape)
      A = np.dot(HT,A)
      A=A**2

      B = np.dot(H,xk1)
      B= np.ones(B.shape)/B
      B=np.dot(HT**2,B)


      #set xk equal to the new vector
      xk = xk1
      A = A[xk1>eps]
      B = B[xk1>eps]
      #If A > k^2*B element by element and all the values of xk > eps
      #or if 1000000 iterations have been done the algorithm is finished
      if (np.all(A<=k**2*B) and A.size is not 0) or i>maxIter:
         if verbose and i>maxIter:
            logging.error('Max number of iterations reached')
         cont =False
      i+=1

   return [xk1, Rvec, H]

def plot_expected_maximiztion(q, I, Rvec, xk, H, volumeDist = True, ax = None, **kwargs):
   """plot_expected_maximization: method used to plot the results of the
   expected maximization algorithm. It plots the data on the left axis and
   the distribution on the right. If the axes are provided it used the
   ones given, if not it creates new ones.
      Args:
         q (numpy.array): vector containing the q values
         I (numpy.array): vector containing the experimental scattering curve
         R_vec (numpy.array): vector containing the radii used by the algorithm
         xk (numpy.array): the probability density vector calculated by the EM method
         H (numpy.array): the matrix containing hte scattering curve for the different
            form factor
         volumeDist (bool): determines whether the distribution is plotted as a volume
            distribution. If set to false then a number distribution is used. Defaults
            to True
         ax (list of plt.axis): The axes which can be used to plot the data. If None
            then the axes are created. Defaults to None
         **kwargs (dict): used to pass preferences regarding the plots
   """
   qConvFact = _UnitsConv[kwargs.get('qunits','nm')]/_UnitsConv[kwargs.get('qorigUnits','nm')]
   IConvFact = _UnitsConv[kwargs.get('Iunits','m')]/_UnitsConv[kwargs.get('IOrigUnits','m')]
   if ax is None:
      ax = []
      fig = plt.figure(figsize=(22,10))
      fig.subplots_adjust(hspace=0.35,wspace = 0.35)
      ax[0] = plt.subplot2grid((1,3),(0,0),rowspan=1,colspan=2)
      ax[1] = plt.subplot2grid((1,3),(0,2),rowspan=1,colspan=1)

   ax[0].loglog(q*qConvFact, I*IConvFact, color = 'k', linestyle = 'None', marker = 'o',\
                markeredgecolor='None', label = 'Experimental Points')

   ax[0].loglog(q*qConvFact, np.dot(H,xk)*IConvFact, color = kwargs.get('color','r'),\
                linewidth = kwargs.get('linewidth',2), label = kwargs.get('label', 'EM Fit'))
   ax[0].legend(loc = 'best', frameon = False, fontsize = 20)
   ax[0].tick_params(size = 6,width = 3,labelsize = 22)
   ax[0].tick_params(which = 'minor', size = 3,width = 2,labelsize = 22)
   ax[0].set_xlabel('Wavevector q ({}$\mathrm{{^{{-1}} }}$)'.format(kwargs.get('qunits','nm')),\
                    fontsize = 22)
   ax[0].set_ylabel('Intensity ({}$\mathrm{{^{{-1}} }}$)'.format(kwargs.get('Iunits','m')),fontsize = 22)

   ax[0].set_xlim(min(q*qConvFact)*0.95,max(q*qConvFact)*1.05)
   ax[0].set_ylim(min(I*IConvFact)*0.9,max(I*IConvFact)*1.1)
   #totVolume = np.sum(xk)
   #Volume = 4./3.*Rvec**3*np.pi
   #numbDist = np.array([x/v for x,v in zip(xk,Volume)])
   #volDist = xk/totVolume
   Volume = np.array(4./3.*Rvec**3*np.pi)
   totVolume = np.sum(Volume*xk[:,0])
   numbDist = np.array([x/v for x,v in zip(xk,Volume)])
   #volDist = xk/totVolume
   #print Volume.shape
   #print xk.shape
   volDist = (Volume*xk[:,0])/totVolume
   if volumeDist:
       ax[1].bar(Rvec/qConvFact,volDist,width = (Rvec[1]-Rvec[0])/qConvFact, linewidth = 0.5)
       ax[1].set_ylabel('Volume Distribution',fontsize = 22)
   else:
       ax[1].bar(Rvec*qConvFact, numbDist, width = (Rvec[1]-Rvec[0])/qConvFact, linewidth = 0.5)
       ax[1].set_ylabel('Number Distribution',fontsize = 22)
       #axDist.set_ylim(0,max(numbDist[10:]))
   ax[1].tick_params(size = 6,width = 3,labelsize = 22)
   ax[1].set_xlabel('Radius ({})'.format(kwargs.get('units','nm')),fontsize = 22)





   return ax
