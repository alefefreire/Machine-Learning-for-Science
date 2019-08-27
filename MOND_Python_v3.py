from __future__ import print_function
import numpy as np
import pandas as pd
import csv
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy.interpolate import interp1d
import scipy.optimize as op
import scipy
from scipy import *
from scipy.special import expi
#import lmfit
#from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit
import emcee
#from emcee import PTSampler
from matplotlib import rcParams
import time
from scipy.interpolate import InterpolatedUnivariateSpline

############################################################################################ MCMC PARAMETERS #######################################################################################################################

imin = 1
imax = 90
Dmin = 0.5
LYd = np.log10(0.5)
LYb = np.log10(0.7)
sigmaY = 0.1
my_moves = emcee.moves.StretchMove ()
alpha = 5000
tau = 100
nsteps = alpha*tau 
burn_in = alpha*tau/10
thin = 1




#############################################################################################  IMPORTING AND DEFINING FUNCTIONS ####################################################################################################
fileCamB=open('/home/alefe/Cluster/SPARC/CamB-Min-MONDRAR-var-Fa0-GYDi-5.txt','r')
file_readCamB=pd.read_table(fileCamB,header=None)
file_arrayCamB=np.array(file_readCamB)
chi2_stringCamB=file_arrayCamB[1,0]
best_fitCamB=file_arrayCamB[4,0]
best_fit_listCamB=best_fitCamB.split()
best_fit_valuesCamB=[eval(best_fit_listCamB[0]),eval(best_fit_listCamB[1]),eval(best_fit_listCamB[2]),eval(best_fit_listCamB[3])]
reference_CamB=file_arrayCamB[7,0]
reference_listCamB=reference_CamB.split()
reference_valuesCamB=[eval(reference_listCamB[0]),eval(reference_listCamB[1]),eval(reference_listCamB[2]),eval(reference_listCamB[3])]

chi2_newCamB=chi2_stringCamB.replace('Sqrt','np.sqrt').replace('Abs','np.abs').replace('Sign','np.sign').replace('Log','np.log').replace('Csc','1/math.sin').replace('Pi','math.pi').replace('Sin','math.sin').replace('E**','np.exp')

def chi2_CamB(a0,YD,df2,Dinc):
   chi2_formula_CamB=eval(chi2_newCamB)
   return chi2_formula_CamB

def log_prior_CamB(theta):

 log10_a0,log10_YD,df2,Dinc=theta
   
 if ((-20<log10_a0<0) and (Dmin/reference_valuesCamB[0]<df2) and (imin-reference_valuesCamB[2]<Dinc<imax-reference_valuesCamB[2])):       
    return 0.0
 return -np.inf

def log_likelihood_CamB(theta):
 log10_a0,log10_YD,df2,Dinc=theta
 
 loglikelihood_CamB = -0.5*chi2_CamB(10**log10_a0,10**log10_YD,df2,Dinc)     
 
 return loglikelihood_CamB

def log_probability_CamB(theta):
    lp = log_prior_CamB(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood_CamB(theta)

my_guess_CamB=[np.log10(best_fit_valuesCamB[0]),np.log10(best_fit_valuesCamB[1]),best_fit_valuesCamB[2],best_fit_valuesCamB[3]]
pos_CamB = my_guess_CamB + 1e-14*np.random.randn(100, len(my_guess_CamB))
nwalkers_CamB, ndim_CamB = pos_CamB.shape

####################################################################################  MPI #####################################################################################################################################################
import sys
import time
from schwimmbad import MPIPool

with MPIPool() as pool:
  if not pool.is_master():
    pool.wait()
    sys.exit(0)
  
  sampler_CamB = emcee.EnsembleSampler(nwalkers_CamB, ndim_CamB, log_probability_CamB, moves= my_moves,pool=pool)
  sampler_CamB.run_mcmc(pos_CamB, nsteps,progress=True)

##################################################### EXPORTING #################################################################################################################################################################
chain_CamB=sampler_CamB.get_chain()[:, :, 0].T
flat_chain_CamB = sampler_CamB.get_chain(discard=burn_in,thin=thin,flat=True)
import os
flat_chain_export_CamB = pd.DataFrame(flat_chain_CamB)
flat_chain_export_CamB.to_csv('/home/alefe/Cluster/MOND_outputs/chain_CamB.csv')