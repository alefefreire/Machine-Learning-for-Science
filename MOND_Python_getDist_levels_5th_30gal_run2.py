from __future__ import print_function
import numpy 
import pandas as pd
import csv
import matplotlib.pyplot as plot
from matplotlib.ticker import MaxNLocator
from scipy.interpolate import interp1d
import scipy.optimize as op
import scipy
from scipy import *
from scipy.special import expi
#import lmfit
#from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit
#import emcee
#from emcee import PTSampler
from matplotlib import rcParams
import time
from scipy.interpolate import InterpolatedUnivariateSpline


import sys
sys.path.insert(0,r'c:\work\dist\git\getdist')
import getdist
from getdist import plots, MCSamples
#import getdist, IPython
import pylab as plt
print('GetDist Version: %s, Matplotlib version: %s'%(getdist.__version__, plt.matplotlib.__version__))
#matplotlib 2 doesn't seem to work well without usetex on
plt.rcParams['text.usetex']=True

labels=[r"$\log_{10} a_0$",r"$\log_{10}\Upsilon_{*d}$",r"$\delta$", r"$\Delta i$"]
labels_with_bulge=[r"$\log_{10} a_0$",r"$\log_{10}\Upsilon_{*d}$",r"$\log_{10}\Upsilon_{*b}$",r"$\delta$", r"$\Delta i$"]
names=["log10_a0","log10_YD","df2","Dinc"]
names_with_bulge=["log10_a0","log10_YD","log10_YB","df2","Dinc"]
ndim=len(labels)
ndim_with_bulge=len(labels_with_bulge)



full_chain_UGC06930 = pd.read_csv('/home/alefe/Cluster/MOND_outputs/chain_UGC06930.csv',usecols=[1,2,3,4])
data_array_chain_UGC06930 = numpy.array(full_chain_UGC06930)
mle_soln_UGC06930=[]

for i in range(ndim):
    mcmc_UGC06930 = numpy.percentile(data_array_chain_UGC06930[:, i], [16, 50, 84])
    q_UGC06930 = numpy.diff(mcmc_UGC06930)
    mle_soln_UGC06930.append(mcmc_UGC06930[1])

log10_a0_sol_UGC06930=mle_soln_UGC06930[0]
log10_YD_sol_UGC06930=mle_soln_UGC06930[1]
df2_sol_UGC06930=mle_soln_UGC06930[2]
Dinc_sol_UGC06930=mle_soln_UGC06930[3]

samp_UGC06930 = getdist.mcsamples.MCSamples(samples=data_array_chain_UGC06930,names = names, labels = labels)
samp_UGC06930.updateSettings({'contours': [0.682689492137086, 0.954499736103642,0.997300203936740,0.999936657516334,0.999999426696856]})
stats_UGC06930 = samp_UGC06930.getMargeStats()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC06930-sigmatab-log10_a0.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_a0_sol_str='%s'%log10_a0_sol_UGC06930
  low=log10_a0_sol_UGC06930-stats_UGC06930.parWithName('log10_a0').limits[i].lower
  low_str='%s'%low
  up=stats_UGC06930.parWithName('log10_a0').limits[i].upper- log10_a0_sol_UGC06930
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC06930-sigmatab-log10_a0.txt', 'a')
  text_file.write(log10_a0_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC06930-sigmatab-log10_YD.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_YD_sol_str='%s'%log10_YD_sol_UGC06930
  low=log10_YD_sol_UGC06930-stats_UGC06930.parWithName('log10_YD').limits[i].lower
  low_str='%s'%low
  up=stats_UGC06930.parWithName('log10_YD').limits[i].upper- log10_YD_sol_UGC06930
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC06930-sigmatab-log10_YD.txt', 'a')
  text_file.write(log10_YD_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC06930-sigmatab-df2.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  df2_sol_str='%s'%df2_sol_UGC06930
  low=df2_sol_UGC06930-stats_UGC06930.parWithName('df2').limits[i].lower
  low_str='%s'%low
  up=stats_UGC06930.parWithName('df2').limits[i].upper- df2_sol_UGC06930
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC06930-sigmatab-df2.txt', 'a')
  text_file.write(df2_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC06930-sigmatab-Dinc.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  Dinc_sol_str='%s'%Dinc_sol_UGC06930
  low=Dinc_sol_UGC06930-stats_UGC06930.parWithName('Dinc').limits[i].lower
  low_str='%s'%low
  up=stats_UGC06930.parWithName('Dinc').limits[i].upper- Dinc_sol_UGC06930
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC06930-sigmatab-Dinc.txt', 'a')
  text_file.write(Dinc_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()


full_chain_UGC06983 = pd.read_csv('/home/alefe/Cluster/MOND_outputs/chain_UGC06983.csv',usecols=[1,2,3,4])
data_array_chain_UGC06983 = numpy.array(full_chain_UGC06983)
mle_soln_UGC06983=[]

for i in range(ndim):
    mcmc_UGC06983 = numpy.percentile(data_array_chain_UGC06983[:, i], [16, 50, 84])
    q_UGC06983 = numpy.diff(mcmc_UGC06983)
    mle_soln_UGC06983.append(mcmc_UGC06983[1])

log10_a0_sol_UGC06983=mle_soln_UGC06983[0]
log10_YD_sol_UGC06983=mle_soln_UGC06983[1]
df2_sol_UGC06983=mle_soln_UGC06983[2]
Dinc_sol_UGC06983=mle_soln_UGC06983[3]

samp_UGC06983 = getdist.mcsamples.MCSamples(samples=data_array_chain_UGC06983,names = names, labels = labels)
samp_UGC06983.updateSettings({'contours': [0.682689492137086, 0.954499736103642,0.997300203936740,0.999936657516334,0.999999426696856]})
stats_UGC06983 = samp_UGC06983.getMargeStats()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC06983-sigmatab-log10_a0.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_a0_sol_str='%s'%log10_a0_sol_UGC06983
  low=log10_a0_sol_UGC06983-stats_UGC06983.parWithName('log10_a0').limits[i].lower
  low_str='%s'%low
  up=stats_UGC06983.parWithName('log10_a0').limits[i].upper- log10_a0_sol_UGC06983
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC06983-sigmatab-log10_a0.txt', 'a')
  text_file.write(log10_a0_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC06983-sigmatab-log10_YD.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_YD_sol_str='%s'%log10_YD_sol_UGC06983
  low=log10_YD_sol_UGC06983-stats_UGC06983.parWithName('log10_YD').limits[i].lower
  low_str='%s'%low
  up=stats_UGC06983.parWithName('log10_YD').limits[i].upper- log10_YD_sol_UGC06983
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC06983-sigmatab-log10_YD.txt', 'a')
  text_file.write(log10_YD_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC06983-sigmatab-df2.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  df2_sol_str='%s'%df2_sol_UGC06983
  low=df2_sol_UGC06983-stats_UGC06983.parWithName('df2').limits[i].lower
  low_str='%s'%low
  up=stats_UGC06983.parWithName('df2').limits[i].upper- df2_sol_UGC06983
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC06983-sigmatab-df2.txt', 'a')
  text_file.write(df2_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC06983-sigmatab-Dinc.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  Dinc_sol_str='%s'%Dinc_sol_UGC06983
  low=Dinc_sol_UGC06983-stats_UGC06983.parWithName('Dinc').limits[i].lower
  low_str='%s'%low
  up=stats_UGC06983.parWithName('Dinc').limits[i].upper- Dinc_sol_UGC06983
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC06983-sigmatab-Dinc.txt', 'a')
  text_file.write(Dinc_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()


full_chain_UGC07089 = pd.read_csv('/home/alefe/Cluster/MOND_outputs/chain_UGC07089.csv',usecols=[1,2,3,4])
data_array_chain_UGC07089 = numpy.array(full_chain_UGC07089)
mle_soln_UGC07089=[]

for i in range(ndim):
    mcmc_UGC07089 = numpy.percentile(data_array_chain_UGC07089[:, i], [16, 50, 84])
    q_UGC07089 = numpy.diff(mcmc_UGC07089)
    mle_soln_UGC07089.append(mcmc_UGC07089[1])

log10_a0_sol_UGC07089=mle_soln_UGC07089[0]
log10_YD_sol_UGC07089=mle_soln_UGC07089[1]
df2_sol_UGC07089=mle_soln_UGC07089[2]
Dinc_sol_UGC07089=mle_soln_UGC07089[3]

samp_UGC07089 = getdist.mcsamples.MCSamples(samples=data_array_chain_UGC07089,names = names, labels = labels)
samp_UGC07089.updateSettings({'contours': [0.682689492137086, 0.954499736103642,0.997300203936740,0.999936657516334,0.999999426696856]})
stats_UGC07089 = samp_UGC07089.getMargeStats()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07089-sigmatab-log10_a0.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_a0_sol_str='%s'%log10_a0_sol_UGC07089
  low=log10_a0_sol_UGC07089-stats_UGC07089.parWithName('log10_a0').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07089.parWithName('log10_a0').limits[i].upper- log10_a0_sol_UGC07089
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07089-sigmatab-log10_a0.txt', 'a')
  text_file.write(log10_a0_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07089-sigmatab-log10_YD.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_YD_sol_str='%s'%log10_YD_sol_UGC07089
  low=log10_YD_sol_UGC07089-stats_UGC07089.parWithName('log10_YD').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07089.parWithName('log10_YD').limits[i].upper- log10_YD_sol_UGC07089
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07089-sigmatab-log10_YD.txt', 'a')
  text_file.write(log10_YD_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07089-sigmatab-df2.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  df2_sol_str='%s'%df2_sol_UGC07089
  low=df2_sol_UGC07089-stats_UGC07089.parWithName('df2').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07089.parWithName('df2').limits[i].upper- df2_sol_UGC07089
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07089-sigmatab-df2.txt', 'a')
  text_file.write(df2_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07089-sigmatab-Dinc.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  Dinc_sol_str='%s'%Dinc_sol_UGC07089
  low=Dinc_sol_UGC07089-stats_UGC07089.parWithName('Dinc').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07089.parWithName('Dinc').limits[i].upper- Dinc_sol_UGC07089
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07089-sigmatab-Dinc.txt', 'a')
  text_file.write(Dinc_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()


full_chain_UGC07125 = pd.read_csv('/home/alefe/Cluster/MOND_outputs/chain_UGC07125.csv',usecols=[1,2,3,4])
data_array_chain_UGC07125 = numpy.array(full_chain_UGC07125)
mle_soln_UGC07125=[]

for i in range(ndim):
    mcmc_UGC07125 = numpy.percentile(data_array_chain_UGC07125[:, i], [16, 50, 84])
    q_UGC07125 = numpy.diff(mcmc_UGC07125)
    mle_soln_UGC07125.append(mcmc_UGC07125[1])

log10_a0_sol_UGC07125=mle_soln_UGC07125[0]
log10_YD_sol_UGC07125=mle_soln_UGC07125[1]
df2_sol_UGC07125=mle_soln_UGC07125[2]
Dinc_sol_UGC07125=mle_soln_UGC07125[3]

samp_UGC07125 = getdist.mcsamples.MCSamples(samples=data_array_chain_UGC07125,names = names, labels = labels)
samp_UGC07125.updateSettings({'contours': [0.682689492137086, 0.954499736103642,0.997300203936740,0.999936657516334,0.999999426696856]})
stats_UGC07125 = samp_UGC07125.getMargeStats()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07125-sigmatab-log10_a0.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_a0_sol_str='%s'%log10_a0_sol_UGC07125
  low=log10_a0_sol_UGC07125-stats_UGC07125.parWithName('log10_a0').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07125.parWithName('log10_a0').limits[i].upper- log10_a0_sol_UGC07125
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07125-sigmatab-log10_a0.txt', 'a')
  text_file.write(log10_a0_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07125-sigmatab-log10_YD.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_YD_sol_str='%s'%log10_YD_sol_UGC07125
  low=log10_YD_sol_UGC07125-stats_UGC07125.parWithName('log10_YD').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07125.parWithName('log10_YD').limits[i].upper- log10_YD_sol_UGC07125
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07125-sigmatab-log10_YD.txt', 'a')
  text_file.write(log10_YD_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07125-sigmatab-df2.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  df2_sol_str='%s'%df2_sol_UGC07125
  low=df2_sol_UGC07125-stats_UGC07125.parWithName('df2').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07125.parWithName('df2').limits[i].upper- df2_sol_UGC07125
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07125-sigmatab-df2.txt', 'a')
  text_file.write(df2_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07125-sigmatab-Dinc.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  Dinc_sol_str='%s'%Dinc_sol_UGC07125
  low=Dinc_sol_UGC07125-stats_UGC07125.parWithName('Dinc').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07125.parWithName('Dinc').limits[i].upper- Dinc_sol_UGC07125
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07125-sigmatab-Dinc.txt', 'a')
  text_file.write(Dinc_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()


full_chain_UGC07151 = pd.read_csv('/home/alefe/Cluster/MOND_outputs/chain_UGC07151.csv',usecols=[1,2,3,4])
data_array_chain_UGC07151 = numpy.array(full_chain_UGC07151)
mle_soln_UGC07151=[]

for i in range(ndim):
    mcmc_UGC07151 = numpy.percentile(data_array_chain_UGC07151[:, i], [16, 50, 84])
    q_UGC07151 = numpy.diff(mcmc_UGC07151)
    mle_soln_UGC07151.append(mcmc_UGC07151[1])

log10_a0_sol_UGC07151=mle_soln_UGC07151[0]
log10_YD_sol_UGC07151=mle_soln_UGC07151[1]
df2_sol_UGC07151=mle_soln_UGC07151[2]
Dinc_sol_UGC07151=mle_soln_UGC07151[3]

samp_UGC07151 = getdist.mcsamples.MCSamples(samples=data_array_chain_UGC07151,names = names, labels = labels)
samp_UGC07151.updateSettings({'contours': [0.682689492137086, 0.954499736103642,0.997300203936740,0.999936657516334,0.999999426696856]})
stats_UGC07151 = samp_UGC07151.getMargeStats()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07151-sigmatab-log10_a0.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_a0_sol_str='%s'%log10_a0_sol_UGC07151
  low=log10_a0_sol_UGC07151-stats_UGC07151.parWithName('log10_a0').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07151.parWithName('log10_a0').limits[i].upper- log10_a0_sol_UGC07151
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07151-sigmatab-log10_a0.txt', 'a')
  text_file.write(log10_a0_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07151-sigmatab-log10_YD.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_YD_sol_str='%s'%log10_YD_sol_UGC07151
  low=log10_YD_sol_UGC07151-stats_UGC07151.parWithName('log10_YD').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07151.parWithName('log10_YD').limits[i].upper- log10_YD_sol_UGC07151
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07151-sigmatab-log10_YD.txt', 'a')
  text_file.write(log10_YD_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07151-sigmatab-df2.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  df2_sol_str='%s'%df2_sol_UGC07151
  low=df2_sol_UGC07151-stats_UGC07151.parWithName('df2').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07151.parWithName('df2').limits[i].upper- df2_sol_UGC07151
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07151-sigmatab-df2.txt', 'a')
  text_file.write(df2_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07151-sigmatab-Dinc.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  Dinc_sol_str='%s'%Dinc_sol_UGC07151
  low=Dinc_sol_UGC07151-stats_UGC07151.parWithName('Dinc').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07151.parWithName('Dinc').limits[i].upper- Dinc_sol_UGC07151
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07151-sigmatab-Dinc.txt', 'a')
  text_file.write(Dinc_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()


full_chain_UGC07232 = pd.read_csv('/home/alefe/Cluster/MOND_outputs/chain_UGC07232.csv',usecols=[1,2,3,4])
data_array_chain_UGC07232 = numpy.array(full_chain_UGC07232)
mle_soln_UGC07232=[]

for i in range(ndim):
    mcmc_UGC07232 = numpy.percentile(data_array_chain_UGC07232[:, i], [16, 50, 84])
    q_UGC07232 = numpy.diff(mcmc_UGC07232)
    mle_soln_UGC07232.append(mcmc_UGC07232[1])

log10_a0_sol_UGC07232=mle_soln_UGC07232[0]
log10_YD_sol_UGC07232=mle_soln_UGC07232[1]
df2_sol_UGC07232=mle_soln_UGC07232[2]
Dinc_sol_UGC07232=mle_soln_UGC07232[3]

samp_UGC07232 = getdist.mcsamples.MCSamples(samples=data_array_chain_UGC07232,names = names, labels = labels)
samp_UGC07232.updateSettings({'contours': [0.682689492137086, 0.954499736103642,0.997300203936740,0.999936657516334,0.999999426696856]})
stats_UGC07232 = samp_UGC07232.getMargeStats()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07232-sigmatab-log10_a0.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_a0_sol_str='%s'%log10_a0_sol_UGC07232
  low=log10_a0_sol_UGC07232-stats_UGC07232.parWithName('log10_a0').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07232.parWithName('log10_a0').limits[i].upper- log10_a0_sol_UGC07232
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07232-sigmatab-log10_a0.txt', 'a')
  text_file.write(log10_a0_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07232-sigmatab-log10_YD.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_YD_sol_str='%s'%log10_YD_sol_UGC07232
  low=log10_YD_sol_UGC07232-stats_UGC07232.parWithName('log10_YD').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07232.parWithName('log10_YD').limits[i].upper- log10_YD_sol_UGC07232
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07232-sigmatab-log10_YD.txt', 'a')
  text_file.write(log10_YD_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07232-sigmatab-df2.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  df2_sol_str='%s'%df2_sol_UGC07232
  low=df2_sol_UGC07232-stats_UGC07232.parWithName('df2').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07232.parWithName('df2').limits[i].upper- df2_sol_UGC07232
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07232-sigmatab-df2.txt', 'a')
  text_file.write(df2_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07232-sigmatab-Dinc.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  Dinc_sol_str='%s'%Dinc_sol_UGC07232
  low=Dinc_sol_UGC07232-stats_UGC07232.parWithName('Dinc').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07232.parWithName('Dinc').limits[i].upper- Dinc_sol_UGC07232
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07232-sigmatab-Dinc.txt', 'a')
  text_file.write(Dinc_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()


full_chain_UGC07261 = pd.read_csv('/home/alefe/Cluster/MOND_outputs/chain_UGC07261.csv',usecols=[1,2,3,4])
data_array_chain_UGC07261 = numpy.array(full_chain_UGC07261)
mle_soln_UGC07261=[]

for i in range(ndim):
    mcmc_UGC07261 = numpy.percentile(data_array_chain_UGC07261[:, i], [16, 50, 84])
    q_UGC07261 = numpy.diff(mcmc_UGC07261)
    mle_soln_UGC07261.append(mcmc_UGC07261[1])

log10_a0_sol_UGC07261=mle_soln_UGC07261[0]
log10_YD_sol_UGC07261=mle_soln_UGC07261[1]
df2_sol_UGC07261=mle_soln_UGC07261[2]
Dinc_sol_UGC07261=mle_soln_UGC07261[3]

samp_UGC07261 = getdist.mcsamples.MCSamples(samples=data_array_chain_UGC07261,names = names, labels = labels)
samp_UGC07261.updateSettings({'contours': [0.682689492137086, 0.954499736103642,0.997300203936740,0.999936657516334,0.999999426696856]})
stats_UGC07261 = samp_UGC07261.getMargeStats()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07261-sigmatab-log10_a0.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_a0_sol_str='%s'%log10_a0_sol_UGC07261
  low=log10_a0_sol_UGC07261-stats_UGC07261.parWithName('log10_a0').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07261.parWithName('log10_a0').limits[i].upper- log10_a0_sol_UGC07261
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07261-sigmatab-log10_a0.txt', 'a')
  text_file.write(log10_a0_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07261-sigmatab-log10_YD.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_YD_sol_str='%s'%log10_YD_sol_UGC07261
  low=log10_YD_sol_UGC07261-stats_UGC07261.parWithName('log10_YD').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07261.parWithName('log10_YD').limits[i].upper- log10_YD_sol_UGC07261
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07261-sigmatab-log10_YD.txt', 'a')
  text_file.write(log10_YD_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07261-sigmatab-df2.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  df2_sol_str='%s'%df2_sol_UGC07261
  low=df2_sol_UGC07261-stats_UGC07261.parWithName('df2').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07261.parWithName('df2').limits[i].upper- df2_sol_UGC07261
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07261-sigmatab-df2.txt', 'a')
  text_file.write(df2_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07261-sigmatab-Dinc.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  Dinc_sol_str='%s'%Dinc_sol_UGC07261
  low=Dinc_sol_UGC07261-stats_UGC07261.parWithName('Dinc').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07261.parWithName('Dinc').limits[i].upper- Dinc_sol_UGC07261
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07261-sigmatab-Dinc.txt', 'a')
  text_file.write(Dinc_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()


full_chain_UGC07323 = pd.read_csv('/home/alefe/Cluster/MOND_outputs/chain_UGC07323.csv',usecols=[1,2,3,4])
data_array_chain_UGC07323 = numpy.array(full_chain_UGC07323)
mle_soln_UGC07323=[]

for i in range(ndim):
    mcmc_UGC07323 = numpy.percentile(data_array_chain_UGC07323[:, i], [16, 50, 84])
    q_UGC07323 = numpy.diff(mcmc_UGC07323)
    mle_soln_UGC07323.append(mcmc_UGC07323[1])

log10_a0_sol_UGC07323=mle_soln_UGC07323[0]
log10_YD_sol_UGC07323=mle_soln_UGC07323[1]
df2_sol_UGC07323=mle_soln_UGC07323[2]
Dinc_sol_UGC07323=mle_soln_UGC07323[3]

samp_UGC07323 = getdist.mcsamples.MCSamples(samples=data_array_chain_UGC07323,names = names, labels = labels)
samp_UGC07323.updateSettings({'contours': [0.682689492137086, 0.954499736103642,0.997300203936740,0.999936657516334,0.999999426696856]})
stats_UGC07323 = samp_UGC07323.getMargeStats()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07323-sigmatab-log10_a0.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_a0_sol_str='%s'%log10_a0_sol_UGC07323
  low=log10_a0_sol_UGC07323-stats_UGC07323.parWithName('log10_a0').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07323.parWithName('log10_a0').limits[i].upper- log10_a0_sol_UGC07323
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07323-sigmatab-log10_a0.txt', 'a')
  text_file.write(log10_a0_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07323-sigmatab-log10_YD.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_YD_sol_str='%s'%log10_YD_sol_UGC07323
  low=log10_YD_sol_UGC07323-stats_UGC07323.parWithName('log10_YD').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07323.parWithName('log10_YD').limits[i].upper- log10_YD_sol_UGC07323
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07323-sigmatab-log10_YD.txt', 'a')
  text_file.write(log10_YD_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07323-sigmatab-df2.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  df2_sol_str='%s'%df2_sol_UGC07323
  low=df2_sol_UGC07323-stats_UGC07323.parWithName('df2').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07323.parWithName('df2').limits[i].upper- df2_sol_UGC07323
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07323-sigmatab-df2.txt', 'a')
  text_file.write(df2_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07323-sigmatab-Dinc.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  Dinc_sol_str='%s'%Dinc_sol_UGC07323
  low=Dinc_sol_UGC07323-stats_UGC07323.parWithName('Dinc').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07323.parWithName('Dinc').limits[i].upper- Dinc_sol_UGC07323
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07323-sigmatab-Dinc.txt', 'a')
  text_file.write(Dinc_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()


full_chain_UGC07399 = pd.read_csv('/home/alefe/Cluster/MOND_outputs/chain_UGC07399.csv',usecols=[1,2,3,4])
data_array_chain_UGC07399 = numpy.array(full_chain_UGC07399)
mle_soln_UGC07399=[]

for i in range(ndim):
    mcmc_UGC07399 = numpy.percentile(data_array_chain_UGC07399[:, i], [16, 50, 84])
    q_UGC07399 = numpy.diff(mcmc_UGC07399)
    mle_soln_UGC07399.append(mcmc_UGC07399[1])

log10_a0_sol_UGC07399=mle_soln_UGC07399[0]
log10_YD_sol_UGC07399=mle_soln_UGC07399[1]
df2_sol_UGC07399=mle_soln_UGC07399[2]
Dinc_sol_UGC07399=mle_soln_UGC07399[3]

samp_UGC07399 = getdist.mcsamples.MCSamples(samples=data_array_chain_UGC07399,names = names, labels = labels)
samp_UGC07399.updateSettings({'contours': [0.682689492137086, 0.954499736103642,0.997300203936740,0.999936657516334,0.999999426696856]})
stats_UGC07399 = samp_UGC07399.getMargeStats()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07399-sigmatab-log10_a0.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_a0_sol_str='%s'%log10_a0_sol_UGC07399
  low=log10_a0_sol_UGC07399-stats_UGC07399.parWithName('log10_a0').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07399.parWithName('log10_a0').limits[i].upper- log10_a0_sol_UGC07399
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07399-sigmatab-log10_a0.txt', 'a')
  text_file.write(log10_a0_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07399-sigmatab-log10_YD.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_YD_sol_str='%s'%log10_YD_sol_UGC07399
  low=log10_YD_sol_UGC07399-stats_UGC07399.parWithName('log10_YD').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07399.parWithName('log10_YD').limits[i].upper- log10_YD_sol_UGC07399
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07399-sigmatab-log10_YD.txt', 'a')
  text_file.write(log10_YD_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07399-sigmatab-df2.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  df2_sol_str='%s'%df2_sol_UGC07399
  low=df2_sol_UGC07399-stats_UGC07399.parWithName('df2').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07399.parWithName('df2').limits[i].upper- df2_sol_UGC07399
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07399-sigmatab-df2.txt', 'a')
  text_file.write(df2_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07399-sigmatab-Dinc.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  Dinc_sol_str='%s'%Dinc_sol_UGC07399
  low=Dinc_sol_UGC07399-stats_UGC07399.parWithName('Dinc').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07399.parWithName('Dinc').limits[i].upper- Dinc_sol_UGC07399
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07399-sigmatab-Dinc.txt', 'a')
  text_file.write(Dinc_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()


full_chain_UGC07524 = pd.read_csv('/home/alefe/Cluster/MOND_outputs/chain_UGC07524.csv',usecols=[1,2,3,4])
data_array_chain_UGC07524 = numpy.array(full_chain_UGC07524)
mle_soln_UGC07524=[]

for i in range(ndim):
    mcmc_UGC07524 = numpy.percentile(data_array_chain_UGC07524[:, i], [16, 50, 84])
    q_UGC07524 = numpy.diff(mcmc_UGC07524)
    mle_soln_UGC07524.append(mcmc_UGC07524[1])

log10_a0_sol_UGC07524=mle_soln_UGC07524[0]
log10_YD_sol_UGC07524=mle_soln_UGC07524[1]
df2_sol_UGC07524=mle_soln_UGC07524[2]
Dinc_sol_UGC07524=mle_soln_UGC07524[3]

samp_UGC07524 = getdist.mcsamples.MCSamples(samples=data_array_chain_UGC07524,names = names, labels = labels)
samp_UGC07524.updateSettings({'contours': [0.682689492137086, 0.954499736103642,0.997300203936740,0.999936657516334,0.999999426696856]})
stats_UGC07524 = samp_UGC07524.getMargeStats()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07524-sigmatab-log10_a0.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_a0_sol_str='%s'%log10_a0_sol_UGC07524
  low=log10_a0_sol_UGC07524-stats_UGC07524.parWithName('log10_a0').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07524.parWithName('log10_a0').limits[i].upper- log10_a0_sol_UGC07524
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07524-sigmatab-log10_a0.txt', 'a')
  text_file.write(log10_a0_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07524-sigmatab-log10_YD.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_YD_sol_str='%s'%log10_YD_sol_UGC07524
  low=log10_YD_sol_UGC07524-stats_UGC07524.parWithName('log10_YD').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07524.parWithName('log10_YD').limits[i].upper- log10_YD_sol_UGC07524
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07524-sigmatab-log10_YD.txt', 'a')
  text_file.write(log10_YD_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07524-sigmatab-df2.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  df2_sol_str='%s'%df2_sol_UGC07524
  low=df2_sol_UGC07524-stats_UGC07524.parWithName('df2').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07524.parWithName('df2').limits[i].upper- df2_sol_UGC07524
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07524-sigmatab-df2.txt', 'a')
  text_file.write(df2_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07524-sigmatab-Dinc.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  Dinc_sol_str='%s'%Dinc_sol_UGC07524
  low=Dinc_sol_UGC07524-stats_UGC07524.parWithName('Dinc').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07524.parWithName('Dinc').limits[i].upper- Dinc_sol_UGC07524
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07524-sigmatab-Dinc.txt', 'a')
  text_file.write(Dinc_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()


full_chain_UGC07559 = pd.read_csv('/home/alefe/Cluster/MOND_outputs/chain_UGC07559.csv',usecols=[1,2,3,4])
data_array_chain_UGC07559 = numpy.array(full_chain_UGC07559)
mle_soln_UGC07559=[]

for i in range(ndim):
    mcmc_UGC07559 = numpy.percentile(data_array_chain_UGC07559[:, i], [16, 50, 84])
    q_UGC07559 = numpy.diff(mcmc_UGC07559)
    mle_soln_UGC07559.append(mcmc_UGC07559[1])

log10_a0_sol_UGC07559=mle_soln_UGC07559[0]
log10_YD_sol_UGC07559=mle_soln_UGC07559[1]
df2_sol_UGC07559=mle_soln_UGC07559[2]
Dinc_sol_UGC07559=mle_soln_UGC07559[3]

samp_UGC07559 = getdist.mcsamples.MCSamples(samples=data_array_chain_UGC07559,names = names, labels = labels)
samp_UGC07559.updateSettings({'contours': [0.682689492137086, 0.954499736103642,0.997300203936740,0.999936657516334,0.999999426696856]})
stats_UGC07559 = samp_UGC07559.getMargeStats()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07559-sigmatab-log10_a0.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_a0_sol_str='%s'%log10_a0_sol_UGC07559
  low=log10_a0_sol_UGC07559-stats_UGC07559.parWithName('log10_a0').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07559.parWithName('log10_a0').limits[i].upper- log10_a0_sol_UGC07559
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07559-sigmatab-log10_a0.txt', 'a')
  text_file.write(log10_a0_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07559-sigmatab-log10_YD.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_YD_sol_str='%s'%log10_YD_sol_UGC07559
  low=log10_YD_sol_UGC07559-stats_UGC07559.parWithName('log10_YD').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07559.parWithName('log10_YD').limits[i].upper- log10_YD_sol_UGC07559
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07559-sigmatab-log10_YD.txt', 'a')
  text_file.write(log10_YD_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07559-sigmatab-df2.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  df2_sol_str='%s'%df2_sol_UGC07559
  low=df2_sol_UGC07559-stats_UGC07559.parWithName('df2').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07559.parWithName('df2').limits[i].upper- df2_sol_UGC07559
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07559-sigmatab-df2.txt', 'a')
  text_file.write(df2_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07559-sigmatab-Dinc.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  Dinc_sol_str='%s'%Dinc_sol_UGC07559
  low=Dinc_sol_UGC07559-stats_UGC07559.parWithName('Dinc').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07559.parWithName('Dinc').limits[i].upper- Dinc_sol_UGC07559
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07559-sigmatab-Dinc.txt', 'a')
  text_file.write(Dinc_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()


full_chain_UGC07577 = pd.read_csv('/home/alefe/Cluster/MOND_outputs/chain_UGC07577.csv',usecols=[1,2,3,4])
data_array_chain_UGC07577 = numpy.array(full_chain_UGC07577)
mle_soln_UGC07577=[]

for i in range(ndim):
    mcmc_UGC07577 = numpy.percentile(data_array_chain_UGC07577[:, i], [16, 50, 84])
    q_UGC07577 = numpy.diff(mcmc_UGC07577)
    mle_soln_UGC07577.append(mcmc_UGC07577[1])

log10_a0_sol_UGC07577=mle_soln_UGC07577[0]
log10_YD_sol_UGC07577=mle_soln_UGC07577[1]
df2_sol_UGC07577=mle_soln_UGC07577[2]
Dinc_sol_UGC07577=mle_soln_UGC07577[3]

samp_UGC07577 = getdist.mcsamples.MCSamples(samples=data_array_chain_UGC07577,names = names, labels = labels)
samp_UGC07577.updateSettings({'contours': [0.682689492137086, 0.954499736103642,0.997300203936740,0.999936657516334,0.999999426696856]})
stats_UGC07577 = samp_UGC07577.getMargeStats()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07577-sigmatab-log10_a0.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_a0_sol_str='%s'%log10_a0_sol_UGC07577
  low=log10_a0_sol_UGC07577-stats_UGC07577.parWithName('log10_a0').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07577.parWithName('log10_a0').limits[i].upper- log10_a0_sol_UGC07577
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07577-sigmatab-log10_a0.txt', 'a')
  text_file.write(log10_a0_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07577-sigmatab-log10_YD.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_YD_sol_str='%s'%log10_YD_sol_UGC07577
  low=log10_YD_sol_UGC07577-stats_UGC07577.parWithName('log10_YD').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07577.parWithName('log10_YD').limits[i].upper- log10_YD_sol_UGC07577
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07577-sigmatab-log10_YD.txt', 'a')
  text_file.write(log10_YD_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07577-sigmatab-df2.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  df2_sol_str='%s'%df2_sol_UGC07577
  low=df2_sol_UGC07577-stats_UGC07577.parWithName('df2').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07577.parWithName('df2').limits[i].upper- df2_sol_UGC07577
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07577-sigmatab-df2.txt', 'a')
  text_file.write(df2_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07577-sigmatab-Dinc.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  Dinc_sol_str='%s'%Dinc_sol_UGC07577
  low=Dinc_sol_UGC07577-stats_UGC07577.parWithName('Dinc').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07577.parWithName('Dinc').limits[i].upper- Dinc_sol_UGC07577
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07577-sigmatab-Dinc.txt', 'a')
  text_file.write(Dinc_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()


full_chain_UGC07603 = pd.read_csv('/home/alefe/Cluster/MOND_outputs/chain_UGC07603.csv',usecols=[1,2,3,4])
data_array_chain_UGC07603 = numpy.array(full_chain_UGC07603)
mle_soln_UGC07603=[]

for i in range(ndim):
    mcmc_UGC07603 = numpy.percentile(data_array_chain_UGC07603[:, i], [16, 50, 84])
    q_UGC07603 = numpy.diff(mcmc_UGC07603)
    mle_soln_UGC07603.append(mcmc_UGC07603[1])

log10_a0_sol_UGC07603=mle_soln_UGC07603[0]
log10_YD_sol_UGC07603=mle_soln_UGC07603[1]
df2_sol_UGC07603=mle_soln_UGC07603[2]
Dinc_sol_UGC07603=mle_soln_UGC07603[3]

samp_UGC07603 = getdist.mcsamples.MCSamples(samples=data_array_chain_UGC07603,names = names, labels = labels)
samp_UGC07603.updateSettings({'contours': [0.682689492137086, 0.954499736103642,0.997300203936740,0.999936657516334,0.999999426696856]})
stats_UGC07603 = samp_UGC07603.getMargeStats()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07603-sigmatab-log10_a0.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_a0_sol_str='%s'%log10_a0_sol_UGC07603
  low=log10_a0_sol_UGC07603-stats_UGC07603.parWithName('log10_a0').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07603.parWithName('log10_a0').limits[i].upper- log10_a0_sol_UGC07603
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07603-sigmatab-log10_a0.txt', 'a')
  text_file.write(log10_a0_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07603-sigmatab-log10_YD.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_YD_sol_str='%s'%log10_YD_sol_UGC07603
  low=log10_YD_sol_UGC07603-stats_UGC07603.parWithName('log10_YD').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07603.parWithName('log10_YD').limits[i].upper- log10_YD_sol_UGC07603
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07603-sigmatab-log10_YD.txt', 'a')
  text_file.write(log10_YD_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07603-sigmatab-df2.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  df2_sol_str='%s'%df2_sol_UGC07603
  low=df2_sol_UGC07603-stats_UGC07603.parWithName('df2').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07603.parWithName('df2').limits[i].upper- df2_sol_UGC07603
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07603-sigmatab-df2.txt', 'a')
  text_file.write(df2_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07603-sigmatab-Dinc.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  Dinc_sol_str='%s'%Dinc_sol_UGC07603
  low=Dinc_sol_UGC07603-stats_UGC07603.parWithName('Dinc').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07603.parWithName('Dinc').limits[i].upper- Dinc_sol_UGC07603
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07603-sigmatab-Dinc.txt', 'a')
  text_file.write(Dinc_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()


full_chain_UGC07690 = pd.read_csv('/home/alefe/Cluster/MOND_outputs/chain_UGC07690.csv',usecols=[1,2,3,4])
data_array_chain_UGC07690 = numpy.array(full_chain_UGC07690)
mle_soln_UGC07690=[]

for i in range(ndim):
    mcmc_UGC07690 = numpy.percentile(data_array_chain_UGC07690[:, i], [16, 50, 84])
    q_UGC07690 = numpy.diff(mcmc_UGC07690)
    mle_soln_UGC07690.append(mcmc_UGC07690[1])

log10_a0_sol_UGC07690=mle_soln_UGC07690[0]
log10_YD_sol_UGC07690=mle_soln_UGC07690[1]
df2_sol_UGC07690=mle_soln_UGC07690[2]
Dinc_sol_UGC07690=mle_soln_UGC07690[3]

samp_UGC07690 = getdist.mcsamples.MCSamples(samples=data_array_chain_UGC07690,names = names, labels = labels)
samp_UGC07690.updateSettings({'contours': [0.682689492137086, 0.954499736103642,0.997300203936740,0.999936657516334,0.999999426696856]})
stats_UGC07690 = samp_UGC07690.getMargeStats()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07690-sigmatab-log10_a0.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_a0_sol_str='%s'%log10_a0_sol_UGC07690
  low=log10_a0_sol_UGC07690-stats_UGC07690.parWithName('log10_a0').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07690.parWithName('log10_a0').limits[i].upper- log10_a0_sol_UGC07690
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07690-sigmatab-log10_a0.txt', 'a')
  text_file.write(log10_a0_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07690-sigmatab-log10_YD.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_YD_sol_str='%s'%log10_YD_sol_UGC07690
  low=log10_YD_sol_UGC07690-stats_UGC07690.parWithName('log10_YD').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07690.parWithName('log10_YD').limits[i].upper- log10_YD_sol_UGC07690
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07690-sigmatab-log10_YD.txt', 'a')
  text_file.write(log10_YD_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07690-sigmatab-df2.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  df2_sol_str='%s'%df2_sol_UGC07690
  low=df2_sol_UGC07690-stats_UGC07690.parWithName('df2').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07690.parWithName('df2').limits[i].upper- df2_sol_UGC07690
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07690-sigmatab-df2.txt', 'a')
  text_file.write(df2_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07690-sigmatab-Dinc.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  Dinc_sol_str='%s'%Dinc_sol_UGC07690
  low=Dinc_sol_UGC07690-stats_UGC07690.parWithName('Dinc').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07690.parWithName('Dinc').limits[i].upper- Dinc_sol_UGC07690
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07690-sigmatab-Dinc.txt', 'a')
  text_file.write(Dinc_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()


full_chain_UGC07866 = pd.read_csv('/home/alefe/Cluster/MOND_outputs/chain_UGC07866.csv',usecols=[1,2,3,4])
data_array_chain_UGC07866 = numpy.array(full_chain_UGC07866)
mle_soln_UGC07866=[]

for i in range(ndim):
    mcmc_UGC07866 = numpy.percentile(data_array_chain_UGC07866[:, i], [16, 50, 84])
    q_UGC07866 = numpy.diff(mcmc_UGC07866)
    mle_soln_UGC07866.append(mcmc_UGC07866[1])

log10_a0_sol_UGC07866=mle_soln_UGC07866[0]
log10_YD_sol_UGC07866=mle_soln_UGC07866[1]
df2_sol_UGC07866=mle_soln_UGC07866[2]
Dinc_sol_UGC07866=mle_soln_UGC07866[3]

samp_UGC07866 = getdist.mcsamples.MCSamples(samples=data_array_chain_UGC07866,names = names, labels = labels)
samp_UGC07866.updateSettings({'contours': [0.682689492137086, 0.954499736103642,0.997300203936740,0.999936657516334,0.999999426696856]})
stats_UGC07866 = samp_UGC07866.getMargeStats()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07866-sigmatab-log10_a0.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_a0_sol_str='%s'%log10_a0_sol_UGC07866
  low=log10_a0_sol_UGC07866-stats_UGC07866.parWithName('log10_a0').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07866.parWithName('log10_a0').limits[i].upper- log10_a0_sol_UGC07866
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07866-sigmatab-log10_a0.txt', 'a')
  text_file.write(log10_a0_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07866-sigmatab-log10_YD.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_YD_sol_str='%s'%log10_YD_sol_UGC07866
  low=log10_YD_sol_UGC07866-stats_UGC07866.parWithName('log10_YD').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07866.parWithName('log10_YD').limits[i].upper- log10_YD_sol_UGC07866
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07866-sigmatab-log10_YD.txt', 'a')
  text_file.write(log10_YD_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07866-sigmatab-df2.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  df2_sol_str='%s'%df2_sol_UGC07866
  low=df2_sol_UGC07866-stats_UGC07866.parWithName('df2').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07866.parWithName('df2').limits[i].upper- df2_sol_UGC07866
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07866-sigmatab-df2.txt', 'a')
  text_file.write(df2_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07866-sigmatab-Dinc.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  Dinc_sol_str='%s'%Dinc_sol_UGC07866
  low=Dinc_sol_UGC07866-stats_UGC07866.parWithName('Dinc').limits[i].lower
  low_str='%s'%low
  up=stats_UGC07866.parWithName('Dinc').limits[i].upper- Dinc_sol_UGC07866
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC07866-sigmatab-Dinc.txt', 'a')
  text_file.write(Dinc_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()


full_chain_UGC08286 = pd.read_csv('/home/alefe/Cluster/MOND_outputs/chain_UGC08286.csv',usecols=[1,2,3,4])
data_array_chain_UGC08286 = numpy.array(full_chain_UGC08286)
mle_soln_UGC08286=[]

for i in range(ndim):
    mcmc_UGC08286 = numpy.percentile(data_array_chain_UGC08286[:, i], [16, 50, 84])
    q_UGC08286 = numpy.diff(mcmc_UGC08286)
    mle_soln_UGC08286.append(mcmc_UGC08286[1])

log10_a0_sol_UGC08286=mle_soln_UGC08286[0]
log10_YD_sol_UGC08286=mle_soln_UGC08286[1]
df2_sol_UGC08286=mle_soln_UGC08286[2]
Dinc_sol_UGC08286=mle_soln_UGC08286[3]

samp_UGC08286 = getdist.mcsamples.MCSamples(samples=data_array_chain_UGC08286,names = names, labels = labels)
samp_UGC08286.updateSettings({'contours': [0.682689492137086, 0.954499736103642,0.997300203936740,0.999936657516334,0.999999426696856]})
stats_UGC08286 = samp_UGC08286.getMargeStats()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC08286-sigmatab-log10_a0.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_a0_sol_str='%s'%log10_a0_sol_UGC08286
  low=log10_a0_sol_UGC08286-stats_UGC08286.parWithName('log10_a0').limits[i].lower
  low_str='%s'%low
  up=stats_UGC08286.parWithName('log10_a0').limits[i].upper- log10_a0_sol_UGC08286
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC08286-sigmatab-log10_a0.txt', 'a')
  text_file.write(log10_a0_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC08286-sigmatab-log10_YD.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_YD_sol_str='%s'%log10_YD_sol_UGC08286
  low=log10_YD_sol_UGC08286-stats_UGC08286.parWithName('log10_YD').limits[i].lower
  low_str='%s'%low
  up=stats_UGC08286.parWithName('log10_YD').limits[i].upper- log10_YD_sol_UGC08286
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC08286-sigmatab-log10_YD.txt', 'a')
  text_file.write(log10_YD_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC08286-sigmatab-df2.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  df2_sol_str='%s'%df2_sol_UGC08286
  low=df2_sol_UGC08286-stats_UGC08286.parWithName('df2').limits[i].lower
  low_str='%s'%low
  up=stats_UGC08286.parWithName('df2').limits[i].upper- df2_sol_UGC08286
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC08286-sigmatab-df2.txt', 'a')
  text_file.write(df2_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC08286-sigmatab-Dinc.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  Dinc_sol_str='%s'%Dinc_sol_UGC08286
  low=Dinc_sol_UGC08286-stats_UGC08286.parWithName('Dinc').limits[i].lower
  low_str='%s'%low
  up=stats_UGC08286.parWithName('Dinc').limits[i].upper- Dinc_sol_UGC08286
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC08286-sigmatab-Dinc.txt', 'a')
  text_file.write(Dinc_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()


full_chain_UGC08490 = pd.read_csv('/home/alefe/Cluster/MOND_outputs/chain_UGC08490.csv',usecols=[1,2,3,4])
data_array_chain_UGC08490 = numpy.array(full_chain_UGC08490)
mle_soln_UGC08490=[]

for i in range(ndim):
    mcmc_UGC08490 = numpy.percentile(data_array_chain_UGC08490[:, i], [16, 50, 84])
    q_UGC08490 = numpy.diff(mcmc_UGC08490)
    mle_soln_UGC08490.append(mcmc_UGC08490[1])

log10_a0_sol_UGC08490=mle_soln_UGC08490[0]
log10_YD_sol_UGC08490=mle_soln_UGC08490[1]
df2_sol_UGC08490=mle_soln_UGC08490[2]
Dinc_sol_UGC08490=mle_soln_UGC08490[3]

samp_UGC08490 = getdist.mcsamples.MCSamples(samples=data_array_chain_UGC08490,names = names, labels = labels)
samp_UGC08490.updateSettings({'contours': [0.682689492137086, 0.954499736103642,0.997300203936740,0.999936657516334,0.999999426696856]})
stats_UGC08490 = samp_UGC08490.getMargeStats()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC08490-sigmatab-log10_a0.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_a0_sol_str='%s'%log10_a0_sol_UGC08490
  low=log10_a0_sol_UGC08490-stats_UGC08490.parWithName('log10_a0').limits[i].lower
  low_str='%s'%low
  up=stats_UGC08490.parWithName('log10_a0').limits[i].upper- log10_a0_sol_UGC08490
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC08490-sigmatab-log10_a0.txt', 'a')
  text_file.write(log10_a0_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC08490-sigmatab-log10_YD.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_YD_sol_str='%s'%log10_YD_sol_UGC08490
  low=log10_YD_sol_UGC08490-stats_UGC08490.parWithName('log10_YD').limits[i].lower
  low_str='%s'%low
  up=stats_UGC08490.parWithName('log10_YD').limits[i].upper- log10_YD_sol_UGC08490
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC08490-sigmatab-log10_YD.txt', 'a')
  text_file.write(log10_YD_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC08490-sigmatab-df2.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  df2_sol_str='%s'%df2_sol_UGC08490
  low=df2_sol_UGC08490-stats_UGC08490.parWithName('df2').limits[i].lower
  low_str='%s'%low
  up=stats_UGC08490.parWithName('df2').limits[i].upper- df2_sol_UGC08490
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC08490-sigmatab-df2.txt', 'a')
  text_file.write(df2_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC08490-sigmatab-Dinc.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  Dinc_sol_str='%s'%Dinc_sol_UGC08490
  low=Dinc_sol_UGC08490-stats_UGC08490.parWithName('Dinc').limits[i].lower
  low_str='%s'%low
  up=stats_UGC08490.parWithName('Dinc').limits[i].upper- Dinc_sol_UGC08490
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC08490-sigmatab-Dinc.txt', 'a')
  text_file.write(Dinc_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()


full_chain_UGC08550 = pd.read_csv('/home/alefe/Cluster/MOND_outputs/chain_UGC08550.csv',usecols=[1,2,3,4])
data_array_chain_UGC08550 = numpy.array(full_chain_UGC08550)
mle_soln_UGC08550=[]

for i in range(ndim):
    mcmc_UGC08550 = numpy.percentile(data_array_chain_UGC08550[:, i], [16, 50, 84])
    q_UGC08550 = numpy.diff(mcmc_UGC08550)
    mle_soln_UGC08550.append(mcmc_UGC08550[1])

log10_a0_sol_UGC08550=mle_soln_UGC08550[0]
log10_YD_sol_UGC08550=mle_soln_UGC08550[1]
df2_sol_UGC08550=mle_soln_UGC08550[2]
Dinc_sol_UGC08550=mle_soln_UGC08550[3]

samp_UGC08550 = getdist.mcsamples.MCSamples(samples=data_array_chain_UGC08550,names = names, labels = labels)
samp_UGC08550.updateSettings({'contours': [0.682689492137086, 0.954499736103642,0.997300203936740,0.999936657516334,0.999999426696856]})
stats_UGC08550 = samp_UGC08550.getMargeStats()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC08550-sigmatab-log10_a0.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_a0_sol_str='%s'%log10_a0_sol_UGC08550
  low=log10_a0_sol_UGC08550-stats_UGC08550.parWithName('log10_a0').limits[i].lower
  low_str='%s'%low
  up=stats_UGC08550.parWithName('log10_a0').limits[i].upper- log10_a0_sol_UGC08550
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC08550-sigmatab-log10_a0.txt', 'a')
  text_file.write(log10_a0_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC08550-sigmatab-log10_YD.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_YD_sol_str='%s'%log10_YD_sol_UGC08550
  low=log10_YD_sol_UGC08550-stats_UGC08550.parWithName('log10_YD').limits[i].lower
  low_str='%s'%low
  up=stats_UGC08550.parWithName('log10_YD').limits[i].upper- log10_YD_sol_UGC08550
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC08550-sigmatab-log10_YD.txt', 'a')
  text_file.write(log10_YD_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC08550-sigmatab-df2.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  df2_sol_str='%s'%df2_sol_UGC08550
  low=df2_sol_UGC08550-stats_UGC08550.parWithName('df2').limits[i].lower
  low_str='%s'%low
  up=stats_UGC08550.parWithName('df2').limits[i].upper- df2_sol_UGC08550
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC08550-sigmatab-df2.txt', 'a')
  text_file.write(df2_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC08550-sigmatab-Dinc.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  Dinc_sol_str='%s'%Dinc_sol_UGC08550
  low=Dinc_sol_UGC08550-stats_UGC08550.parWithName('Dinc').limits[i].lower
  low_str='%s'%low
  up=stats_UGC08550.parWithName('Dinc').limits[i].upper- Dinc_sol_UGC08550
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC08550-sigmatab-Dinc.txt', 'a')
  text_file.write(Dinc_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()


full_chain_UGC08699 = pd.read_csv('/home/alefe/Cluster/MOND_outputs/chain_UGC08699.csv',usecols=[1,2,3,4,5])
data_array_chain_UGC08699 = numpy.array(full_chain_UGC08699)
mle_soln_UGC08699=[]

for i in range(ndim_with_bulge):
    mcmc_UGC08699 = numpy.percentile(data_array_chain_UGC08699[:, i], [16, 50, 84])
    q_UGC08699 = numpy.diff(mcmc_UGC08699)
    mle_soln_UGC08699.append(mcmc_UGC08699[1])

log10_a0_sol_UGC08699=mle_soln_UGC08699[0]
log10_YD_sol_UGC08699=mle_soln_UGC08699[1]
log10_YB_sol_UGC08699=mle_soln_UGC08699[2]
df2_sol_UGC08699=mle_soln_UGC08699[3]
Dinc_sol_UGC08699=mle_soln_UGC08699[4]

samp_UGC08699 = getdist.mcsamples.MCSamples(samples=data_array_chain_UGC08699,names = names_with_bulge, labels = labels_with_bulge)
samp_UGC08699.updateSettings({'contours': [0.682689492137086, 0.954499736103642,0.997300203936740,0.999936657516334,0.999999426696856]})
stats_UGC08699 = samp_UGC08699.getMargeStats()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC08699-sigmatab-log10_a0.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_a0_sol_str='%s'%log10_a0_sol_UGC08699
  low=log10_a0_sol_UGC08699-stats_UGC08699.parWithName('log10_a0').limits[i].lower
  low_str='%s'%low
  up=stats_UGC08699.parWithName('log10_a0').limits[i].upper- log10_a0_sol_UGC08699
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC08699-sigmatab-log10_a0.txt', 'a')
  text_file.write(log10_a0_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC08699-sigmatab-log10_YD.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_YD_sol_str='%s'%log10_YD_sol_UGC08699
  low=log10_YD_sol_UGC08699-stats_UGC08699.parWithName('log10_YD').limits[i].lower
  low_str='%s'%low
  up=stats_UGC08699.parWithName('log10_YD').limits[i].upper- log10_YD_sol_UGC08699
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC08699-sigmatab-log10_YD.txt', 'a')
  text_file.write(log10_YD_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC08699-sigmatab-log10_YB.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_YB_sol_str='%s'%log10_YB_sol_UGC08699
  low=log10_YB_sol_UGC08699-stats_UGC08699.parWithName('log10_YB').limits[i].lower
  low_str='%s'%low
  up=stats_UGC08699.parWithName('log10_YB').limits[i].upper- log10_YB_sol_UGC08699
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC08699-sigmatab-log10_YB.txt', 'a')
  text_file.write(log10_YB_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC08699-sigmatab-df2.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  df2_sol_str='%s'%df2_sol_UGC08699
  low=df2_sol_UGC08699-stats_UGC08699.parWithName('df2').limits[i].lower
  low_str='%s'%low
  up=stats_UGC08699.parWithName('df2').limits[i].upper- df2_sol_UGC08699
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC08699-sigmatab-df2.txt', 'a')
  text_file.write(df2_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC08699-sigmatab-Dinc.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  Dinc_sol_str='%s'%Dinc_sol_UGC08699
  low=Dinc_sol_UGC08699-stats_UGC08699.parWithName('Dinc').limits[i].lower
  low_str='%s'%low
  up=stats_UGC08699.parWithName('Dinc').limits[i].upper- Dinc_sol_UGC08699
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC08699-sigmatab-Dinc.txt', 'a')
  text_file.write(Dinc_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()



full_chain_UGC08837 = pd.read_csv('/home/alefe/Cluster/MOND_outputs/chain_UGC08837.csv',usecols=[1,2,3,4])
data_array_chain_UGC08837 = numpy.array(full_chain_UGC08837)
mle_soln_UGC08837=[]

for i in range(ndim):
    mcmc_UGC08837 = numpy.percentile(data_array_chain_UGC08837[:, i], [16, 50, 84])
    q_UGC08837 = numpy.diff(mcmc_UGC08837)
    mle_soln_UGC08837.append(mcmc_UGC08837[1])

log10_a0_sol_UGC08837=mle_soln_UGC08837[0]
log10_YD_sol_UGC08837=mle_soln_UGC08837[1]
df2_sol_UGC08837=mle_soln_UGC08837[2]
Dinc_sol_UGC08837=mle_soln_UGC08837[3]

samp_UGC08837 = getdist.mcsamples.MCSamples(samples=data_array_chain_UGC08837,names = names, labels = labels)
samp_UGC08837.updateSettings({'contours': [0.682689492137086, 0.954499736103642,0.997300203936740,0.999936657516334,0.999999426696856]})
stats_UGC08837 = samp_UGC08837.getMargeStats()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC08837-sigmatab-log10_a0.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_a0_sol_str='%s'%log10_a0_sol_UGC08837
  low=log10_a0_sol_UGC08837-stats_UGC08837.parWithName('log10_a0').limits[i].lower
  low_str='%s'%low
  up=stats_UGC08837.parWithName('log10_a0').limits[i].upper- log10_a0_sol_UGC08837
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC08837-sigmatab-log10_a0.txt', 'a')
  text_file.write(log10_a0_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC08837-sigmatab-log10_YD.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_YD_sol_str='%s'%log10_YD_sol_UGC08837
  low=log10_YD_sol_UGC08837-stats_UGC08837.parWithName('log10_YD').limits[i].lower
  low_str='%s'%low
  up=stats_UGC08837.parWithName('log10_YD').limits[i].upper- log10_YD_sol_UGC08837
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC08837-sigmatab-log10_YD.txt', 'a')
  text_file.write(log10_YD_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC08837-sigmatab-df2.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  df2_sol_str='%s'%df2_sol_UGC08837
  low=df2_sol_UGC08837-stats_UGC08837.parWithName('df2').limits[i].lower
  low_str='%s'%low
  up=stats_UGC08837.parWithName('df2').limits[i].upper- df2_sol_UGC08837
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC08837-sigmatab-df2.txt', 'a')
  text_file.write(df2_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC08837-sigmatab-Dinc.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  Dinc_sol_str='%s'%Dinc_sol_UGC08837
  low=Dinc_sol_UGC08837-stats_UGC08837.parWithName('Dinc').limits[i].lower
  low_str='%s'%low
  up=stats_UGC08837.parWithName('Dinc').limits[i].upper- Dinc_sol_UGC08837
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC08837-sigmatab-Dinc.txt', 'a')
  text_file.write(Dinc_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()


full_chain_UGC09037 = pd.read_csv('/home/alefe/Cluster/MOND_outputs/chain_UGC09037.csv',usecols=[1,2,3,4])
data_array_chain_UGC09037 = numpy.array(full_chain_UGC09037)
mle_soln_UGC09037=[]

for i in range(ndim):
    mcmc_UGC09037 = numpy.percentile(data_array_chain_UGC09037[:, i], [16, 50, 84])
    q_UGC09037 = numpy.diff(mcmc_UGC09037)
    mle_soln_UGC09037.append(mcmc_UGC09037[1])

log10_a0_sol_UGC09037=mle_soln_UGC09037[0]
log10_YD_sol_UGC09037=mle_soln_UGC09037[1]
df2_sol_UGC09037=mle_soln_UGC09037[2]
Dinc_sol_UGC09037=mle_soln_UGC09037[3]

samp_UGC09037 = getdist.mcsamples.MCSamples(samples=data_array_chain_UGC09037,names = names, labels = labels)
samp_UGC09037.updateSettings({'contours': [0.682689492137086, 0.954499736103642,0.997300203936740,0.999936657516334,0.999999426696856]})
stats_UGC09037 = samp_UGC09037.getMargeStats()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC09037-sigmatab-log10_a0.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_a0_sol_str='%s'%log10_a0_sol_UGC09037
  low=log10_a0_sol_UGC09037-stats_UGC09037.parWithName('log10_a0').limits[i].lower
  low_str='%s'%low
  up=stats_UGC09037.parWithName('log10_a0').limits[i].upper- log10_a0_sol_UGC09037
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC09037-sigmatab-log10_a0.txt', 'a')
  text_file.write(log10_a0_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC09037-sigmatab-log10_YD.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_YD_sol_str='%s'%log10_YD_sol_UGC09037
  low=log10_YD_sol_UGC09037-stats_UGC09037.parWithName('log10_YD').limits[i].lower
  low_str='%s'%low
  up=stats_UGC09037.parWithName('log10_YD').limits[i].upper- log10_YD_sol_UGC09037
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC09037-sigmatab-log10_YD.txt', 'a')
  text_file.write(log10_YD_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC09037-sigmatab-df2.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  df2_sol_str='%s'%df2_sol_UGC09037
  low=df2_sol_UGC09037-stats_UGC09037.parWithName('df2').limits[i].lower
  low_str='%s'%low
  up=stats_UGC09037.parWithName('df2').limits[i].upper- df2_sol_UGC09037
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC09037-sigmatab-df2.txt', 'a')
  text_file.write(df2_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC09037-sigmatab-Dinc.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  Dinc_sol_str='%s'%Dinc_sol_UGC09037
  low=Dinc_sol_UGC09037-stats_UGC09037.parWithName('Dinc').limits[i].lower
  low_str='%s'%low
  up=stats_UGC09037.parWithName('Dinc').limits[i].upper- Dinc_sol_UGC09037
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC09037-sigmatab-Dinc.txt', 'a')
  text_file.write(Dinc_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()


full_chain_UGC09133 = pd.read_csv('/home/alefe/Cluster/MOND_outputs/chain_UGC09133.csv',usecols=[1,2,3,4,5])
data_array_chain_UGC09133 = numpy.array(full_chain_UGC09133)
mle_soln_UGC09133=[]

for i in range(ndim_with_bulge):
    mcmc_UGC09133 = numpy.percentile(data_array_chain_UGC09133[:, i], [16, 50, 84])
    q_UGC09133 = numpy.diff(mcmc_UGC09133)
    mle_soln_UGC09133.append(mcmc_UGC09133[1])

log10_a0_sol_UGC09133=mle_soln_UGC09133[0]
log10_YD_sol_UGC09133=mle_soln_UGC09133[1]
log10_YB_sol_UGC09133=mle_soln_UGC09133[2]
df2_sol_UGC09133=mle_soln_UGC09133[3]
Dinc_sol_UGC09133=mle_soln_UGC09133[4]

samp_UGC09133 = getdist.mcsamples.MCSamples(samples=data_array_chain_UGC09133,names = names_with_bulge, labels = labels_with_bulge)
samp_UGC09133.updateSettings({'contours': [0.682689492137086, 0.954499736103642,0.997300203936740,0.999936657516334,0.999999426696856]})
stats_UGC09133 = samp_UGC09133.getMargeStats()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC09133-sigmatab-log10_a0.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_a0_sol_str='%s'%log10_a0_sol_UGC09133
  low=log10_a0_sol_UGC09133-stats_UGC09133.parWithName('log10_a0').limits[i].lower
  low_str='%s'%low
  up=stats_UGC09133.parWithName('log10_a0').limits[i].upper- log10_a0_sol_UGC09133
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC09133-sigmatab-log10_a0.txt', 'a')
  text_file.write(log10_a0_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC09133-sigmatab-log10_YD.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_YD_sol_str='%s'%log10_YD_sol_UGC09133
  low=log10_YD_sol_UGC09133-stats_UGC09133.parWithName('log10_YD').limits[i].lower
  low_str='%s'%low
  up=stats_UGC09133.parWithName('log10_YD').limits[i].upper- log10_YD_sol_UGC09133
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC09133-sigmatab-log10_YD.txt', 'a')
  text_file.write(log10_YD_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC09133-sigmatab-log10_YB.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_YB_sol_str='%s'%log10_YB_sol_UGC09133
  low=log10_YB_sol_UGC09133-stats_UGC09133.parWithName('log10_YB').limits[i].lower
  low_str='%s'%low
  up=stats_UGC09133.parWithName('log10_YB').limits[i].upper- log10_YB_sol_UGC09133
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC09133-sigmatab-log10_YB.txt', 'a')
  text_file.write(log10_YB_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC09133-sigmatab-df2.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  df2_sol_str='%s'%df2_sol_UGC09133
  low=df2_sol_UGC09133-stats_UGC09133.parWithName('df2').limits[i].lower
  low_str='%s'%low
  up=stats_UGC09133.parWithName('df2').limits[i].upper- df2_sol_UGC09133
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC09133-sigmatab-df2.txt', 'a')
  text_file.write(df2_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC09133-sigmatab-Dinc.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  Dinc_sol_str='%s'%Dinc_sol_UGC09133
  low=Dinc_sol_UGC09133-stats_UGC09133.parWithName('Dinc').limits[i].lower
  low_str='%s'%low
  up=stats_UGC09133.parWithName('Dinc').limits[i].upper- Dinc_sol_UGC09133
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC09133-sigmatab-Dinc.txt', 'a')
  text_file.write(Dinc_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()



full_chain_UGC09992 = pd.read_csv('/home/alefe/Cluster/MOND_outputs/chain_UGC09992.csv',usecols=[1,2,3,4])
data_array_chain_UGC09992 = numpy.array(full_chain_UGC09992)
mle_soln_UGC09992=[]

for i in range(ndim):
    mcmc_UGC09992 = numpy.percentile(data_array_chain_UGC09992[:, i], [16, 50, 84])
    q_UGC09992 = numpy.diff(mcmc_UGC09992)
    mle_soln_UGC09992.append(mcmc_UGC09992[1])

log10_a0_sol_UGC09992=mle_soln_UGC09992[0]
log10_YD_sol_UGC09992=mle_soln_UGC09992[1]
df2_sol_UGC09992=mle_soln_UGC09992[2]
Dinc_sol_UGC09992=mle_soln_UGC09992[3]

samp_UGC09992 = getdist.mcsamples.MCSamples(samples=data_array_chain_UGC09992,names = names, labels = labels)
samp_UGC09992.updateSettings({'contours': [0.682689492137086, 0.954499736103642,0.997300203936740,0.999936657516334,0.999999426696856]})
stats_UGC09992 = samp_UGC09992.getMargeStats()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC09992-sigmatab-log10_a0.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_a0_sol_str='%s'%log10_a0_sol_UGC09992
  low=log10_a0_sol_UGC09992-stats_UGC09992.parWithName('log10_a0').limits[i].lower
  low_str='%s'%low
  up=stats_UGC09992.parWithName('log10_a0').limits[i].upper- log10_a0_sol_UGC09992
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC09992-sigmatab-log10_a0.txt', 'a')
  text_file.write(log10_a0_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC09992-sigmatab-log10_YD.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_YD_sol_str='%s'%log10_YD_sol_UGC09992
  low=log10_YD_sol_UGC09992-stats_UGC09992.parWithName('log10_YD').limits[i].lower
  low_str='%s'%low
  up=stats_UGC09992.parWithName('log10_YD').limits[i].upper- log10_YD_sol_UGC09992
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC09992-sigmatab-log10_YD.txt', 'a')
  text_file.write(log10_YD_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC09992-sigmatab-df2.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  df2_sol_str='%s'%df2_sol_UGC09992
  low=df2_sol_UGC09992-stats_UGC09992.parWithName('df2').limits[i].lower
  low_str='%s'%low
  up=stats_UGC09992.parWithName('df2').limits[i].upper- df2_sol_UGC09992
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC09992-sigmatab-df2.txt', 'a')
  text_file.write(df2_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC09992-sigmatab-Dinc.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  Dinc_sol_str='%s'%Dinc_sol_UGC09992
  low=Dinc_sol_UGC09992-stats_UGC09992.parWithName('Dinc').limits[i].lower
  low_str='%s'%low
  up=stats_UGC09992.parWithName('Dinc').limits[i].upper- Dinc_sol_UGC09992
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC09992-sigmatab-Dinc.txt', 'a')
  text_file.write(Dinc_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()


full_chain_UGC10310 = pd.read_csv('/home/alefe/Cluster/MOND_outputs/chain_UGC10310.csv',usecols=[1,2,3,4])
data_array_chain_UGC10310 = numpy.array(full_chain_UGC10310)
mle_soln_UGC10310=[]

for i in range(ndim):
    mcmc_UGC10310 = numpy.percentile(data_array_chain_UGC10310[:, i], [16, 50, 84])
    q_UGC10310 = numpy.diff(mcmc_UGC10310)
    mle_soln_UGC10310.append(mcmc_UGC10310[1])

log10_a0_sol_UGC10310=mle_soln_UGC10310[0]
log10_YD_sol_UGC10310=mle_soln_UGC10310[1]
df2_sol_UGC10310=mle_soln_UGC10310[2]
Dinc_sol_UGC10310=mle_soln_UGC10310[3]

samp_UGC10310 = getdist.mcsamples.MCSamples(samples=data_array_chain_UGC10310,names = names, labels = labels)
samp_UGC10310.updateSettings({'contours': [0.682689492137086, 0.954499736103642,0.997300203936740,0.999936657516334,0.999999426696856]})
stats_UGC10310 = samp_UGC10310.getMargeStats()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC10310-sigmatab-log10_a0.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_a0_sol_str='%s'%log10_a0_sol_UGC10310
  low=log10_a0_sol_UGC10310-stats_UGC10310.parWithName('log10_a0').limits[i].lower
  low_str='%s'%low
  up=stats_UGC10310.parWithName('log10_a0').limits[i].upper- log10_a0_sol_UGC10310
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC10310-sigmatab-log10_a0.txt', 'a')
  text_file.write(log10_a0_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC10310-sigmatab-log10_YD.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_YD_sol_str='%s'%log10_YD_sol_UGC10310
  low=log10_YD_sol_UGC10310-stats_UGC10310.parWithName('log10_YD').limits[i].lower
  low_str='%s'%low
  up=stats_UGC10310.parWithName('log10_YD').limits[i].upper- log10_YD_sol_UGC10310
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC10310-sigmatab-log10_YD.txt', 'a')
  text_file.write(log10_YD_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC10310-sigmatab-df2.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  df2_sol_str='%s'%df2_sol_UGC10310
  low=df2_sol_UGC10310-stats_UGC10310.parWithName('df2').limits[i].lower
  low_str='%s'%low
  up=stats_UGC10310.parWithName('df2').limits[i].upper- df2_sol_UGC10310
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC10310-sigmatab-df2.txt', 'a')
  text_file.write(df2_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC10310-sigmatab-Dinc.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  Dinc_sol_str='%s'%Dinc_sol_UGC10310
  low=Dinc_sol_UGC10310-stats_UGC10310.parWithName('Dinc').limits[i].lower
  low_str='%s'%low
  up=stats_UGC10310.parWithName('Dinc').limits[i].upper- Dinc_sol_UGC10310
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC10310-sigmatab-Dinc.txt', 'a')
  text_file.write(Dinc_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()


full_chain_UGC11455 = pd.read_csv('/home/alefe/Cluster/MOND_outputs/chain_UGC11455.csv',usecols=[1,2,3,4])
data_array_chain_UGC11455 = numpy.array(full_chain_UGC11455)
mle_soln_UGC11455=[]

for i in range(ndim):
    mcmc_UGC11455 = numpy.percentile(data_array_chain_UGC11455[:, i], [16, 50, 84])
    q_UGC11455 = numpy.diff(mcmc_UGC11455)
    mle_soln_UGC11455.append(mcmc_UGC11455[1])

log10_a0_sol_UGC11455=mle_soln_UGC11455[0]
log10_YD_sol_UGC11455=mle_soln_UGC11455[1]
df2_sol_UGC11455=mle_soln_UGC11455[2]
Dinc_sol_UGC11455=mle_soln_UGC11455[3]

samp_UGC11455 = getdist.mcsamples.MCSamples(samples=data_array_chain_UGC11455,names = names, labels = labels)
samp_UGC11455.updateSettings({'contours': [0.682689492137086, 0.954499736103642,0.997300203936740,0.999936657516334,0.999999426696856]})
stats_UGC11455 = samp_UGC11455.getMargeStats()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC11455-sigmatab-log10_a0.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_a0_sol_str='%s'%log10_a0_sol_UGC11455
  low=log10_a0_sol_UGC11455-stats_UGC11455.parWithName('log10_a0').limits[i].lower
  low_str='%s'%low
  up=stats_UGC11455.parWithName('log10_a0').limits[i].upper- log10_a0_sol_UGC11455
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC11455-sigmatab-log10_a0.txt', 'a')
  text_file.write(log10_a0_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC11455-sigmatab-log10_YD.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_YD_sol_str='%s'%log10_YD_sol_UGC11455
  low=log10_YD_sol_UGC11455-stats_UGC11455.parWithName('log10_YD').limits[i].lower
  low_str='%s'%low
  up=stats_UGC11455.parWithName('log10_YD').limits[i].upper- log10_YD_sol_UGC11455
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC11455-sigmatab-log10_YD.txt', 'a')
  text_file.write(log10_YD_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC11455-sigmatab-df2.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  df2_sol_str='%s'%df2_sol_UGC11455
  low=df2_sol_UGC11455-stats_UGC11455.parWithName('df2').limits[i].lower
  low_str='%s'%low
  up=stats_UGC11455.parWithName('df2').limits[i].upper- df2_sol_UGC11455
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC11455-sigmatab-df2.txt', 'a')
  text_file.write(df2_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC11455-sigmatab-Dinc.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  Dinc_sol_str='%s'%Dinc_sol_UGC11455
  low=Dinc_sol_UGC11455-stats_UGC11455.parWithName('Dinc').limits[i].lower
  low_str='%s'%low
  up=stats_UGC11455.parWithName('Dinc').limits[i].upper- Dinc_sol_UGC11455
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC11455-sigmatab-Dinc.txt', 'a')
  text_file.write(Dinc_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()


full_chain_UGC11557 = pd.read_csv('/home/alefe/Cluster/MOND_outputs/chain_UGC11557.csv',usecols=[1,2,3,4])
data_array_chain_UGC11557 = numpy.array(full_chain_UGC11557)
mle_soln_UGC11557=[]

for i in range(ndim):
    mcmc_UGC11557 = numpy.percentile(data_array_chain_UGC11557[:, i], [16, 50, 84])
    q_UGC11557 = numpy.diff(mcmc_UGC11557)
    mle_soln_UGC11557.append(mcmc_UGC11557[1])

log10_a0_sol_UGC11557=mle_soln_UGC11557[0]
log10_YD_sol_UGC11557=mle_soln_UGC11557[1]
df2_sol_UGC11557=mle_soln_UGC11557[2]
Dinc_sol_UGC11557=mle_soln_UGC11557[3]

samp_UGC11557 = getdist.mcsamples.MCSamples(samples=data_array_chain_UGC11557,names = names, labels = labels)
samp_UGC11557.updateSettings({'contours': [0.682689492137086, 0.954499736103642,0.997300203936740,0.999936657516334,0.999999426696856]})
stats_UGC11557 = samp_UGC11557.getMargeStats()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC11557-sigmatab-log10_a0.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_a0_sol_str='%s'%log10_a0_sol_UGC11557
  low=log10_a0_sol_UGC11557-stats_UGC11557.parWithName('log10_a0').limits[i].lower
  low_str='%s'%low
  up=stats_UGC11557.parWithName('log10_a0').limits[i].upper- log10_a0_sol_UGC11557
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC11557-sigmatab-log10_a0.txt', 'a')
  text_file.write(log10_a0_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC11557-sigmatab-log10_YD.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_YD_sol_str='%s'%log10_YD_sol_UGC11557
  low=log10_YD_sol_UGC11557-stats_UGC11557.parWithName('log10_YD').limits[i].lower
  low_str='%s'%low
  up=stats_UGC11557.parWithName('log10_YD').limits[i].upper- log10_YD_sol_UGC11557
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC11557-sigmatab-log10_YD.txt', 'a')
  text_file.write(log10_YD_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC11557-sigmatab-df2.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  df2_sol_str='%s'%df2_sol_UGC11557
  low=df2_sol_UGC11557-stats_UGC11557.parWithName('df2').limits[i].lower
  low_str='%s'%low
  up=stats_UGC11557.parWithName('df2').limits[i].upper- df2_sol_UGC11557
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC11557-sigmatab-df2.txt', 'a')
  text_file.write(df2_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC11557-sigmatab-Dinc.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  Dinc_sol_str='%s'%Dinc_sol_UGC11557
  low=Dinc_sol_UGC11557-stats_UGC11557.parWithName('Dinc').limits[i].lower
  low_str='%s'%low
  up=stats_UGC11557.parWithName('Dinc').limits[i].upper- Dinc_sol_UGC11557
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC11557-sigmatab-Dinc.txt', 'a')
  text_file.write(Dinc_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()


full_chain_UGC11820 = pd.read_csv('/home/alefe/Cluster/MOND_outputs/chain_UGC11820.csv',usecols=[1,2,3,4])
data_array_chain_UGC11820 = numpy.array(full_chain_UGC11820)
mle_soln_UGC11820=[]

for i in range(ndim):
    mcmc_UGC11820 = numpy.percentile(data_array_chain_UGC11820[:, i], [16, 50, 84])
    q_UGC11820 = numpy.diff(mcmc_UGC11820)
    mle_soln_UGC11820.append(mcmc_UGC11820[1])

log10_a0_sol_UGC11820=mle_soln_UGC11820[0]
log10_YD_sol_UGC11820=mle_soln_UGC11820[1]
df2_sol_UGC11820=mle_soln_UGC11820[2]
Dinc_sol_UGC11820=mle_soln_UGC11820[3]

samp_UGC11820 = getdist.mcsamples.MCSamples(samples=data_array_chain_UGC11820,names = names, labels = labels)
samp_UGC11820.updateSettings({'contours': [0.682689492137086, 0.954499736103642,0.997300203936740,0.999936657516334,0.999999426696856]})
stats_UGC11820 = samp_UGC11820.getMargeStats()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC11820-sigmatab-log10_a0.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_a0_sol_str='%s'%log10_a0_sol_UGC11820
  low=log10_a0_sol_UGC11820-stats_UGC11820.parWithName('log10_a0').limits[i].lower
  low_str='%s'%low
  up=stats_UGC11820.parWithName('log10_a0').limits[i].upper- log10_a0_sol_UGC11820
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC11820-sigmatab-log10_a0.txt', 'a')
  text_file.write(log10_a0_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC11820-sigmatab-log10_YD.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_YD_sol_str='%s'%log10_YD_sol_UGC11820
  low=log10_YD_sol_UGC11820-stats_UGC11820.parWithName('log10_YD').limits[i].lower
  low_str='%s'%low
  up=stats_UGC11820.parWithName('log10_YD').limits[i].upper- log10_YD_sol_UGC11820
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC11820-sigmatab-log10_YD.txt', 'a')
  text_file.write(log10_YD_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC11820-sigmatab-df2.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  df2_sol_str='%s'%df2_sol_UGC11820
  low=df2_sol_UGC11820-stats_UGC11820.parWithName('df2').limits[i].lower
  low_str='%s'%low
  up=stats_UGC11820.parWithName('df2').limits[i].upper- df2_sol_UGC11820
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC11820-sigmatab-df2.txt', 'a')
  text_file.write(df2_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC11820-sigmatab-Dinc.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  Dinc_sol_str='%s'%Dinc_sol_UGC11820
  low=Dinc_sol_UGC11820-stats_UGC11820.parWithName('Dinc').limits[i].lower
  low_str='%s'%low
  up=stats_UGC11820.parWithName('Dinc').limits[i].upper- Dinc_sol_UGC11820
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC11820-sigmatab-Dinc.txt', 'a')
  text_file.write(Dinc_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()


full_chain_UGC11914 = pd.read_csv('/home/alefe/Cluster/MOND_outputs/chain_UGC11914.csv',usecols=[1,2,3,4,5])
data_array_chain_UGC11914 = numpy.array(full_chain_UGC11914)
mle_soln_UGC11914=[]

for i in range(ndim_with_bulge):
    mcmc_UGC11914 = numpy.percentile(data_array_chain_UGC11914[:, i], [16, 50, 84])
    q_UGC11914 = numpy.diff(mcmc_UGC11914)
    mle_soln_UGC11914.append(mcmc_UGC11914[1])

log10_a0_sol_UGC11914=mle_soln_UGC11914[0]
log10_YD_sol_UGC11914=mle_soln_UGC11914[1]
log10_YB_sol_UGC11914=mle_soln_UGC11914[2]
df2_sol_UGC11914=mle_soln_UGC11914[3]
Dinc_sol_UGC11914=mle_soln_UGC11914[4]

samp_UGC11914 = getdist.mcsamples.MCSamples(samples=data_array_chain_UGC11914,names = names_with_bulge, labels = labels_with_bulge)
samp_UGC11914.updateSettings({'contours': [0.682689492137086, 0.954499736103642,0.997300203936740,0.999936657516334,0.999999426696856]})
stats_UGC11914 = samp_UGC11914.getMargeStats()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC11914-sigmatab-log10_a0.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_a0_sol_str='%s'%log10_a0_sol_UGC11914
  low=log10_a0_sol_UGC11914-stats_UGC11914.parWithName('log10_a0').limits[i].lower
  low_str='%s'%low
  up=stats_UGC11914.parWithName('log10_a0').limits[i].upper- log10_a0_sol_UGC11914
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC11914-sigmatab-log10_a0.txt', 'a')
  text_file.write(log10_a0_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC11914-sigmatab-log10_YD.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_YD_sol_str='%s'%log10_YD_sol_UGC11914
  low=log10_YD_sol_UGC11914-stats_UGC11914.parWithName('log10_YD').limits[i].lower
  low_str='%s'%low
  up=stats_UGC11914.parWithName('log10_YD').limits[i].upper- log10_YD_sol_UGC11914
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC11914-sigmatab-log10_YD.txt', 'a')
  text_file.write(log10_YD_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC11914-sigmatab-log10_YB.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_YB_sol_str='%s'%log10_YB_sol_UGC11914
  low=log10_YB_sol_UGC11914-stats_UGC11914.parWithName('log10_YB').limits[i].lower
  low_str='%s'%low
  up=stats_UGC11914.parWithName('log10_YB').limits[i].upper- log10_YB_sol_UGC11914
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC11914-sigmatab-log10_YB.txt', 'a')
  text_file.write(log10_YB_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC11914-sigmatab-df2.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  df2_sol_str='%s'%df2_sol_UGC11914
  low=df2_sol_UGC11914-stats_UGC11914.parWithName('df2').limits[i].lower
  low_str='%s'%low
  up=stats_UGC11914.parWithName('df2').limits[i].upper- df2_sol_UGC11914
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC11914-sigmatab-df2.txt', 'a')
  text_file.write(df2_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC11914-sigmatab-Dinc.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  Dinc_sol_str='%s'%Dinc_sol_UGC11914
  low=Dinc_sol_UGC11914-stats_UGC11914.parWithName('Dinc').limits[i].lower
  low_str='%s'%low
  up=stats_UGC11914.parWithName('Dinc').limits[i].upper- Dinc_sol_UGC11914
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC11914-sigmatab-Dinc.txt', 'a')
  text_file.write(Dinc_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()



full_chain_UGC12506 = pd.read_csv('/home/alefe/Cluster/MOND_outputs/chain_UGC12506.csv',usecols=[1,2,3,4])
data_array_chain_UGC12506 = numpy.array(full_chain_UGC12506)
mle_soln_UGC12506=[]

for i in range(ndim):
    mcmc_UGC12506 = numpy.percentile(data_array_chain_UGC12506[:, i], [16, 50, 84])
    q_UGC12506 = numpy.diff(mcmc_UGC12506)
    mle_soln_UGC12506.append(mcmc_UGC12506[1])

log10_a0_sol_UGC12506=mle_soln_UGC12506[0]
log10_YD_sol_UGC12506=mle_soln_UGC12506[1]
df2_sol_UGC12506=mle_soln_UGC12506[2]
Dinc_sol_UGC12506=mle_soln_UGC12506[3]

samp_UGC12506 = getdist.mcsamples.MCSamples(samples=data_array_chain_UGC12506,names = names, labels = labels)
samp_UGC12506.updateSettings({'contours': [0.682689492137086, 0.954499736103642,0.997300203936740,0.999936657516334,0.999999426696856]})
stats_UGC12506 = samp_UGC12506.getMargeStats()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC12506-sigmatab-log10_a0.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_a0_sol_str='%s'%log10_a0_sol_UGC12506
  low=log10_a0_sol_UGC12506-stats_UGC12506.parWithName('log10_a0').limits[i].lower
  low_str='%s'%low
  up=stats_UGC12506.parWithName('log10_a0').limits[i].upper- log10_a0_sol_UGC12506
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC12506-sigmatab-log10_a0.txt', 'a')
  text_file.write(log10_a0_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC12506-sigmatab-log10_YD.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_YD_sol_str='%s'%log10_YD_sol_UGC12506
  low=log10_YD_sol_UGC12506-stats_UGC12506.parWithName('log10_YD').limits[i].lower
  low_str='%s'%low
  up=stats_UGC12506.parWithName('log10_YD').limits[i].upper- log10_YD_sol_UGC12506
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC12506-sigmatab-log10_YD.txt', 'a')
  text_file.write(log10_YD_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC12506-sigmatab-df2.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  df2_sol_str='%s'%df2_sol_UGC12506
  low=df2_sol_UGC12506-stats_UGC12506.parWithName('df2').limits[i].lower
  low_str='%s'%low
  up=stats_UGC12506.parWithName('df2').limits[i].upper- df2_sol_UGC12506
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC12506-sigmatab-df2.txt', 'a')
  text_file.write(df2_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC12506-sigmatab-Dinc.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  Dinc_sol_str='%s'%Dinc_sol_UGC12506
  low=Dinc_sol_UGC12506-stats_UGC12506.parWithName('Dinc').limits[i].lower
  low_str='%s'%low
  up=stats_UGC12506.parWithName('Dinc').limits[i].upper- Dinc_sol_UGC12506
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC12506-sigmatab-Dinc.txt', 'a')
  text_file.write(Dinc_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()


full_chain_UGC12632 = pd.read_csv('/home/alefe/Cluster/MOND_outputs/chain_UGC12632.csv',usecols=[1,2,3,4])
data_array_chain_UGC12632 = numpy.array(full_chain_UGC12632)
mle_soln_UGC12632=[]

for i in range(ndim):
    mcmc_UGC12632 = numpy.percentile(data_array_chain_UGC12632[:, i], [16, 50, 84])
    q_UGC12632 = numpy.diff(mcmc_UGC12632)
    mle_soln_UGC12632.append(mcmc_UGC12632[1])

log10_a0_sol_UGC12632=mle_soln_UGC12632[0]
log10_YD_sol_UGC12632=mle_soln_UGC12632[1]
df2_sol_UGC12632=mle_soln_UGC12632[2]
Dinc_sol_UGC12632=mle_soln_UGC12632[3]

samp_UGC12632 = getdist.mcsamples.MCSamples(samples=data_array_chain_UGC12632,names = names, labels = labels)
samp_UGC12632.updateSettings({'contours': [0.682689492137086, 0.954499736103642,0.997300203936740,0.999936657516334,0.999999426696856]})
stats_UGC12632 = samp_UGC12632.getMargeStats()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC12632-sigmatab-log10_a0.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_a0_sol_str='%s'%log10_a0_sol_UGC12632
  low=log10_a0_sol_UGC12632-stats_UGC12632.parWithName('log10_a0').limits[i].lower
  low_str='%s'%low
  up=stats_UGC12632.parWithName('log10_a0').limits[i].upper- log10_a0_sol_UGC12632
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC12632-sigmatab-log10_a0.txt', 'a')
  text_file.write(log10_a0_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC12632-sigmatab-log10_YD.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  log10_YD_sol_str='%s'%log10_YD_sol_UGC12632
  low=log10_YD_sol_UGC12632-stats_UGC12632.parWithName('log10_YD').limits[i].lower
  low_str='%s'%low
  up=stats_UGC12632.parWithName('log10_YD').limits[i].upper- log10_YD_sol_UGC12632
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC12632-sigmatab-log10_YD.txt', 'a')
  text_file.write(log10_YD_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC12632-sigmatab-df2.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  df2_sol_str='%s'%df2_sol_UGC12632
  low=df2_sol_UGC12632-stats_UGC12632.parWithName('df2').limits[i].lower
  low_str='%s'%low
  up=stats_UGC12632.parWithName('df2').limits[i].upper- df2_sol_UGC12632
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC12632-sigmatab-df2.txt', 'a')
  text_file.write(df2_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()

text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC12632-sigmatab-Dinc.txt', 'w')
text_file.write('max s- s+   -  rows are 1,2,3,4,5 sigmas\n')
text_file.close()

for i in range(5):
  Dinc_sol_str='%s'%Dinc_sol_UGC12632
  low=Dinc_sol_UGC12632-stats_UGC12632.parWithName('Dinc').limits[i].lower
  low_str='%s'%low
  up=stats_UGC12632.parWithName('Dinc').limits[i].upper- Dinc_sol_UGC12632
  up_str='%s'%up
  text_file = open('/home/alefe/Cluster/MOND_GetDist_corner/UGC12632-sigmatab-Dinc.txt', 'a')
  text_file.write(Dinc_sol_str + ',' + low_str + ',' + up_str + '\n')
  text_file.close()
