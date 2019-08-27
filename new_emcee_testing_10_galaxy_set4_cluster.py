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


#importing data
#importing data
filename_NGC7331='NGC7331_rotmod.dat'
filename_NGC7793='NGC7793_rotmod.dat'
filename_NGC7814='NGC7814_rotmod.dat'
filename_UGC02259='UGC02259_rotmod.dat'
filename_UGC03546='UGC03546_rotmod.dat'
filename_UGC06446='UGC06446_rotmod.dat'
filename_UGC06930='UGC06930_rotmod.dat'
filename_UGC06983='UGC06983_rotmod.dat'
filename_UGC07261='UGC07261_rotmod.dat'
filename_UGC07690='UGC07690_rotmod.dat'

data_NGC7331=pd.read_table(filename_NGC7331,header=None)
data_NGC7793=pd.read_table(filename_NGC7793,header=None)
data_NGC7814=pd.read_table(filename_NGC7814,header=None)
data_UGC02259=pd.read_table(filename_UGC02259,header=None)
data_UGC03546=pd.read_table(filename_UGC03546,header=None)
data_UGC06446=pd.read_table(filename_UGC06446,header=None)
data_UGC06930=pd.read_table(filename_UGC06930,header=None)
data_UGC06983=pd.read_table(filename_UGC06983,header=None)
data_UGC07261=pd.read_table(filename_UGC07261,header=None)
data_UGC07690=pd.read_table(filename_UGC07690,header=None)

data_array_NGC7331=np.array(data_NGC7331)
data_array_NGC7793=np.array(data_NGC7793)
data_array_NGC7814=np.array(data_NGC7814)
data_array_UGC02259=np.array(data_UGC02259)
data_array_UGC03546=np.array(data_UGC03546)
data_array_UGC06446=np.array(data_UGC06446)
data_array_UGC06930=np.array(data_UGC06930)
data_array_UGC06983=np.array(data_UGC06983)
data_array_UGC07261=np.array(data_UGC07261)
data_array_UGC07690=np.array(data_UGC07690)

radius_NGC7331=data_array_NGC7331[:,0]
radius_NGC7793=data_array_NGC7793[:,0]
radius_NGC7814=data_array_NGC7814[:,0]
radius_UGC02259=data_array_UGC02259[:,0]
radius_UGC03546=data_array_UGC03546[:,0]
radius_UGC06446=data_array_UGC06446[:,0]
radius_UGC06930=data_array_UGC06930[:,0]
radius_UGC06983=data_array_UGC06983[:,0]
radius_UGC07261=data_array_UGC07261[:,0]
radius_UGC07690=data_array_UGC07690[:,0]

velocity_NGC7331=data_array_NGC7331[:,1]
velocity_NGC7793=data_array_NGC7793[:,1]
velocity_NGC7814=data_array_NGC7814[:,1]
velocity_UGC02259=data_array_UGC02259[:,1]
velocity_UGC03546=data_array_UGC03546[:,1]
velocity_UGC06446=data_array_UGC06446[:,1]
velocity_UGC06930=data_array_UGC06930[:,1]
velocity_UGC06983=data_array_UGC06983[:,1]
velocity_UGC07261=data_array_UGC07261[:,1]
velocity_UGC07690=data_array_UGC07690[:,1]

error_velocity_NGC7331=data_array_NGC7331[:,2]
error_velocity_NGC7793=data_array_NGC7793[:,2]
error_velocity_NGC7814=data_array_NGC7814[:,2]
error_velocity_UGC02259=data_array_UGC02259[:,2]
error_velocity_UGC03546=data_array_UGC03546[:,2]
error_velocity_UGC06446=data_array_UGC06446[:,2]
error_velocity_UGC06930=data_array_UGC06930[:,2]
error_velocity_UGC06983=data_array_UGC06983[:,2]
error_velocity_UGC07261=data_array_UGC07261[:,2]
error_velocity_UGC07690=data_array_UGC07690[:,2]

gas_NGC7331=data_array_NGC7331[:,3]
gas_NGC7793=data_array_NGC7793[:,3]
gas_NGC7814=data_array_NGC7814[:,3]
gas_UGC02259=data_array_UGC02259[:,3]
gas_UGC03546=data_array_UGC03546[:,3]
gas_UGC06446=data_array_UGC06446[:,3]
gas_UGC06930=data_array_UGC06930[:,3]
gas_UGC06983=data_array_UGC06983[:,3]
gas_UGC07261=data_array_UGC07261[:,3]
gas_UGC07690=data_array_UGC07690[:,3]

disk_NGC7331=data_array_NGC7331[:,4]
disk_NGC7793=data_array_NGC7793[:,4]
disk_NGC7814=data_array_NGC7814[:,4]
disk_UGC02259=data_array_UGC02259[:,4]
disk_UGC03546=data_array_UGC03546[:,4]
disk_UGC06446=data_array_UGC06446[:,4]
disk_UGC06930=data_array_UGC06930[:,4]
disk_UGC06983=data_array_UGC06983[:,4]
disk_UGC07261=data_array_UGC07261[:,4]
disk_UGC07690=data_array_UGC07690[:,4]


bulge_NGC7331=data_array_NGC7331[:,5]
bulge_NGC7814=data_array_NGC7814[:,5]
bulge_UGC03546=data_array_UGC03546[:,5]


#Rotation curves of gas and disk

vdisk_NGC7331=InterpolatedUnivariateSpline(radius_NGC7331,disk_NGC7331)

vgas_NGC7331=InterpolatedUnivariateSpline(radius_NGC7331,gas_NGC7331)

vbulge_NGC7331=InterpolatedUnivariateSpline(radius_NGC7331,bulge_NGC7331)

vdisk_NGC7793=InterpolatedUnivariateSpline(radius_NGC7793,disk_NGC7793)

vgas_NGC7793=InterpolatedUnivariateSpline(radius_NGC7793,gas_NGC7793)

#no bulge

vdisk_NGC7814=InterpolatedUnivariateSpline(radius_NGC7814,disk_NGC7814)

vgas_NGC7814=InterpolatedUnivariateSpline(radius_NGC7814,gas_NGC7814)

vbulge_NGC7814=InterpolatedUnivariateSpline(radius_NGC7814,bulge_NGC7814)

vdisk_UGC02259=InterpolatedUnivariateSpline(radius_UGC02259,disk_UGC02259)

vgas_UGC02259=InterpolatedUnivariateSpline(radius_UGC02259,gas_UGC02259)

#no bulge

vdisk_UGC03546=InterpolatedUnivariateSpline(radius_UGC03546,disk_UGC03546)

vgas_UGC03546=InterpolatedUnivariateSpline(radius_UGC03546,gas_UGC03546)

vbulge_UGC03546=InterpolatedUnivariateSpline(radius_UGC03546,bulge_UGC03546)

vdisk_UGC06446=InterpolatedUnivariateSpline(radius_UGC06446,disk_UGC06446)

vgas_UGC06446=InterpolatedUnivariateSpline(radius_UGC06446,gas_UGC06446)

#no bulge

vdisk_UGC06930=InterpolatedUnivariateSpline(radius_UGC06930,disk_UGC06930)

vgas_UGC06930=InterpolatedUnivariateSpline(radius_UGC06930,gas_UGC06930)

#no bulge

vdisk_UGC06983=InterpolatedUnivariateSpline(radius_UGC06983,disk_UGC06983)

vgas_UGC06983=InterpolatedUnivariateSpline(radius_UGC06983,gas_UGC06983)

#no bulge

vdisk_UGC07261=InterpolatedUnivariateSpline(radius_UGC07261,disk_UGC07261)

vgas_UGC07261=InterpolatedUnivariateSpline(radius_UGC07261,gas_UGC07261)

#no bulge

vdisk_UGC07690=InterpolatedUnivariateSpline(radius_UGC07690,disk_UGC07690)

vgas_UGC07690=InterpolatedUnivariateSpline(radius_UGC07690,gas_UGC07690)

#no bulge

#Setting some important constants
import math

kpc = 3.08568*10**16 #To convert from kpc to km
G0 = 4.51824*10**-39 # in kpc^3/Msun*s^2 units
H0 = 2.33336*10**-18 # in s**-1 units
rho_c = (3*H0**2)/(8*math.pi*4.51824*10**-39)# critical density Msun/kpc^3
h=0.671# h

#The rotation curves formula:
#NGC7331

rs_NGC7331=lambda m200_NGC7331:28.8*(m200_NGC7331/(10**12)/h)**0.43

c_NGC7331=lambda m200_NGC7331:(10**0.905)*(m200_NGC7331/(10**12)/h)**(-0.101)

rho_s_NGC7331=lambda m200_NGC7331:(200/3)*(rho_c*c_NGC7331(m200_NGC7331)**3)/(np.log(1+c_NGC7331(m200_NGC7331))-c_NGC7331(m200_NGC7331)/(1+c_NGC7331(m200_NGC7331)))

mnfw_NGC7331=lambda r,m200_NGC7331:4*math.pi*(rs_NGC7331(m200_NGC7331)**3)*rho_s_NGC7331(m200_NGC7331)*(-(r/(r+rs_NGC7331(m200_NGC7331)))-np.log(rs_NGC7331(m200_NGC7331))+np.log(r+rs_NGC7331(m200_NGC7331)))

vvnfwinmg_NGC7331= lambda r,m200_NGC7331,b,l:((-2*math.pi*G0*(kpc**2)*b*rho_s_NGC7331(m200_NGC7331)*rs_NGC7331(m200_NGC7331)**3)/r)*(2*r/(r+rs_NGC7331(m200_NGC7331))+np.exp((r+rs_NGC7331(m200_NGC7331))/l)*(r/l-1)*expi(-(r+rs_NGC7331(m200_NGC7331))/l)
  +np.exp(-(r+rs_NGC7331(m200_NGC7331))/l)*(1+r/l)*(np.exp(2*rs_NGC7331(m200_NGC7331)/l)*expi(-rs_NGC7331(m200_NGC7331)/l)+expi(rs_NGC7331(m200_NGC7331)/l)-expi((r+rs_NGC7331(m200_NGC7331))/l)))

vvnfw_NGC7331= lambda r,m200_NGC7331:G0*(kpc**2)*mnfw_NGC7331(r,m200_NGC7331)/r

vv_NGC7331=lambda r,Y_NGC7331,Yb_NGC7331,m200_NGC7331,b,l: vgas_NGC7331(r)**2+G0*(kpc**2)*mnfw_NGC7331(r,m200_NGC7331)/r+vvnfwinmg_NGC7331(r,m200_NGC7331,b,l)+Y_NGC7331*vdisk_NGC7331(r)**2+Yb_NGC7331*vbulge_NGC7331(r)**2

#NGC7793

rs_NGC7793=lambda m200_NGC7793:28.8*(m200_NGC7793/(10**12)/h)**0.43

c_NGC7793=lambda m200_NGC7793:(10**0.905)*(m200_NGC7793/(10**12)/h)**(-0.101)

rho_s_NGC7793=lambda m200_NGC7793:(200/3)*(rho_c*c_NGC7793(m200_NGC7793)**3)/(np.log(1+c_NGC7793(m200_NGC7793))-c_NGC7793(m200_NGC7793)/(1+c_NGC7793(m200_NGC7793)))

mnfw_NGC7793=lambda r,m200_NGC7793:4*math.pi*(rs_NGC7793(m200_NGC7793)**3)*rho_s_NGC7793(m200_NGC7793)*(-(r/(r+rs_NGC7793(m200_NGC7793)))-np.log(rs_NGC7793(m200_NGC7793))+np.log(r+rs_NGC7793(m200_NGC7793)))

vvnfwinmg_NGC7793= lambda r,m200_NGC7793,b,l:((-2*math.pi*G0*(kpc**2)*b*rho_s_NGC7793(m200_NGC7793)*rs_NGC7793(m200_NGC7793)**3)/r)*(2*r/(r+rs_NGC7793(m200_NGC7793))+np.exp((r+rs_NGC7793(m200_NGC7793))/l)*(r/l-1)*expi(-(r+rs_NGC7793(m200_NGC7793))/l)
  +np.exp(-(r+rs_NGC7793(m200_NGC7793))/l)*(1+r/l)*(np.exp(2*rs_NGC7793(m200_NGC7793)/l)*expi(-rs_NGC7793(m200_NGC7793)/l)+expi(rs_NGC7793(m200_NGC7793)/l)-expi((r+rs_NGC7793(m200_NGC7793))/l)))

vvnfw_NGC7793= lambda r,m200_NGC7793:G0*(kpc**2)*mnfw_NGC7793(r,m200_NGC7793)/r

vv_NGC7793=lambda r,Y_NGC7793,m200_NGC7793,b,l: vgas_NGC7793(r)**2+G0*(kpc**2)*mnfw_NGC7793(r,m200_NGC7793)/r+vvnfwinmg_NGC7793(r,m200_NGC7793,b,l)+Y_NGC7793*vdisk_NGC7793(r)**2

#NGC7814

rs_NGC7814=lambda m200_NGC7814:28.8*(m200_NGC7814/(10**12)/h)**0.43

c_NGC7814=lambda m200_NGC7814:(10**0.905)*(m200_NGC7814/(10**12)/h)**(-0.101)

rho_s_NGC7814=lambda m200_NGC7814:(200/3)*(rho_c*c_NGC7814(m200_NGC7814)**3)/(np.log(1+c_NGC7814(m200_NGC7814))-c_NGC7814(m200_NGC7814)/(1+c_NGC7814(m200_NGC7814)))

mnfw_NGC7814=lambda r,m200_NGC7814:4*math.pi*(rs_NGC7814(m200_NGC7814)**3)*rho_s_NGC7814(m200_NGC7814)*(-(r/(r+rs_NGC7814(m200_NGC7814)))-np.log(rs_NGC7814(m200_NGC7814))+np.log(r+rs_NGC7814(m200_NGC7814)))

vvnfwinmg_NGC7814= lambda r,m200_NGC7814,b,l:((-2*math.pi*G0*(kpc**2)*b*rho_s_NGC7814(m200_NGC7814)*rs_NGC7814(m200_NGC7814)**3)/r)*(2*r/(r+rs_NGC7814(m200_NGC7814))+np.exp((r+rs_NGC7814(m200_NGC7814))/l)*(r/l-1)*expi(-(r+rs_NGC7814(m200_NGC7814))/l)
  +np.exp(-(r+rs_NGC7814(m200_NGC7814))/l)*(1+r/l)*(np.exp(2*rs_NGC7814(m200_NGC7814)/l)*expi(-rs_NGC7814(m200_NGC7814)/l)+expi(rs_NGC7814(m200_NGC7814)/l)-expi((r+rs_NGC7814(m200_NGC7814))/l)))

vvnfw_NGC7814= lambda r,m200_NGC7814:G0*(kpc**2)*mnfw_NGC7814(r,m200_NGC7814)/r

vv_NGC7814=lambda r,Y_NGC7814,Yb_NGC7814,m200_NGC7814,b,l: vgas_NGC7814(r)**2+G0*(kpc**2)*mnfw_NGC7814(r,m200_NGC7814)/r+vvnfwinmg_NGC7814(r,m200_NGC7814,b,l)+Y_NGC7814*vdisk_NGC7814(r)**2+Yb_NGC7814*vbulge_NGC7814(r)**2

#UGC02259

rs_UGC02259=lambda m200_UGC02259:28.8*(m200_UGC02259/(10**12)/h)**0.43

c_UGC02259=lambda m200_UGC02259:(10**0.905)*(m200_UGC02259/(10**12)/h)**(-0.101)

rho_s_UGC02259=lambda m200_UGC02259:(200/3)*(rho_c*c_UGC02259(m200_UGC02259)**3)/(np.log(1+c_UGC02259(m200_UGC02259))-c_UGC02259(m200_UGC02259)/(1+c_UGC02259(m200_UGC02259)))

mnfw_UGC02259=lambda r,m200_UGC02259:4*math.pi*(rs_UGC02259(m200_UGC02259)**3)*rho_s_UGC02259(m200_UGC02259)*(-(r/(r+rs_UGC02259(m200_UGC02259)))-np.log(rs_UGC02259(m200_UGC02259))+np.log(r+rs_UGC02259(m200_UGC02259)))

vvnfwinmg_UGC02259= lambda r,m200_UGC02259,b,l:((-2*math.pi*G0*(kpc**2)*b*rho_s_UGC02259(m200_UGC02259)*rs_UGC02259(m200_UGC02259)**3)/r)*(2*r/(r+rs_UGC02259(m200_UGC02259))+np.exp((r+rs_UGC02259(m200_UGC02259))/l)*(r/l-1)*expi(-(r+rs_UGC02259(m200_UGC02259))/l)
  +np.exp(-(r+rs_UGC02259(m200_UGC02259))/l)*(1+r/l)*(np.exp(2*rs_UGC02259(m200_UGC02259)/l)*expi(-rs_UGC02259(m200_UGC02259)/l)+expi(rs_UGC02259(m200_UGC02259)/l)-expi((r+rs_UGC02259(m200_UGC02259))/l)))

vvnfw_UGC02259= lambda r,m200_UGC02259:G0*(kpc**2)*mnfw_UGC02259(r,m200_UGC02259)/r

vv_UGC02259=lambda r,Y_UGC02259,m200_UGC02259,b,l: vgas_UGC02259(r)**2+G0*(kpc**2)*mnfw_UGC02259(r,m200_UGC02259)/r+vvnfwinmg_UGC02259(r,m200_UGC02259,b,l)+Y_UGC02259*vdisk_UGC02259(r)**2

#UGC03546

rs_UGC03546=lambda m200_UGC03546:28.8*(m200_UGC03546/(10**12)/h)**0.43

c_UGC03546=lambda m200_UGC03546:(10**0.905)*(m200_UGC03546/(10**12)/h)**(-0.101)

rho_s_UGC03546=lambda m200_UGC03546:(200/3)*(rho_c*c_UGC03546(m200_UGC03546)**3)/(np.log(1+c_UGC03546(m200_UGC03546))-c_UGC03546(m200_UGC03546)/(1+c_UGC03546(m200_UGC03546)))

mnfw_UGC03546=lambda r,m200_UGC03546:4*math.pi*(rs_UGC03546(m200_UGC03546)**3)*rho_s_UGC03546(m200_UGC03546)*(-(r/(r+rs_UGC03546(m200_UGC03546)))-np.log(rs_UGC03546(m200_UGC03546))+np.log(r+rs_UGC03546(m200_UGC03546)))

vvnfwinmg_UGC03546= lambda r,m200_UGC03546,b,l:((-2*math.pi*G0*(kpc**2)*b*rho_s_UGC03546(m200_UGC03546)*rs_UGC03546(m200_UGC03546)**3)/r)*(2*r/(r+rs_UGC03546(m200_UGC03546))+np.exp((r+rs_UGC03546(m200_UGC03546))/l)*(r/l-1)*expi(-(r+rs_UGC03546(m200_UGC03546))/l)
  +np.exp(-(r+rs_UGC03546(m200_UGC03546))/l)*(1+r/l)*(np.exp(2*rs_UGC03546(m200_UGC03546)/l)*expi(-rs_UGC03546(m200_UGC03546)/l)+expi(rs_UGC03546(m200_UGC03546)/l)-expi((r+rs_UGC03546(m200_UGC03546))/l)))

vvnfw_UGC03546= lambda r,m200_UGC03546:G0*(kpc**2)*mnfw_UGC03546(r,m200_UGC03546)/r

vv_UGC03546=lambda r,Y_UGC03546,Yb_UGC03546,m200_UGC03546,b,l: vgas_UGC03546(r)**2+G0*(kpc**2)*mnfw_UGC03546(r,m200_UGC03546)/r+vvnfwinmg_UGC03546(r,m200_UGC03546,b,l)+Y_UGC03546*vdisk_UGC03546(r)**2+Yb_UGC03546*vbulge_UGC03546(r)**2

#UGC06446

rs_UGC06446=lambda m200_UGC06446:28.8*(m200_UGC06446/(10**12)/h)**0.43

c_UGC06446=lambda m200_UGC06446:(10**0.905)*(m200_UGC06446/(10**12)/h)**(-0.101)

rho_s_UGC06446=lambda m200_UGC06446:(200/3)*(rho_c*c_UGC06446(m200_UGC06446)**3)/(np.log(1+c_UGC06446(m200_UGC06446))-c_UGC06446(m200_UGC06446)/(1+c_UGC06446(m200_UGC06446)))

mnfw_UGC06446=lambda r,m200_UGC06446:4*math.pi*(rs_UGC06446(m200_UGC06446)**3)*rho_s_UGC06446(m200_UGC06446)*(-(r/(r+rs_UGC06446(m200_UGC06446)))-np.log(rs_UGC06446(m200_UGC06446))+np.log(r+rs_UGC06446(m200_UGC06446)))

vvnfwinmg_UGC06446= lambda r,m200_UGC06446,b,l:((-2*math.pi*G0*(kpc**2)*b*rho_s_UGC06446(m200_UGC06446)*rs_UGC06446(m200_UGC06446)**3)/r)*(2*r/(r+rs_UGC06446(m200_UGC06446))+np.exp((r+rs_UGC06446(m200_UGC06446))/l)*(r/l-1)*expi(-(r+rs_UGC06446(m200_UGC06446))/l)
  +np.exp(-(r+rs_UGC06446(m200_UGC06446))/l)*(1+r/l)*(np.exp(2*rs_UGC06446(m200_UGC06446)/l)*expi(-rs_UGC06446(m200_UGC06446)/l)+expi(rs_UGC06446(m200_UGC06446)/l)-expi((r+rs_UGC06446(m200_UGC06446))/l)))

vvnfw_UGC06446= lambda r,m200_UGC06446:G0*(kpc**2)*mnfw_UGC06446(r,m200_UGC06446)/r

vv_UGC06446=lambda r,Y_UGC06446,m200_UGC06446,b,l: vgas_UGC06446(r)**2+G0*(kpc**2)*mnfw_UGC06446(r,m200_UGC06446)/r+vvnfwinmg_UGC06446(r,m200_UGC06446,b,l)+Y_UGC06446*vdisk_UGC06446(r)**2

#UGC06930

rs_UGC06930=lambda m200_UGC06930:28.8*(m200_UGC06930/(10**12)/h)**0.43

c_UGC06930=lambda m200_UGC06930:(10**0.905)*(m200_UGC06930/(10**12)/h)**(-0.101)

rho_s_UGC06930=lambda m200_UGC06930:(200/3)*(rho_c*c_UGC06930(m200_UGC06930)**3)/(np.log(1+c_UGC06930(m200_UGC06930))-c_UGC06930(m200_UGC06930)/(1+c_UGC06930(m200_UGC06930)))

mnfw_UGC06930=lambda r,m200_UGC06930:4*math.pi*(rs_UGC06930(m200_UGC06930)**3)*rho_s_UGC06930(m200_UGC06930)*(-(r/(r+rs_UGC06930(m200_UGC06930)))-np.log(rs_UGC06930(m200_UGC06930))+np.log(r+rs_UGC06930(m200_UGC06930)))

vvnfwinmg_UGC06930= lambda r,m200_UGC06930,b,l:((-2*math.pi*G0*(kpc**2)*b*rho_s_UGC06930(m200_UGC06930)*rs_UGC06930(m200_UGC06930)**3)/r)*(2*r/(r+rs_UGC06930(m200_UGC06930))+np.exp((r+rs_UGC06930(m200_UGC06930))/l)*(r/l-1)*expi(-(r+rs_UGC06930(m200_UGC06930))/l)
  +np.exp(-(r+rs_UGC06930(m200_UGC06930))/l)*(1+r/l)*(np.exp(2*rs_UGC06930(m200_UGC06930)/l)*expi(-rs_UGC06930(m200_UGC06930)/l)+expi(rs_UGC06930(m200_UGC06930)/l)-expi((r+rs_UGC06930(m200_UGC06930))/l)))

vvnfw_UGC06930= lambda r,m200_UGC06930:G0*(kpc**2)*mnfw_UGC06930(r,m200_UGC06930)/r

vv_UGC06930=lambda r,Y_UGC06930,m200_UGC06930,b,l: vgas_UGC06930(r)**2+G0*(kpc**2)*mnfw_UGC06930(r,m200_UGC06930)/r+vvnfwinmg_UGC06930(r,m200_UGC06930,b,l)+Y_UGC06930*vdisk_UGC06930(r)**2

#UGC06983

rs_UGC06983=lambda m200_UGC06983:28.8*(m200_UGC06983/(10**12)/h)**0.43

c_UGC06983=lambda m200_UGC06983:(10**0.905)*(m200_UGC06983/(10**12)/h)**(-0.101)

rho_s_UGC06983=lambda m200_UGC06983:(200/3)*(rho_c*c_UGC06983(m200_UGC06983)**3)/(np.log(1+c_UGC06983(m200_UGC06983))-c_UGC06983(m200_UGC06983)/(1+c_UGC06983(m200_UGC06983)))

mnfw_UGC06983=lambda r,m200_UGC06983:4*math.pi*(rs_UGC06983(m200_UGC06983)**3)*rho_s_UGC06983(m200_UGC06983)*(-(r/(r+rs_UGC06983(m200_UGC06983)))-np.log(rs_UGC06983(m200_UGC06983))+np.log(r+rs_UGC06983(m200_UGC06983)))

vvnfwinmg_UGC06983= lambda r,m200_UGC06983,b,l:((-2*math.pi*G0*(kpc**2)*b*rho_s_UGC06983(m200_UGC06983)*rs_UGC06983(m200_UGC06983)**3)/r)*(2*r/(r+rs_UGC06983(m200_UGC06983))+np.exp((r+rs_UGC06983(m200_UGC06983))/l)*(r/l-1)*expi(-(r+rs_UGC06983(m200_UGC06983))/l)
  +np.exp(-(r+rs_UGC06983(m200_UGC06983))/l)*(1+r/l)*(np.exp(2*rs_UGC06983(m200_UGC06983)/l)*expi(-rs_UGC06983(m200_UGC06983)/l)+expi(rs_UGC06983(m200_UGC06983)/l)-expi((r+rs_UGC06983(m200_UGC06983))/l)))

vvnfw_UGC06983= lambda r,m200_UGC06983:G0*(kpc**2)*mnfw_UGC06983(r,m200_UGC06983)/r

vv_UGC06983=lambda r,Y_UGC06983,m200_UGC06983,b,l: vgas_UGC06983(r)**2+G0*(kpc**2)*mnfw_UGC06983(r,m200_UGC06983)/r+vvnfwinmg_UGC06983(r,m200_UGC06983,b,l)+Y_UGC06983*vdisk_UGC06983(r)**2

#UGC07261

rs_UGC07261=lambda m200_UGC07261:28.8*(m200_UGC07261/(10**12)/h)**0.43

c_UGC07261=lambda m200_UGC07261:(10**0.905)*(m200_UGC07261/(10**12)/h)**(-0.101)

rho_s_UGC07261=lambda m200_UGC07261:(200/3)*(rho_c*c_UGC07261(m200_UGC07261)**3)/(np.log(1+c_UGC07261(m200_UGC07261))-c_UGC07261(m200_UGC07261)/(1+c_UGC07261(m200_UGC07261)))

mnfw_UGC07261=lambda r,m200_UGC07261:4*math.pi*(rs_UGC07261(m200_UGC07261)**3)*rho_s_UGC07261(m200_UGC07261)*(-(r/(r+rs_UGC07261(m200_UGC07261)))-np.log(rs_UGC07261(m200_UGC07261))+np.log(r+rs_UGC07261(m200_UGC07261)))

vvnfwinmg_UGC07261= lambda r,m200_UGC07261,b,l:((-2*math.pi*G0*(kpc**2)*b*rho_s_UGC07261(m200_UGC07261)*rs_UGC07261(m200_UGC07261)**3)/r)*(2*r/(r+rs_UGC07261(m200_UGC07261))+np.exp((r+rs_UGC07261(m200_UGC07261))/l)*(r/l-1)*expi(-(r+rs_UGC07261(m200_UGC07261))/l)
  +np.exp(-(r+rs_UGC07261(m200_UGC07261))/l)*(1+r/l)*(np.exp(2*rs_UGC07261(m200_UGC07261)/l)*expi(-rs_UGC07261(m200_UGC07261)/l)+expi(rs_UGC07261(m200_UGC07261)/l)-expi((r+rs_UGC07261(m200_UGC07261))/l)))

vvnfw_UGC07261= lambda r,m200_UGC07261:G0*(kpc**2)*mnfw_UGC07261(r,m200_UGC07261)/r

vv_UGC07261=lambda r,Y_UGC07261,m200_UGC07261,b,l: vgas_UGC07261(r)**2+G0*(kpc**2)*mnfw_UGC07261(r,m200_UGC07261)/r+vvnfwinmg_UGC07261(r,m200_UGC07261,b,l)+Y_UGC07261*vdisk_UGC07261(r)**2

#UGC07690

rs_UGC07690=lambda m200_UGC07690:28.8*(m200_UGC07690/(10**12)/h)**0.43

c_UGC07690=lambda m200_UGC07690:(10**0.905)*(m200_UGC07690/(10**12)/h)**(-0.101)

rho_s_UGC07690=lambda m200_UGC07690:(200/3)*(rho_c*c_UGC07690(m200_UGC07690)**3)/(np.log(1+c_UGC07690(m200_UGC07690))-c_UGC07690(m200_UGC07690)/(1+c_UGC07690(m200_UGC07690)))

mnfw_UGC07690=lambda r,m200_UGC07690:4*math.pi*(rs_UGC07690(m200_UGC07690)**3)*rho_s_UGC07690(m200_UGC07690)*(-(r/(r+rs_UGC07690(m200_UGC07690)))-np.log(rs_UGC07690(m200_UGC07690))+np.log(r+rs_UGC07690(m200_UGC07690)))

vvnfwinmg_UGC07690= lambda r,m200_UGC07690,b,l:((-2*math.pi*G0*(kpc**2)*b*rho_s_UGC07690(m200_UGC07690)*rs_UGC07690(m200_UGC07690)**3)/r)*(2*r/(r+rs_UGC07690(m200_UGC07690))+np.exp((r+rs_UGC07690(m200_UGC07690))/l)*(r/l-1)*expi(-(r+rs_UGC07690(m200_UGC07690))/l)
  +np.exp(-(r+rs_UGC07690(m200_UGC07690))/l)*(1+r/l)*(np.exp(2*rs_UGC07690(m200_UGC07690)/l)*expi(-rs_UGC07690(m200_UGC07690)/l)+expi(rs_UGC07690(m200_UGC07690)/l)-expi((r+rs_UGC07690(m200_UGC07690))/l)))

vvnfw_UGC07690= lambda r,m200_UGC07690:G0*(kpc**2)*mnfw_UGC07690(r,m200_UGC07690)/r

vv_UGC07690=lambda r,Y_UGC07690,m200_UGC07690,b,l: vgas_UGC07690(r)**2+G0*(kpc**2)*mnfw_UGC07690(r,m200_UGC07690)/r+vvnfwinmg_UGC07690(r,m200_UGC07690,b,l)+Y_UGC07690*vdisk_UGC07690(r)**2

########################################################################################################################################################

#Running MCMC
#############################
print('emcee version=',emcee.__version__)

x=radius_NGC7331,radius_NGC7793,radius_NGC7814,radius_UGC02259,radius_UGC03546,radius_UGC06446,radius_UGC06930,radius_UGC06983,radius_UGC07261,radius_UGC07690

y=velocity_NGC7331,velocity_NGC7793,velocity_NGC7814,velocity_UGC02259,velocity_UGC03546,velocity_UGC06446,velocity_UGC06930,velocity_UGC06983,velocity_UGC07261,velocity_UGC07690

yerr=error_velocity_NGC7331,error_velocity_NGC7793,error_velocity_NGC7814,error_velocity_UGC02259,error_velocity_UGC03546,error_velocity_UGC06446,error_velocity_UGC06930,error_velocity_UGC06983,error_velocity_UGC07261,error_velocity_UGC07690

def log_likelihood(theta, x, y, yerr):
    
 Y_NGC7331,Y_NGC7793,Y_NGC7814,Y_UGC02259,Y_UGC03546,Y_UGC06446,Y_UGC06930,Y_UGC06983,Y_UGC07261,Y_UGC07690,m200_NGC7331,m200_NGC7793,m200_NGC7814,m200_UGC02259,m200_UGC03546,m200_UGC06446,m200_UGC06930,m200_UGC06983,m200_UGC07261,m200_UGC07690,Yb_NGC7331,Yb_NGC7814,Yb_UGC03546,b,l = theta    

 model_NGC7331=vgas_NGC7331(radius_NGC7331)**2+G0*(kpc**2)*mnfw_NGC7331(radius_NGC7331,m200_NGC7331)/radius_NGC7331+vvnfwinmg_NGC7331(radius_NGC7331,m200_NGC7331,b,l)+Y_NGC7331*vdisk_NGC7331(radius_NGC7331)**2+Yb_NGC7331*vbulge_NGC7331(radius_NGC7331)**2
 inv_sigma2_NGC7331=1.0/(4*velocity_NGC7331**2*error_velocity_NGC7331**2)

 model_NGC7793=vgas_NGC7793(radius_NGC7793)**2+G0*(kpc**2)*mnfw_NGC7793(radius_NGC7793,m200_NGC7793)/radius_NGC7793+vvnfwinmg_NGC7793(radius_NGC7793,m200_NGC7793,b,l)+Y_NGC7793*vdisk_NGC7793(radius_NGC7793)**2
 inv_sigma2_NGC7793=1.0/(4*velocity_NGC7793**2*error_velocity_NGC7793**2)

 model_NGC7814=vgas_NGC7814(radius_NGC7814)**2+G0*(kpc**2)*mnfw_NGC7814(radius_NGC7814,m200_NGC7814)/radius_NGC7814+vvnfwinmg_NGC7814(radius_NGC7814,m200_NGC7814,b,l)+Y_NGC7814*vdisk_NGC7814(radius_NGC7814)**2+Yb_NGC7814*vbulge_NGC7814(radius_NGC7814)**2
 inv_sigma2_NGC7814=1.0/(4*velocity_NGC7814**2*error_velocity_NGC7814**2)

 model_UGC02259=vgas_UGC02259(radius_UGC02259)**2+G0*(kpc**2)*mnfw_UGC02259(radius_UGC02259,m200_UGC02259)/radius_UGC02259+vvnfwinmg_UGC02259(radius_UGC02259,m200_UGC02259,b,l)+Y_UGC02259*vdisk_UGC02259(radius_UGC02259)**2
 inv_sigma2_UGC02259=1.0/(4*velocity_UGC02259**2*error_velocity_UGC02259**2)

 model_UGC03546=vgas_UGC03546(radius_UGC03546)**2+G0*(kpc**2)*mnfw_UGC03546(radius_UGC03546,m200_UGC03546)/radius_UGC03546+vvnfwinmg_UGC03546(radius_UGC03546,m200_UGC03546,b,l)+Y_UGC03546*vdisk_UGC03546(radius_UGC03546)**2+Yb_UGC03546*vbulge_UGC03546(radius_UGC03546)**2
 inv_sigma2_UGC03546=1.0/(4*velocity_UGC03546**2*error_velocity_UGC03546**2)

 model_UGC06446=vgas_UGC06446(radius_UGC06446)**2+G0*(kpc**2)*mnfw_UGC06446(radius_UGC06446,m200_UGC06446)/radius_UGC06446+vvnfwinmg_UGC06446(radius_UGC06446,m200_UGC06446,b,l)+Y_UGC06446*vdisk_UGC06446(radius_UGC06446)**2
 inv_sigma2_UGC06446=1.0/(4*velocity_UGC06446**2*error_velocity_UGC06446**2)

 model_UGC06930=vgas_UGC06930(radius_UGC06930)**2+G0*(kpc**2)*mnfw_UGC06930(radius_UGC06930,m200_UGC06930)/radius_UGC06930+vvnfwinmg_UGC06930(radius_UGC06930,m200_UGC06930,b,l)+Y_UGC06930*vdisk_UGC06930(radius_UGC06930)**2
 inv_sigma2_UGC06930=1.0/(4*velocity_UGC06930**2*error_velocity_UGC06930**2)

 model_UGC06983=vgas_UGC06983(radius_UGC06983)**2+G0*(kpc**2)*mnfw_UGC06983(radius_UGC06983,m200_UGC06983)/radius_UGC06983+vvnfwinmg_UGC06983(radius_UGC06983,m200_UGC06983,b,l)+Y_UGC06983*vdisk_UGC06983(radius_UGC06983)**2
 inv_sigma2_UGC06983=1.0/(4*velocity_UGC06983**2*error_velocity_UGC06983**2)

 model_UGC07261=vgas_UGC07261(radius_UGC07261)**2+G0*(kpc**2)*mnfw_UGC07261(radius_UGC07261,m200_UGC07261)/radius_UGC07261+vvnfwinmg_UGC07261(radius_UGC07261,m200_UGC07261,b,l)+Y_UGC07261*vdisk_UGC07261(radius_UGC07261)**2
 inv_sigma2_UGC07261=1.0/(4*velocity_UGC07261**2*error_velocity_UGC07261**2)

 model_UGC07690=vgas_UGC07690(radius_UGC07690)**2+G0*(kpc**2)*mnfw_UGC07690(radius_UGC07690,m200_UGC07690)/radius_UGC07690+vvnfwinmg_UGC07690(radius_UGC07690,m200_UGC07690,b,l)+Y_UGC07690*vdisk_UGC07690(radius_UGC07690)**2
 inv_sigma2_UGC07690=1.0/(4*velocity_UGC07690**2*error_velocity_UGC07690**2)
     
 loglikelihood = -0.5*(np.sum((velocity_NGC7331**2-model_NGC7331)**2*inv_sigma2_NGC7331))-0.5*(np.sum((velocity_NGC7793**2-model_NGC7793)**2*inv_sigma2_NGC7793))-0.5*(np.sum((velocity_NGC7814**2-model_NGC7814)**2*inv_sigma2_NGC7814))-0.5*(np.sum((velocity_UGC02259**2-model_UGC02259)**2*inv_sigma2_UGC02259))-0.5*(np.sum((velocity_UGC03546**2-model_UGC03546)**2*inv_sigma2_UGC03546))-0.5*(np.sum((velocity_UGC06446**2-model_UGC06446)**2*inv_sigma2_UGC06446))-0.5*(np.sum((velocity_UGC06930**2-model_UGC06930)**2*inv_sigma2_UGC06930))-0.5*(np.sum((velocity_UGC06983**2-model_UGC06983)**2*inv_sigma2_UGC06983))-0.5*(np.sum((velocity_UGC07261**2-model_UGC07261)**2*inv_sigma2_UGC07261))-0.5*(np.sum((velocity_UGC07690**2-model_UGC07690)**2*inv_sigma2_UGC07690))      
 
 return loglikelihood

#minimum values for lambda:

l_min=(radius_NGC7331[0]+radius_NGC7793[0]+radius_NGC7814[0]+radius_UGC02259[0]+radius_UGC03546[0]+radius_UGC06446[0]+radius_UGC06930[0]+radius_UGC06983[0]+radius_UGC07261[0]+radius_UGC07690[0])/10

def log_prior(theta):
 
 Y_NGC7331,Y_NGC7793,Y_NGC7814,Y_UGC02259,Y_UGC03546,Y_UGC06446,Y_UGC06930,Y_UGC06983,Y_UGC07261,Y_UGC07690,m200_NGC7331,m200_NGC7793,m200_NGC7814,m200_UGC02259,m200_UGC03546,m200_UGC06446,m200_UGC06930,m200_UGC06983,m200_UGC07261,m200_UGC07690,Yb_NGC7331,Yb_NGC7814,Yb_UGC03546,b,l = theta
   
 if ((0.3<Y_NGC7331<0.8) and
 (0.3<Y_NGC7793<0.8) and
 (0.3<Y_NGC7814<0.8) and
 (0.3<Y_UGC02259<0.8) and
 (0.3<Y_UGC03546<0.8) and
 (0.3<Y_UGC06446<0.8) and
 (0.3<Y_UGC06930<0.8) and
 (0.3<Y_UGC06983<0.8) and
 (0.3<Y_UGC07261<0.8) and
 (0.3<Y_UGC07690<0.8) and
 (10**9<m200_NGC7331<10**14) and
 (10**9<m200_NGC7793<10**14) and
 (10**9<m200_NGC7814<10**14) and
 (10**9<m200_UGC02259<10**14) and
 (10**9<m200_UGC03546<10**14) and
 (10**9<m200_UGC06446<10**14) and
 (10**9<m200_UGC06930<10**14) and
 (10**9<m200_UGC06983<10**14) and
 (10**9<m200_UGC07261<10**14) and
 (10**9<m200_UGC07690<10**14) and
 (0.3<Yb_NGC7331<0.8) and
 (0.3<Yb_NGC7814<0.8) and
 (0.3<Yb_UGC03546<0.8) and
 (-2 < b < 2) and  (l_min < l < 100)):       
     
    return 0.0
 return -np.inf

def log_probability(theta, x, y, yerr):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, x, y, yerr)

my_guess=[0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5*10**12,0.5*10**12,0.5*10**12,0.5*10**12,0.5*10**12,0.5*10**12,0.5*10**12,0.5*10**12,0.5*10**12,0.5*10**12,0.5,0.5,0.5,0,l_min+1]

pos = my_guess + 1e-4*np.random.randn(100, len(my_guess))
nwalkers, ndim = pos.shape
burn_in=2*10**4 # based on autocorrelation time analysis burn_in=tau
nsteps =70*burn_in # recommended by the emcee developers N > 50*tau
thin=burn_in/2

from multiprocessing import Pool
from contextlib import closing
from multiprocessing import cpu_count
ncpu = cpu_count()
print("{0} CPUs".format(ncpu))

# Here, if I want to change the stretch move 
my_moves= emcee.moves.StretchMove()

with closing(Pool()) as pool:
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, moves= my_moves,pool=pool,args=(x, y, yerr))
    start = time.time()
    sampler.run_mcmc(pos, nsteps, progress=True)
    end = time.time()
    multi_time = end - start
    pool.terminate()
    print("Multiprocessing took {0:.1f} seconds".format(multi_time))
    
    
chain=sampler.get_chain()[:, :, 0].T

accep_frac=sampler.acceptance_fraction
print('Acceptance fraction=',accep_frac.sum()/nwalkers)

tau=emcee.autocorr.integrated_time(chain)
print('Integrated autocorrelation time=',tau)

flat_chain = sampler.get_chain(discard= burn_in, thin=thin, flat=True)

import os
flat_chain_export = pd.DataFrame(flat_chain)
flat_chain_export.to_csv("/remote/pi213b/almeida/Python/new_emcee_testing_10_galaxy/flatchain_testing_10_galaxy_set4.csv")

#Maximum likelihod estimators
############################################
print('Best fit values')
print('-------------------------------------')
labels=["Y_NGC7331", "Y_NGC7793", "Y_NGC7814", "Y_UGC02259", "Y_UGC03546", "Y_UGC06446", "Y_UGC06930", "Y_UGC06983", "Y_UGC07261", "Y_UGC07690", "m200_NGC7331", "m200_NGC7793", "m200_NGC7814", "m200_UGC02259", 
        "m200_UGC03546", "m200_UGC06446", "m200_UGC06930", "m200_UGC06983", "m200_UGC07261", "m200_UGC07690", "Yb_NGC7331", "Yb_NGC7814", "Yb_UGC03546",'beta','lambda']

mle_soln=[]

for i in range(ndim):
    mcmc = np.percentile(flat_chain[:, i], [16, 50, 84])
    q = np.diff(mcmc)
    mle_soln.append(mcmc[1])
    print(labels[i],'=',mcmc[1],'+',q[0],'-',q[1])

#Best fit values
##############################################
Y_NGC7331_sol=mle_soln[0]
Y_NGC7793_sol=mle_soln[1]
Y_NGC7814_sol=mle_soln[2]
Y_UGC02259_sol=mle_soln[3]
Y_UGC03546_sol=mle_soln[4]
Y_UGC06446_sol=mle_soln[5]
Y_UGC06930_sol=mle_soln[6]
Y_UGC06983_sol=mle_soln[7]
Y_UGC07261_sol=mle_soln[8]
Y_UGC07690_sol=mle_soln[9]
m200_NGC7331_sol=mle_soln[10]
m200_NGC7793_sol=mle_soln[11]
m200_NGC7814_sol=mle_soln[12]
m200_UGC02259_sol=mle_soln[13]
m200_UGC03546_sol=mle_soln[14]
m200_UGC06446_sol=mle_soln[15]
m200_UGC06930_sol=mle_soln[16]
m200_UGC06983_sol=mle_soln[17]
m200_UGC07261_sol=mle_soln[18]
m200_UGC07690_sol=mle_soln[19]
Yb_NGC7331_sol=mle_soln[20]
Yb_NGC7814_sol=mle_soln[21]
Yb_UGC03546_sol=mle_soln[22]
bsol=mle_soln[23]
lsol=mle_soln[24]

print('---------------------------------------------')
print('Goodness of fit')
print('---------------------------------------------')

sigma2_NGC7331=1.0/(error_velocity_NGC7331**2)
sigma2_NGC7793=1.0/(error_velocity_NGC7793**2)
sigma2_NGC7814=1.0/(error_velocity_NGC7814**2)
sigma2_UGC02259=1.0/(error_velocity_UGC02259**2)
sigma2_UGC03546=1.0/(error_velocity_UGC03546**2)
sigma2_UGC06446=1.0/(error_velocity_UGC06446**2)
sigma2_UGC06930=1.0/(error_velocity_UGC06930**2)
sigma2_UGC06983=1.0/(error_velocity_UGC06983**2)
sigma2_UGC07261=1.0/(error_velocity_UGC07261**2)
sigma2_UGC07690=1.0/(error_velocity_UGC07690**2)

dof_NGC7331=len(velocity_NGC7331)-5
dof_NGC7793=len(velocity_NGC7793)-4
dof_NGC7814=len(velocity_NGC7814)-5
dof_UGC02259=len(velocity_UGC02259)-4
dof_UGC03546=len(velocity_UGC03546)-5
dof_UGC06446=len(velocity_UGC06446)-4
dof_UGC06930=len(velocity_UGC06930)-4
dof_UGC06983=len(velocity_UGC06983)-4
dof_UGC07261=len(velocity_UGC07261)-4
dof_UGC07690=len(velocity_UGC07690)-4

SchiNGC7331=np.sum((velocity_NGC7331-np.sqrt(vv_NGC7331(radius_NGC7331,Y_NGC7331_sol,Yb_NGC7331_sol,m200_NGC7331_sol,bsol,lsol)))**2*sigma2_NGC7331)
SchiNGC7793=np.sum((velocity_NGC7793-np.sqrt(vv_NGC7793(radius_NGC7793,Y_NGC7793_sol,m200_NGC7793_sol,bsol,lsol)))**2*sigma2_NGC7793)
SchiNGC7814=np.sum((velocity_NGC7814-np.sqrt(vv_NGC7814(radius_NGC7814,Y_NGC7814_sol,Yb_NGC7814_sol,m200_NGC7814_sol,bsol,lsol)))**2*sigma2_NGC7814)
SchiUGC02259=np.sum((velocity_UGC02259-np.sqrt(vv_UGC02259(radius_UGC02259,Y_UGC02259_sol,m200_UGC02259_sol,bsol,lsol)))**2*sigma2_UGC02259)
SchiUGC03546=np.sum((velocity_UGC03546-np.sqrt(vv_UGC03546(radius_UGC03546,Y_UGC03546_sol,Yb_UGC03546_sol,m200_UGC03546_sol,bsol,lsol)))**2*sigma2_UGC03546)
SchiUGC06446=np.sum((velocity_UGC06446-np.sqrt(vv_UGC06446(radius_UGC06446,Y_UGC06446_sol,m200_UGC06446_sol,bsol,lsol)))**2*sigma2_UGC06446)
SchiUGC06930=np.sum((velocity_UGC06930-np.sqrt(vv_UGC06930(radius_UGC06930,Y_UGC06930_sol,m200_UGC06930_sol,bsol,lsol)))**2*sigma2_UGC06930)
SchiUGC06983=np.sum((velocity_UGC06983-np.sqrt(vv_UGC06983(radius_UGC06983,Y_UGC06983_sol,m200_UGC06983_sol,bsol,lsol)))**2*sigma2_UGC06983)
SchiUGC07261=np.sum((velocity_UGC07261-np.sqrt(vv_UGC07261(radius_UGC07261,Y_UGC07261_sol,m200_UGC07261_sol,bsol,lsol)))**2*sigma2_UGC07261)
SchiUGC07690=np.sum((velocity_UGC07690-np.sqrt(vv_UGC07690(radius_UGC07690,Y_UGC07690_sol,m200_UGC07690_sol,bsol,lsol)))**2*sigma2_UGC07690)

Schi=SchiNGC7331+SchiNGC7793+SchiNGC7814+SchiUGC02259+SchiUGC03546+SchiUGC06446+SchiUGC06930+SchiUGC06983+SchiUGC07261+SchiUGC07690

data_total=len(velocity_NGC7331)+len(velocity_NGC7793)+len(velocity_NGC7814)+len(velocity_UGC02259)+len(velocity_UGC03546)+len(velocity_UGC06446)+len(velocity_UGC06930)+len(velocity_UGC06983)+len(velocity_UGC07261)+len(velocity_UGC07690)

dof=data_total-len(my_guess)

print('Schired_NGC7331:',(1.0/dof_NGC7331)*(np.sum((velocity_NGC7331-np.sqrt(vv_NGC7331(radius_NGC7331,Y_NGC7331_sol,Yb_NGC7331_sol,m200_NGC7331_sol,bsol,lsol)))**2*sigma2_NGC7331)))
print('Schired_NGC7793:',(1.0/dof_NGC7793)*(np.sum((velocity_NGC7793-np.sqrt(vv_NGC7793(radius_NGC7793,Y_NGC7793_sol,m200_NGC7793_sol,bsol,lsol)))**2*sigma2_NGC7793)))
print('Schired_NGC7814:',(1.0/dof_NGC7814)*(np.sum((velocity_NGC7814-np.sqrt(vv_NGC7814(radius_NGC7814,Y_NGC7814_sol,Yb_NGC7814_sol,m200_NGC7814_sol,bsol,lsol)))**2*sigma2_NGC7814)))
print('Schired_UGC02259:',(1.0/dof_UGC02259)*(np.sum((velocity_UGC02259-np.sqrt(vv_UGC02259(radius_UGC02259,Y_UGC02259_sol,m200_UGC02259_sol,bsol,lsol)))**2*sigma2_UGC02259)))
print('Schired_UGC03546:',(1.0/dof_UGC03546)*(np.sum((velocity_UGC03546-np.sqrt(vv_UGC03546(radius_UGC03546,Y_UGC03546_sol,Yb_UGC03546_sol,m200_UGC03546_sol,bsol,lsol)))**2*sigma2_UGC03546)))
print('Schired_UGC06446:',(1.0/dof_UGC06446)*(np.sum((velocity_UGC06446-np.sqrt(vv_UGC06446(radius_UGC06446,Y_UGC06446_sol,m200_UGC06446_sol,bsol,lsol)))**2*sigma2_UGC06446)))
print('Schired_UGC06930:',(1.0/dof_UGC06930)*(np.sum((velocity_UGC06930-np.sqrt(vv_UGC06930(radius_UGC06930,Y_UGC06930_sol,m200_UGC06930_sol,bsol,lsol)))**2*sigma2_UGC06930)))
print('Schired_UGC06983:',(1.0/dof_UGC06983)*(np.sum((velocity_UGC06983-np.sqrt(vv_UGC06983(radius_UGC06983,Y_UGC06983_sol,m200_UGC06983_sol,bsol,lsol)))**2*sigma2_UGC06983)))
print('Schired_UGC07261:',(1.0/dof_UGC07261)*(np.sum((velocity_UGC07261-np.sqrt(vv_UGC07261(radius_UGC07261,Y_UGC07261_sol,m200_UGC07261_sol,bsol,lsol)))**2*sigma2_UGC07261)))
print('Schired_UGC07690:',(1.0/dof_UGC07690)*(np.sum((velocity_UGC07690-np.sqrt(vv_UGC07690(radius_UGC07690,Y_UGC07690_sol,m200_UGC07690_sol,bsol,lsol)))**2*sigma2_UGC07690)))
print('Schired:',Schi/dof)

#NGC7331

plt.xlim([0,radius_NGC7331[-1]+1])

plt.errorbar(radius_NGC7331,velocity_NGC7331,error_velocity_NGC7331,fmt='ro',lw=2)

plt.plot(radius_NGC7331,vgas_NGC7331(radius_NGC7331),'y',lw=2)

plt.plot(radius_NGC7331,np.sqrt(Y_NGC7331_sol*vdisk_NGC7331(radius_NGC7331)**2),'g',lw=2)

plt.plot(radius_NGC7331,np.sqrt(Yb_NGC7331_sol*vbulge_NGC7331(radius_NGC7331)**2),'r',lw=2)

plt.plot(radius_NGC7331,np.sqrt(vvnfw_NGC7331(radius_NGC7331,m200_NGC7331_sol)+vvnfwinmg_NGC7331(radius_NGC7331,m200_NGC7331_sol,bsol,lsol)),'b',lw=2)

plt.plot(radius_NGC7331,np.sqrt(vv_NGC7331(radius_NGC7331,Y_NGC7331_sol,Yb_NGC7331_sol,m200_NGC7331_sol,bsol,lsol)),'k',lw=2)

plt.xlabel('Radius (kpc)',fontsize=18)

plt.ylabel('Velocity (km/s)',fontsize=18)

plt.title('NGC7331')

import os

filenameNGC7331='result_NGC7331_new_emcee_testing_10_galaxy_set4.pdf'

pathNGC7331='/remote/pi213b/almeida/Python/new_emcee_testing_10_galaxy'

fullpathNGC7331=os.path.join(pathNGC7331,filenameNGC7331)

plt.savefig(fullpathNGC7331)

plt.close()

#NGC7793

plt.xlim([0,radius_NGC7793[-1]+1])

plt.errorbar(radius_NGC7793,velocity_NGC7793,error_velocity_NGC7793,fmt='ro',lw=2)

plt.plot(radius_NGC7793,vgas_NGC7793(radius_NGC7793),'y',lw=2)

plt.plot(radius_NGC7793,np.sqrt(Y_NGC7793_sol*vdisk_NGC7793(radius_NGC7793)**2),'g',lw=2)

#no bulge

plt.plot(radius_NGC7793,np.sqrt(vvnfw_NGC7793(radius_NGC7793,m200_NGC7793_sol)+vvnfwinmg_NGC7793(radius_NGC7793,m200_NGC7793_sol,bsol,lsol)),'b',lw=2)

plt.plot(radius_NGC7793,np.sqrt(vv_NGC7793(radius_NGC7793,Y_NGC7793_sol,m200_NGC7793_sol,bsol,lsol)),'k',lw=2)

plt.xlabel('Radius (kpc)',fontsize=18)

plt.ylabel('Velocity (km/s)',fontsize=18)

plt.title('NGC7793')

import os

filenameNGC7793='result_NGC7793_new_emcee_testing_10_galaxy_set4.pdf'

pathNGC7793='/remote/pi213b/almeida/Python/new_emcee_testing_10_galaxy'

fullpathNGC7793=os.path.join(pathNGC7793,filenameNGC7793)

plt.savefig(fullpathNGC7793)

plt.close()

#NGC7814

plt.xlim([0,radius_NGC7814[-1]+1])

plt.errorbar(radius_NGC7814,velocity_NGC7814,error_velocity_NGC7814,fmt='ro',lw=2)

plt.plot(radius_NGC7814,vgas_NGC7814(radius_NGC7814),'y',lw=2)

plt.plot(radius_NGC7814,np.sqrt(Y_NGC7814_sol*vdisk_NGC7814(radius_NGC7814)**2),'g',lw=2)

plt.plot(radius_NGC7814,np.sqrt(Yb_NGC7814_sol*vbulge_NGC7814(radius_NGC7814)**2),'r',lw=2)

plt.plot(radius_NGC7814,np.sqrt(vvnfw_NGC7814(radius_NGC7814,m200_NGC7814_sol)+vvnfwinmg_NGC7814(radius_NGC7814,m200_NGC7814_sol,bsol,lsol)),'b',lw=2)

plt.plot(radius_NGC7814,np.sqrt(vv_NGC7814(radius_NGC7814,Y_NGC7814_sol,Yb_NGC7814_sol,m200_NGC7814_sol,bsol,lsol)),'k',lw=2)

plt.xlabel('Radius (kpc)',fontsize=18)

plt.ylabel('Velocity (km/s)',fontsize=18)

plt.title('NGC7814')

import os

filenameNGC7814='result_NGC7814_new_emcee_testing_10_galaxy_set4.pdf'

pathNGC7814='/remote/pi213b/almeida/Python/new_emcee_testing_10_galaxy'

fullpathNGC7814=os.path.join(pathNGC7814,filenameNGC7814)

plt.savefig(fullpathNGC7814)

plt.close()

#UGC02259

plt.xlim([0,radius_UGC02259[-1]+1])

plt.errorbar(radius_UGC02259,velocity_UGC02259,error_velocity_UGC02259,fmt='ro',lw=2)

plt.plot(radius_UGC02259,vgas_UGC02259(radius_UGC02259),'y',lw=2)

plt.plot(radius_UGC02259,np.sqrt(Y_UGC02259_sol*vdisk_UGC02259(radius_UGC02259)**2),'g',lw=2)

#no bulge

plt.plot(radius_UGC02259,np.sqrt(vvnfw_UGC02259(radius_UGC02259,m200_UGC02259_sol)+vvnfwinmg_UGC02259(radius_UGC02259,m200_UGC02259_sol,bsol,lsol)),'b',lw=2)

plt.plot(radius_UGC02259,np.sqrt(vv_UGC02259(radius_UGC02259,Y_UGC02259_sol,m200_UGC02259_sol,bsol,lsol)),'k',lw=2)

plt.xlabel('Radius (kpc)',fontsize=18)

plt.ylabel('Velocity (km/s)',fontsize=18)

plt.title('UGC02259')

import os

filenameUGC02259='result_UGC02259_new_emcee_testing_10_galaxy_set4.pdf'

pathUGC02259='/remote/pi213b/almeida/Python/new_emcee_testing_10_galaxy'

fullpathUGC02259=os.path.join(pathUGC02259,filenameUGC02259)

plt.savefig(fullpathUGC02259)

plt.close()

#UGC03546

plt.xlim([0,radius_UGC03546[-1]+1])

plt.errorbar(radius_UGC03546,velocity_UGC03546,error_velocity_UGC03546,fmt='ro',lw=2)

plt.plot(radius_UGC03546,vgas_UGC03546(radius_UGC03546),'y',lw=2)

plt.plot(radius_UGC03546,np.sqrt(Y_UGC03546_sol*vdisk_UGC03546(radius_UGC03546)**2),'g',lw=2)

plt.plot(radius_UGC03546,np.sqrt(Yb_UGC03546_sol*vbulge_UGC03546(radius_UGC03546)**2),'r',lw=2)

plt.plot(radius_UGC03546,np.sqrt(vvnfw_UGC03546(radius_UGC03546,m200_UGC03546_sol)+vvnfwinmg_UGC03546(radius_UGC03546,m200_UGC03546_sol,bsol,lsol)),'b',lw=2)

plt.plot(radius_UGC03546,np.sqrt(vv_UGC03546(radius_UGC03546,Y_UGC03546_sol,Yb_UGC03546_sol,m200_UGC03546_sol,bsol,lsol)),'k',lw=2)

plt.xlabel('Radius (kpc)',fontsize=18)

plt.ylabel('Velocity (km/s)',fontsize=18)

plt.title('UGC03546')

import os

filenameUGC03546='result_UGC03546_new_emcee_testing_10_galaxy_set4.pdf'

pathUGC03546='/remote/pi213b/almeida/Python/new_emcee_testing_10_galaxy'

fullpathUGC03546=os.path.join(pathUGC03546,filenameUGC03546)

plt.savefig(fullpathUGC03546)

plt.close()

#UGC06446

plt.xlim([0,radius_UGC06446[-1]+1])

plt.errorbar(radius_UGC06446,velocity_UGC06446,error_velocity_UGC06446,fmt='ro',lw=2)

plt.plot(radius_UGC06446,vgas_UGC06446(radius_UGC06446),'y',lw=2)

plt.plot(radius_UGC06446,np.sqrt(Y_UGC06446_sol*vdisk_UGC06446(radius_UGC06446)**2),'g',lw=2)

#no bulge

plt.plot(radius_UGC06446,np.sqrt(vvnfw_UGC06446(radius_UGC06446,m200_UGC06446_sol)+vvnfwinmg_UGC06446(radius_UGC06446,m200_UGC06446_sol,bsol,lsol)),'b',lw=2)

plt.plot(radius_UGC06446,np.sqrt(vv_UGC06446(radius_UGC06446,Y_UGC06446_sol,m200_UGC06446_sol,bsol,lsol)),'k',lw=2)

plt.xlabel('Radius (kpc)',fontsize=18)

plt.ylabel('Velocity (km/s)',fontsize=18)

plt.title('UGC06446')

import os

filenameUGC06446='result_UGC06446_new_emcee_testing_10_galaxy_set4.pdf'

pathUGC06446='/remote/pi213b/almeida/Python/new_emcee_testing_10_galaxy'

fullpathUGC06446=os.path.join(pathUGC06446,filenameUGC06446)

plt.savefig(fullpathUGC06446)

plt.close()

#UGC06930

plt.xlim([0,radius_UGC06930[-1]+1])

plt.errorbar(radius_UGC06930,velocity_UGC06930,error_velocity_UGC06930,fmt='ro',lw=2)

plt.plot(radius_UGC06930,vgas_UGC06930(radius_UGC06930),'y',lw=2)

plt.plot(radius_UGC06930,np.sqrt(Y_UGC06930_sol*vdisk_UGC06930(radius_UGC06930)**2),'g',lw=2)

#no bulge

plt.plot(radius_UGC06930,np.sqrt(vvnfw_UGC06930(radius_UGC06930,m200_UGC06930_sol)+vvnfwinmg_UGC06930(radius_UGC06930,m200_UGC06930_sol,bsol,lsol)),'b',lw=2)

plt.plot(radius_UGC06930,np.sqrt(vv_UGC06930(radius_UGC06930,Y_UGC06930_sol,m200_UGC06930_sol,bsol,lsol)),'k',lw=2)

plt.xlabel('Radius (kpc)',fontsize=18)

plt.ylabel('Velocity (km/s)',fontsize=18)

plt.title('UGC06930')

import os

filenameUGC06930='result_UGC06930_new_emcee_testing_10_galaxy_set4.pdf'

pathUGC06930='/remote/pi213b/almeida/Python/new_emcee_testing_10_galaxy'

fullpathUGC06930=os.path.join(pathUGC06930,filenameUGC06930)

plt.savefig(fullpathUGC06930)

plt.close()

#UGC06983

plt.xlim([0,radius_UGC06983[-1]+1])

plt.errorbar(radius_UGC06983,velocity_UGC06983,error_velocity_UGC06983,fmt='ro',lw=2)

plt.plot(radius_UGC06983,vgas_UGC06983(radius_UGC06983),'y',lw=2)

plt.plot(radius_UGC06983,np.sqrt(Y_UGC06983_sol*vdisk_UGC06983(radius_UGC06983)**2),'g',lw=2)

#no bulge

plt.plot(radius_UGC06983,np.sqrt(vvnfw_UGC06983(radius_UGC06983,m200_UGC06983_sol)+vvnfwinmg_UGC06983(radius_UGC06983,m200_UGC06983_sol,bsol,lsol)),'b',lw=2)

plt.plot(radius_UGC06983,np.sqrt(vv_UGC06983(radius_UGC06983,Y_UGC06983_sol,m200_UGC06983_sol,bsol,lsol)),'k',lw=2)

plt.xlabel('Radius (kpc)',fontsize=18)

plt.ylabel('Velocity (km/s)',fontsize=18)

plt.title('UGC06983')

import os

filenameUGC06983='result_UGC06983_new_emcee_testing_10_galaxy_set4.pdf'

pathUGC06983='/remote/pi213b/almeida/Python/new_emcee_testing_10_galaxy'

fullpathUGC06983=os.path.join(pathUGC06983,filenameUGC06983)

plt.savefig(fullpathUGC06983)

plt.close()

#UGC07261

plt.xlim([0,radius_UGC07261[-1]+1])

plt.errorbar(radius_UGC07261,velocity_UGC07261,error_velocity_UGC07261,fmt='ro',lw=2)

plt.plot(radius_UGC07261,vgas_UGC07261(radius_UGC07261),'y',lw=2)

plt.plot(radius_UGC07261,np.sqrt(Y_UGC07261_sol*vdisk_UGC07261(radius_UGC07261)**2),'g',lw=2)

#no bulge

plt.plot(radius_UGC07261,np.sqrt(vvnfw_UGC07261(radius_UGC07261,m200_UGC07261_sol)+vvnfwinmg_UGC07261(radius_UGC07261,m200_UGC07261_sol,bsol,lsol)),'b',lw=2)

plt.plot(radius_UGC07261,np.sqrt(vv_UGC07261(radius_UGC07261,Y_UGC07261_sol,m200_UGC07261_sol,bsol,lsol)),'k',lw=2)

plt.xlabel('Radius (kpc)',fontsize=18)

plt.ylabel('Velocity (km/s)',fontsize=18)

plt.title('UGC07261')

import os

filenameUGC07261='result_UGC07261_new_emcee_testing_10_galaxy_set4.pdf'

pathUGC07261='/remote/pi213b/almeida/Python/new_emcee_testing_10_galaxy'

fullpathUGC07261=os.path.join(pathUGC07261,filenameUGC07261)

plt.savefig(fullpathUGC07261)

plt.close()

#UGC07690

plt.xlim([0,radius_UGC07690[-1]+1])

plt.errorbar(radius_UGC07690,velocity_UGC07690,error_velocity_UGC07690,fmt='ro',lw=2)

plt.plot(radius_UGC07690,vgas_UGC07690(radius_UGC07690),'y',lw=2)

plt.plot(radius_UGC07690,np.sqrt(Y_UGC07690_sol*vdisk_UGC07690(radius_UGC07690)**2),'g',lw=2)

#no bulge

plt.plot(radius_UGC07690,np.sqrt(vvnfw_UGC07690(radius_UGC07690,m200_UGC07690_sol)+vvnfwinmg_UGC07690(radius_UGC07690,m200_UGC07690_sol,bsol,lsol)),'b',lw=2)

plt.plot(radius_UGC07690,np.sqrt(vv_UGC07690(radius_UGC07690,Y_UGC07690_sol,m200_UGC07690_sol,bsol,lsol)),'k',lw=2)

plt.xlabel('Radius (kpc)',fontsize=18)

plt.ylabel('Velocity (km/s)',fontsize=18)

plt.title('UGC07690')

import os

filenameUGC07690='result_UGC07690_new_emcee_testing_10_galaxy_set4.pdf'

pathUGC07690='/remote/pi213b/almeida/Python/new_emcee_testing_10_galaxy'

fullpathUGC07690=os.path.join(pathUGC07690,filenameUGC07690)

plt.savefig(fullpathUGC07690)

plt.close()