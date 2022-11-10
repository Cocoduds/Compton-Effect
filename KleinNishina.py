# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 10:40:50 2022

@author: dudle
"""
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from Calibration import p1 as Calibrate
from scipy.signal import savgol_filter
import ComptonEffect
from ComptonEffect import _1gaussian
from ComptonEffect import _2gaussian
#import os
#os.system('cls')

electronRadius = 7.94*10**(-30)
bins = np.linspace(0,510,511)

def DFC(eR, wavelengthRatio, angle):
    return 0.5*eR**2*wavelengthRatio**2*(wavelengthRatio+(1/wavelengthRatio)-2*(np.sin(angle))**2)

angles = ComptonEffect.Deg()
crossSections = []

for i in range(0,8):
    crossSections.append(DFC(electronRadius, ComptonEffect.WavelengthRatio(i), angles[i]))




data90_1 = np.genfromtxt('Cs-137_90deg_1.txt', skip_header = 5, skip_footer = 1, autostrip = True)
data90_1[:,1] = data90_1[:,1]-data90_1[:,2]
plt.figure(11)
popt_2gauss, pcov_2gauss = curve_fit(_2gaussian, bins, data90_1[:,1], p0=[30, 20, 10, 30, 250, 10])
perr_2gauss = np.sqrt(np.diag(pcov_2gauss))
pars_1 = popt_2gauss[0:3]
pars_2 = popt_2gauss[3:6]
gauss_peak_1 = _1gaussian(bins, *pars_1)
gauss_peak_2 = _1gaussian(bins, *pars_2)
plt.plot(Calibrate(bins), _2gaussian(bins, *popt_2gauss),label="45 deg")
#GRAPHICS
plt.plot(Calibrate(bins), data90_1[:,1],label="90 deg")
plt.xlabel('Energy (KeV)')
plt.ylabel('Counts')
plt.legend()
ax = plt.gca()
# ax.set_ylim([0, 80])
wR90_1=pars_2[1]/661.7
#error.append(Calibrate(perr_2gauss[1]))    

data90_2 = np.genfromtxt('Cs-137_90deg_2.txt', skip_header = 5, skip_footer = 1, autostrip = True)
data90_2[:,1] = data90_2[:,1]-data90_2[:,2]
plt.figure(11)
popt_2gauss, pcov_2gauss = curve_fit(_2gaussian, bins, data90_2[:,1], p0=[30, 20, 10, 30, 250, 10])
perr_2gauss = np.sqrt(np.diag(pcov_2gauss))
pars_1 = popt_2gauss[0:3]
pars_2 = popt_2gauss[3:6]
gauss_peak_1 = _1gaussian(bins, *pars_1)
gauss_peak_2 = _1gaussian(bins, *pars_2)
plt.plot(Calibrate(bins), _2gaussian(bins, *popt_2gauss),label="45 deg")
#GRAPHICS
plt.plot(Calibrate(bins), data90_2[:,1],label="90 deg")
plt.xlabel('Energy (KeV)')
plt.ylabel('Counts')
plt.legend()
ax = plt.gca()
# ax.set_ylim([0, 80])
wR90_2=pars_2[1]/661.7
#error.append(Calibrate(perr_2gauss[1]))    


plt.figure(10)
plt.scatter(90*np.pi/180, DFC(electronRadius, wR90_1, 90*np.pi/180), color = 'blue', label = 'ploarized 0') 
plt.scatter(90*np.pi/180, DFC(electronRadius, wR90_2, 90*np.pi/180), color = 'green', label = 'polarized 90')       
plt.scatter(angles, crossSections, color = 'orange', label = 'unpolarized')
plt.xlabel('Angle (Radians)')
plt.ylabel('Differential Cross Section')
plt.legend()
ax = plt.gca()
ax.set_ylim([-0.1*10e-59, 0.5*10e-59])
plt.show()

