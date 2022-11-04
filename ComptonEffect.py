# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 14:16:00 2022

@author: dudle
"""
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from Calibration import p1 as Calibrate
from scipy.signal import savgol_filter

plt.close("all")
data=np.genfromtxt('Cs-137_2.txt', skip_header = 5, skip_footer = 1, autostrip = True)
noisedata = np.genfromtxt('Cs-137_norod.txt', skip_header = 5, skip_footer = 1, autostrip = True)

#REMOVING BACKGORUND NOISE
for i in range(1,6):
    data[:,i] = data[:,i]-noisedata[:,i+3]




def _1gaussian(x, amp1, cen1, sigma1):
    return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen1)/sigma1)**2)))


def _2gaussian(x, amp1,cen1,sigma1, amp2,cen2,sigma2):
    return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen1)/sigma1)**2))) + \
            amp2*(1/(sigma2*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen2)/sigma2)**2)))



#plotting
data[:,0] = np.linspace(0,510,511)
bins = data[:,0]
deg=[30,60,70,80,90]              
for i in range(1,6):
    #seperate figures
    plt.figure(i)
    # noise reduction
    data_smoothed = savgol_filter(data[:,i], 51, 2)
    
    # SINGLE GAUSSIAN FIT
    # popt_gauss, pcov_gauss = curve_fit(_1gaussian, bins, data[:,i], p0=[60, 150, 10])
    # perr_gauss = np.sqrt(np.diag(pcov_gauss))
    # plt.plot(bins, _1gaussian(bins, *popt_gauss),label=str(deg[1-i])+" deg")
    
    # TWO GAUSSIAN FITS
    popt_2gauss, pcov_2gauss = curve_fit(_2gaussian, bins, data[:,i], p0=[30, 20, 10, 30,250, 10])
    perr_2gauss = np.sqrt(np.diag(pcov_2gauss))
    pars_1 = popt_2gauss[0:3]
    pars_2 = popt_2gauss[3:6]
    gauss_peak_1 = _1gaussian(bins, *pars_1)
    gauss_peak_2 = _1gaussian(bins, *pars_2)
    plt.plot(bins, _2gaussian(bins, *popt_2gauss),label=str(deg[1-i])+" deg")
    
    
    #GRAPHICS
    plt.plot(bins, data[:,i],label=str(deg[1-i])+" deg")
    plt.xlabel('Energy (KeV)')
    plt.ylabel('Counts')
    plt.legend()
    ax = plt.gca()
    # ax.set_ylim([0, 80])

