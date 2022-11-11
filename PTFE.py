# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 09:21:05 2022

@author: dudle
"""
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from Calibration import p1 as Calibrate
from scipy.signal import savgol_filter
from Smoother import Smooth

plt.close("all")
data=np.genfromtxt('Cs-137_PTFE.txt', skip_header = 5, skip_footer = 1, autostrip = True)
peaks=[]
error=[]


binsnumber = 100




def _1gaussian(x, amp1, cen1, sigma1):
    return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen1)/sigma1)**2)))


def _2gaussian(x, amp1,cen1,sigma1, amp2,cen2,sigma2):
    return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen1)/sigma1)**2))) + \
            amp2*(1/(sigma2*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen2)/sigma2)**2)))

def Compton(energy,x):
    return(energy/(1+energy/(511)*(1-np.cos(x))))

#plotting
data[:,0] = np.linspace(0,510,511)
bins = data[:,0]
deg=[30,60,90,120]
heightguess=[30,25,35,40]
p1guess = [125,100,75,50]
hg1=[15,10,15,15]  

for i in range(0,4):
    deg[i-1] = deg[i-1]*np.pi/180             
for i in range(1,5):
    #seperate figures
    plt.figure(i)
    # noise reduction
    data_smoothed = savgol_filter(data[:,2*i-1], 51, 2)
    
    # SINGLE GAUSSIAN FIT
    # popt_gauss, pcov_gauss = curve_fit(_1gaussian, bins, data[:,2*i-1], p0=[heightguess[i-1], Compton(661,deg[i-1]), 150])
    # perr_gauss = np.sqrt(np.diag(pcov_gauss))
    # plt.plot(Calibrate(bins), _1gaussian(bins, *popt_gauss),label=str(deg[i-1])+" deg")
    
    # TWO GAUSSIAN FITS
    err_count=np.sqrt(np.abs(Smooth(data[:,2*i-1]-data[:,2*i], binsnumber))+1)
    popt_2gauss, pcov_2gauss = curve_fit(_2gaussian, Smooth(bins,binsnumber), Smooth(data[:,2*i-1]-data[:,2*i], binsnumber), p0=[hg1[i-1], p1guess[i-1], 100, heightguess[i-1], (Compton(661,deg[i-1])+13)/1.6, 150], maxfev=50000)
    perr_2gauss = np.sqrt(np.diag(pcov_2gauss))
    pars_1 = popt_2gauss[0:3]
    pars_2 = popt_2gauss[3:6]
    gauss_peak_1 = _1gaussian(bins, *pars_1)
    gauss_peak_2 = _1gaussian(bins, *pars_2)
    plt.plot(Smooth(Calibrate(bins),binsnumber), _2gaussian(Smooth(bins,binsnumber), *popt_2gauss),label=str(deg[i-1])+" deg")
    
    
    #GRAPHICS
    plt.plot(Smooth(Calibrate(bins),binsnumber), Smooth(data[:,2*i-1]-data[:,2*i], binsnumber),label=str(deg[i-1])+" deg")
    plt.xlabel('Energy (KeV)')
    plt.ylabel('Counts')
    plt.legend()
    ax = plt.gca()
    # ax.set_ylim([0, 80])
    if i == 3:
        peaks.append(Calibrate(pars_1[1]))
        error.append(Calibrate(np.sqrt(perr_2gauss[1])))
    else:
        peaks.append(Calibrate(pars_2[1]))
        error.append(Calibrate(np.sqrt(perr_2gauss[4])))



plt.figure(8)    
plt.errorbar(deg,peaks, yerr = error, xerr = 5*np.pi/180,  color = 'blue', ls='none',label = 'Data', fmt='.')
x=np.linspace(0,np.pi, 1000)
plt.plot(x, Compton(661,x), label = 'Expected', color = 'orange')
plt.legend()
plt.xlabel('Radians')
plt.ylabel('Energy(KeV)')

def Deg():
    return(deg)
    
def WavelengthRatio(n):
    return(peaks[n]/661.7)