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
data120 = np.genfromtxt('Cs-137_120deg.txt', skip_header = 5, skip_footer = 1, autostrip = True)
data45 = np.genfromtxt('Cs-137_45deg(andBKR).txt', skip_header = 5, skip_footer = 1, autostrip = True)
data[:,6]= noisedata[:,1]
peaks=[]
error=[]

#REMOVING BACKGORUND NOISE
for i in range(1,6):
    data[:,i] = data[:,i]-noisedata[:,i+3]
#NOISE DATA FROM 0 DEG IS BIGGER THAN DATA
#data[:,6] = data [:,6]-noisedata[:,2]
data120[:,1] = data120[:,1]-data120[:,2]
data45[:,1] = data45[:,1]-data45[:,2]



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
deg=[30,60,70,80,90]

for i in range(0,5):
    deg[i-1] = deg[i-1]*np.pi/180             
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
    err_count=np.sqrt(np.abs(data[:,i])+1)
    popt_2gauss, pcov_2gauss = curve_fit(_2gaussian, bins, data[:,i], p0=[30, 20, 10, 30, 250, 10], sigma=err_count, absolute_sigma=True)
    perr_2gauss = np.sqrt(np.diag(pcov_2gauss))
    pars_1 = popt_2gauss[0:3]
    pars_2 = popt_2gauss[3:6]
    gauss_peak_1 = _1gaussian(bins, *pars_1)
    gauss_peak_2 = _1gaussian(bins, *pars_2)
    plt.plot(Calibrate(bins), _2gaussian(bins, *popt_2gauss),label=str(deg[1-i])+" deg")
    
    
    #GRAPHICS
    plt.plot(Calibrate(bins), data[:,i],label=str(deg[1-i])+" deg")
    plt.xlabel('Energy (KeV)')
    plt.ylabel('Counts')
    plt.legend()
    ax = plt.gca()
    # ax.set_ylim([0, 80])
    peaks.append(Calibrate(pars_2[1]))
    error.append(Calibrate(np.sqrt(perr_2gauss[4])))


#0           THIS IS POINTLESS BECAUSE YOU CANT SHEILD
deg.append(0)
plt.figure(6)
err_count=[]
for i in data[:,6]:
    if i<20:
        err_count.append(np.sqrt(np.abs(i)+1))
    else:
        err_count.append(np.sqrt(np.abs(i)))
popt_2gauss, pcov_2gauss = curve_fit(_2gaussian, bins, data[:,6], p0=[500, 250, 250, 2500, 661, 50], sigma=err_count, absolute_sigma=True, maxfev = 500000)
perr_2gauss = np.sqrt(np.diag(pcov_2gauss))
pars_1 = popt_2gauss[0:3]
pars_2 = popt_2gauss[3:6]
gauss_peak_1 = _1gaussian(bins, *pars_1)
gauss_peak_2 = _1gaussian(bins, *pars_2)
plt.plot(Calibrate(bins), _2gaussian(bins, *popt_2gauss),label="0 deg")
#GRAPHICS
plt.plot(Calibrate(bins), data[:,6],label="0 deg")
plt.xlabel('Energy (KeV)')
plt.ylabel('Counts')
plt.legend()
ax = plt.gca()
# ax.set_ylim([0, 80])
peaks.append(Calibrate(pars_1[1]))
error.append(Calibrate(np.sqrt(perr_2gauss[4])))


#120 Degrees
deg.append(120*np.pi/180)
plt.figure(7)
err_count=np.sqrt(np.abs(data120[:,1])+1)
popt_2gauss, pcov_2gauss = curve_fit(_2gaussian, bins, data120[:,1], p0=[500, 250, 250, 2500, 661, 50], sigma=err_count, absolute_sigma=True)
perr_2gauss = np.sqrt(np.diag(pcov_2gauss))
pars_1 = popt_2gauss[0:3]
pars_2 = popt_2gauss[3:6]
gauss_peak_1 = _1gaussian(bins, *pars_1)
gauss_peak_2 = _1gaussian(bins, *pars_2)
plt.plot(Calibrate(bins), _2gaussian(bins, *popt_2gauss),label="120 deg")
#GRAPHICS
plt.plot(Calibrate(bins), data120[:,1],label="120 deg")
plt.xlabel('Energy (KeV)')
plt.ylabel('Counts')
plt.legend()
ax = plt.gca()
# ax.set_ylim([0, 80])
peaks.append(Calibrate(pars_1[1]))
error.append(Calibrate(np.sqrt(perr_2gauss[1])))


#45 Degrees
deg.append(45*np.pi/180)
plt.figure(8)
err_count=np.sqrt(np.abs(data45[:,1])+1)
popt_2gauss, pcov_2gauss = curve_fit(_2gaussian, bins, data45[:,1], p0=[30, 20, 10, 30, 250, 10], sigma=err_count, absolute_sigma=True)
perr_2gauss = np.sqrt(np.diag(pcov_2gauss))
pars_1 = popt_2gauss[0:3]
pars_2 = popt_2gauss[3:6]
gauss_peak_1 = _1gaussian(bins, *pars_1)
gauss_peak_2 = _1gaussian(bins, *pars_2)
plt.plot(Calibrate(bins), _2gaussian(bins, *popt_2gauss),label="45 deg")
#GRAPHICS
plt.plot(Calibrate(bins), data45[:,1],label="45 deg")
plt.xlabel('Energy (KeV)')
plt.ylabel('Counts')
plt.legend()
ax = plt.gca()
# ax.set_ylim([0, 80])
peaks.append(Calibrate(pars_2[1]))
error.append(Calibrate(np.sqrt(perr_2gauss[4])))


plt.figure(9)    
plt.errorbar(deg,peaks, yerr = error, xerr = 5*np.pi/180,  color = 'blue', ls='none',label = 'Data', fmt='.')
x=np.linspace(0,np.pi, 1000)
plt.plot(x, Compton(661,x), label = 'Expected', color = 'orange')
plt.legend()
plt.xlabel('Radians')
plt.ylabel('Energy(KeV)')

def Deg():
    return(deg)
    
def WavelengthRatio(n):
    return(peaks[i]/661.7)
