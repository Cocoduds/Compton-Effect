# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 10:08:13 2022

@author: dudle
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from Calibration import p1 as Calibrate
from scipy.signal import savgol_filter
from scipy.stats import chisquare
from Smoother import Smooth
import ComptonEffect


plt.close("all")
data=np.genfromtxt('Cs-137_PTFE.txt', skip_header = 5, skip_footer = 1, autostrip = True)
fluxData = np.genfromtxt('Cs-137_Flux.txt', skip_header = 5, skip_footer = 1, autostrip = True)
backgroundFlux = fluxData[:,1]-fluxData[:,2]
peaks=[]
error=[]
fluxlist=[]
totalcounts = []
thickness = 1/2*np.pi*0.02013/2


 
def _1gaussian(x, amp1, cen1, sigma1):
    return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen1)/sigma1)**2)))

 
def _2gaussian(x, amp1,cen1,sigma1, amp2,cen2,sigma2):
    return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen1)/sigma1)**2))) + \
            amp2*(1/(sigma2*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen2)/sigma2)**2)))

def Compton(energy,x):
    return(energy/(1+energy/(511)*(1-np.cos(x))))

def integer(Array):
    array=[]
    for k in range(0,len(Array)):
        d=Array[k]
        array.append(int(d))
    return array

def Flux(peak,std,data):
    bins_number=2*int(std)+1
    bins_peak=np.linspace(int(peak)-int(std),int(peak)+int(std),bins_number)
    bins_peak=integer(bins_peak)
    bins_area=[]
    for j in bins_peak:
        bins_area.append(data[j])
    area=np.sum(bins_area)
    print(len(bins_area))
    return area
    
def data_fix(Data):
    data=[]
    for k in range(0,len(Data[0])):
        data.append(Data[0][k])
    return data

#plotting
data[:,0] = np.linspace(0,510,511)
bins = data[:,0]
deg=[30,60,90,120,0]
chi2=[]

for i in range(1,5):
    deg[i-1] = deg[i-1]*np.pi/180            
for i in range(1,5):
    #seperate figures
    plt.figure(i)
    # noise reduction
    data_smoothed = savgol_filter(data[:,i], 51, 2)
    currentData = data[:,2*i-1]-data[:,2*i]
   
    # SINGLE GAUSSIAN FIT
    # popt_gauss, pcov_gauss = curve_fit(_1gaussian, bins, data[:,i], p0=[60, 150, 10])
    # perr_gauss = np.sqrt(np.diag(pcov_gauss))
    # plt.plot(bins, _1gaussian(bins, *popt_gauss),label=str(deg[1-i])+" deg")
   
    # TWO GAUSSIAN FITS
    popt_2gauss, pcov_2gauss = curve_fit(_2gaussian, bins, currentData, p0=[30, 20, 10, 30, 250, 10])
    perr_2gauss = np.sqrt(np.diag(pcov_2gauss))
    pars_1 = popt_2gauss[0:3]
    pars_2 = popt_2gauss[3:6]
    gauss_peak_1 = _1gaussian(bins, *pars_1)
    gauss_peak_2 = _1gaussian(bins, *pars_2)
    plt.plot(Calibrate(bins), _2gaussian(bins, *popt_2gauss),label=str(deg[i-2])+" deg")
   
    
    #GRAPHICS
    plt.plot(Calibrate(bins), currentData,label=str(deg[i-1])+" deg")
    plt.xlabel('Energy (KeV)')
    plt.ylabel('Counts')
    plt.legend()
    ax = plt.gca()
    # ax.set_ylim([0, 80])
    bins_number=2*int(popt_2gauss[5])+1
    bins_peak=np.linspace(int(pars_2[1])-int(popt_2gauss[5]),int(pars_2[1])+int(popt_2gauss[5]),bins_number)
    data_BP=[]
    for j in bins_peak:
        data_BP.append(currentData[int(j)])
    chi2.append(chisquare(data_BP,_2gaussian(bins_peak, *popt_2gauss),ddof=bins_number-3))
    peaks.append(Calibrate(pars_2[1]))
    error.append(Calibrate(perr_2gauss[1]))
    
    dat1=np.array([currentData])
    dat2=(np.ceil(dat1)).astype(int)
    dat3=data_fix(dat2)
    peak=pars_2[1]
    std=popt_2gauss[5]
    flux=Flux(peak,std,dat3)
    print(flux)
    fluxlist.append(flux)
    totalcounts.append(sum(currentData))
    
    
#%%
#CALCULATING FLUX IN
#seperate figures
plt.figure(7)
# noise reduction
currentData = backgroundFlux
   
# SINGLE GAUSSIAN FIT
popt_gauss, pcov_gauss = curve_fit(_1gaussian, bins, currentData, p0=[6000, 400, 50])
perr_gauss = np.sqrt(np.diag(pcov_gauss))
plt.plot(Calibrate(bins), _1gaussian(bins, *popt_gauss),label="background flux")
   
   

#GRAPHICS
plt.plot(Calibrate(bins), currentData,label=str(deg[i-1])+" deg")
plt.xlabel('Energy (KeV)')
plt.ylabel('Counts')
plt.legend()
ax = plt.gca()
# ax.set_ylim([0, 80])
bins_number=2*int(popt_gauss[2])+1
bins_peak=np.linspace(int(popt_gauss[1])-int(popt_gauss[2]),int(popt_gauss[1])+int(popt_gauss[2]),bins_number)
data_BP=[]
for j in bins_peak:
    data_BP.append(currentData[int(j)])
chi2.append(chisquare(data_BP,_1gaussian(bins_peak, *popt_gauss),ddof=bins_number-3))
peaks.append(Calibrate(popt_gauss[1]))
error.append(Calibrate(perr_gauss[1]))

dat1=np.array([currentData])
dat2=(np.ceil(dat1)).astype(int)
dat3=data_fix(dat2)
peak=popt_gauss[1]
std=popt_gauss[2]
fluxin=Flux(peak,std,dat3)
print(fluxin)

 
#%%
#PLOTTING THE PEAKS

plt.figure(9)   
plt.errorbar([0,30*np.pi/180,60*np.pi/180,90*np.pi/180,120*np.pi/180],peaks, yerr = error, xerr = 5*np.pi/180,  color = 'blue', ls='none',label = 'Data', fmt='.')
x=np.linspace(0,np.pi, 1000)
plt.plot(x, Compton(661,x), label = 'Expected', color = 'orange')
plt.legend()
plt.xlabel('Radians')
plt.ylabel('Energy(KeV)')


#%%
#PLOTTING THE FLUX UNDER GAUSSIAN
xaxis = [30,60,90,120]
x=[]
for i in xaxis:
    x.append(i*np.pi/180)
plt.figure(10)
plt.scatter(x, 1/(18e28*0.0267*thickness)*(fluxlist/fluxin), label = "under gaussian")

#%%
fluxin = np.sum(currentData)
print(fluxin)
plt.scatter(x, 1/(18e28*0.0267*thickness)*(totalcounts/fluxin), label = "all data")
plt.legend()


#%%
#KLEIN NISHINA
electronRadius = 7.94*10**(-30)
bins = np.linspace(0,510,511)

def DFC(eR, wavelengthRatio, angle):
    return 0.5*eR*wavelengthRatio**2*(wavelengthRatio+(1/wavelengthRatio)-(np.sin(angle))**2)

angles = ComptonEffect.Deg()
crossSections = []

for i in range(0,8):
    crossSections.append(DFC(electronRadius, ComptonEffect.WavelengthRatio(i), angles[i]))
 


plt.figure(10)
plt.scatter(angles, crossSections, color = 'orange', label = 'unpolarized')
plt.xlabel('Angle (Radians)')
plt.ylabel('Differential Cross Section')
plt.legend()
ax = plt.gca()
ax.set_ylim([-0.1*10e-59, 1e-28])
plt.show()
#%%
fig, (ax1,ax2,ax3,ax4,ax5,ax6) = plt.subplots(6, sharex=True)
fig.suptitle('Energy spectrum of CS-17 at different scattering angles')
ax1.plot(Calibrate(bins), savgol_filter(data[:,1],51, 2), label= "30 deg")
# ax2.plot(Calibrate(bins), savgol_filter(data45[:,1],51, 2), label= "45 deg")
ax3.plot(Calibrate(bins), savgol_filter(data[:,2],51, 2), label= "60 deg") 
ax4.plot(Calibrate(bins), savgol_filter(data[:,3],51, 2), label= "70 deg")
ax5.plot(Calibrate(bins), savgol_filter(data[:,4],51, 2), label= "80 deg")
ax6.plot(Calibrate(bins), savgol_filter(data[:,5],51, 2), label= "90 deg")  
ax1.legend()
ax2.legend()
ax3.legend()
ax4.legend()
ax5.legend()
ax6.legend()
