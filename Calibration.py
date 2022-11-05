# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 14:06:39 2022

@author: dudle
"""


import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

plt.close("all")

Am=np.loadtxt("Am_251cal.csv", delimiter = ',', encoding='utf-8-sig')
Cs=np.loadtxt("Cs_137cal.csv", delimiter = ',', encoding='utf-8-sig')
Co=np.loadtxt("Co_57cal.csv", delimiter = ',', encoding='utf-8-sig')


expected = [59.6, 662, 122]
peaks = []


#x is energy
#y is counts

def fit(x,a,mu,sigma):
    y_gaus=a*np.exp(-(x-mu)**2 / (2*sigma**2))
    return y_gaus 



def FindMaxima(numbers):
    maxima = []
    length = len(numbers)
    if length >= 2:
        if numbers[0] > numbers[1]:
            maxima.append(numbers[0])
        if length > 3:
            for i in range(1, length-1):    
                if numbers[i] > numbers[i-1] and numbers[i] > numbers[i+1]:
                    maxima.append(numbers[i])
        if numbers[length-1] > numbers[length-2]:     
            maxima.append(numbers[length-1])       
    return maxima


x = np.linspace(0,511,512)


xloads= np.linspace(0,511,5110)


plt.figure(1)
mu=np.mean(x)
sigma= np.std(x)
a= 1/(sigma*np.sqrt(2*np.pi))
initial_guess=[a,mu,sigma]
gaus,gaus_cov = curve_fit(fit,x,Am,initial_guess,maxfev=5000)
plt.plot(x,Am)
y=fit(xloads,*gaus).tolist()
plt.plot(xloads,y)
peaks.append(gaus[1])
plt.show()



Cs_res=[]
xres=[]
for i in range(350,511):
    Cs_res.append(Cs[i])
    xres.append(x[i])

xloads= np.linspace(350,511,(511-350)*10)
plt.figure(2)
gaus,gaus_cov = curve_fit(fit,xres,Cs_res,initial_guess,maxfev=5000)
plt.plot(x,Cs)
y=fit(xloads,*gaus).tolist()
plt.plot(xloads,y)
peaks.append(xloads[y.index(max(y))])
plt.show()

xloads= np.linspace(0,511,5110)
plt.figure(3)
gaus,gaus_cov = curve_fit(fit,x,Co,initial_guess,maxfev=5000)
plt.plot(x,Co)
y=fit(xloads,*gaus).tolist()
plt.plot(xloads,y)
peaks.append(xloads[y.index(max(y))])
plt.show()


plt.figure(4)
plt.scatter(peaks,expected)
fit,cov = np.polyfit(peaks,expected,1,cov=True)
sig = np.sqrt(cov[0,0])

p1=np.poly1d(fit)
print('Fit polynomial (for corr. fac.)')
print(p1)
plt.plot(expected,p1(expected))
