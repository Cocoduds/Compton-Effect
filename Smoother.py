# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 10:47:19 2022

@author: dudle
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from Calibration import p1 as Calibrate
from scipy.signal import savgol_filter
from statistics import mean

plt.close("all")
data=np.genfromtxt('Cs-137_PTFE.txt', skip_header = 5, skip_footer = 1, autostrip = True)
data[:,0] = np.linspace(0,510,511)

def Smooth(data, bins):
    dataNew = []
    stepsize = int(len(data)/bins)
    step = 0
    for i in range(0,bins):
        currentstep=[]
        for i in range(step,step+stepsize):
            currentstep.append(data[step])
            step = step + 1
        dataNew.append(mean(currentstep))
    return (dataNew)

plt.plot(Smooth(data[:,0], 250), Smooth(data[:,1],250))