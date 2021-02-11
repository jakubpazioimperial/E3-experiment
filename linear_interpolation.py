#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 24 19:51:43 2021

@author: jakubpazio
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit

# Fitting function
def func(x, a):
    # return a*np.exp(b*x) + c
    return a * x 
    #return a*x+b
 
#Experimental x and y data points    
xData = np.array([149, 356, 238, 191,10.55,21.72])
yData = np.array([511, 1200, 834, 662, 26.3, 59.6])
 
#Plot experimental data points
plt.plot(xData, yData, 'bo', label='experimental-data')
 
# Initial guess for the parameters
initialGuess = [1.0]    
 
#Perform the curve-fit
popt, pcov = curve_fit(func, xData, yData, initialGuess)
print(popt)
 
#x values for the fitted function
xFit = np.arange(0.0, 400.0, 0.01)
 
#Plot the fitted function


plt.plot(xFit, func(xFit, *popt), 'r', label='fit params: a=%5.3f' % tuple(popt))
 
plt.xlabel('Bin Number')
plt.ylabel(' Energy [keV]')
plt.legend()
plt.show()