#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 23 21:31:28 2021

@author: jakubpazio
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit

#%% 
'Importing Data Section'

alldata = pd.read_csv('Experiment_Data.csv',header=6)
# print(alldata)
# type(alldata)

ALLDATA = np.array(alldata)
Bin_No = []
No_Source = []
Na_22 = []
Mn_54 = []
Cs_137 = []
Am_241 = []

for i in range(len(ALLDATA)):
    Bin_No.append(ALLDATA[i][0])
    No_Source.append(ALLDATA[i][1])
    Na_22.append(ALLDATA[i][2])
    Mn_54.append(ALLDATA[i][3])
    Cs_137.append(ALLDATA[i][4])
    Am_241.append(ALLDATA[i][5])
    
#%%
'Gaussian Fit Section'

'First we trim range Bin_No and for a specific source'
left= 170
right = 220
bin_range  =Bin_No[left:right]
Na_22_range = Na_22[left:right]
Mn_54_range = Mn_54[left:right]
Cs_137_range = Cs_137[left:right]
Am_241_range = Am_241[left:right]
nosource= No_Source[left:right]
yerr=[]
for i in range(len(Cs_137_range)):
    yerr.append(math.sqrt(Cs_137_range[i]))
    
# plt.plot(bin_range,Na_22_range,'bo',label='Sodium 22')
# plt.plot(bin_range,Mn_54_range,'bo',label='Manganese 54')
# plt.plot(bin_range,Cs_137_range,'bo',label='Caesium 137')
plt.errorbar(bin_range,Cs_137_range, yerr=yerr,label='Caesium 137', capsize=4,marker='o', linestyle="",capthick=1)
# plt.plot(bin_range,Am_241_range,'bo',label='Americium 241')
# plt.plot(bin_range,nosource,'go',label='no source')


# Fitting function
def gaussian(x, amp, cen, wid):
    return amp * np.exp(-(x-cen)**2 / wid)

initialGuess = [1200,190,15]   

# popt, pcov = curve_fit(gaussian, bin_range , Na_22_range, initialGuess)
# popt, pcov = curve_fit(gaussian, bin_range , Mn_54_range, initialGuess)
popt, pcov = curve_fit(gaussian, bin_range , Cs_137_range, initialGuess)
# popt, pcov = curve_fit(gaussian, bin_range , Am_241_range, initialGuess)
print(popt)
print(pcov)


 
#x values for the fitted function
xFit = np.arange(left, right , 0.01)
 
#Plot the fitted function
plt.plot(xFit, gaussian(xFit, *popt), 'r', label='fit params: amp=%5.3f, cen=%5.3f , wid=%5.3f' % tuple(popt))
# plt.ylim([0, 2000])




#%%
'Plotting Section'

# plt.plot(Bin_No,No_Source,label='No Source')
# plt.plot(Bin_No,Na_22,label='Sodium 22')
# plt.plot(Bin_No,Mn_54,label='Manganese 54')
# plt.plot(Bin_No,Cs_137,label='Caesium 137')
# plt.plot(Bin_No,Am_241,label='Americium 241')
plt.xlabel('Bin Number [A.u]',fontsize="large",fontweight='demi')
plt.ylabel('Count',fontsize="large",fontweight='demi')
plt.title("Count Number vs Bin Number Cs-137 Isotope",fontweight='demi')
plt.legend()
plt.ylim(0, 1800)

plt.show()

# Bin_Number = ALLDATA[0]
# print(Bin_Number)