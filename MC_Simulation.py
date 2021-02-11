#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 00:09:11 2021

@author: jakubpazio
"""

from scipy.stats import norm
import math
from random import seed
import numpy as np
from random import random
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#  Source
Photon_en = 662 #(keV)

# Target
Electron_r_mass = 511 #(keV)
Scattering_an = 0 #(deg)

# Detector
Gain = 0.005 # (volts per keV)
MCA_b_d = 9 # bit depth
MCA_m_s = 5 #max signal  (volts)
Detector_r = 0.075 #resolution  at 662 keV 


Scattering_Photon_e = Photon_en/(1 + (Photon_en/Electron_r_mass)*(1-math.cos(math.pi*Scattering_an/180)))
# print('Scattering_Photon_e',Scattering_Photon_e)
Analog_signal = Scattering_Photon_e * Gain
# print('Analog_signal',Analog_signal)

sc = Detector_r * Analog_signal / 2



#%%
' Random seed and Noise section'

rand = []
noise = []
Outcome_bin = []
N = 20000 # number of itterations
for i in range(N):
    rand.append(random())
    noise.append(norm.ppf(rand[i], loc=0, scale=sc))
    Outcome_bin.append(math.floor((2**MCA_b_d)*(Analog_signal + noise[i])/MCA_m_s))

'NORM.INV is norm.ppf'

# print(sum(Outcome_bin))

# print('rand',rand)
# print('noise',noise)
# print('Outcome_bin',Outcome_bin)




#%%
'Bins and Count'

bin_width = 5
bins=[]
range_bin = 600
Count = np.zeros(int(range_bin/bin_width))

for j in range(int(range_bin/bin_width)):
    bins.append(j*bin_width)
    Count_1=0
    Count_2=0
    for i in range(len(Outcome_bin)):

        if Outcome_bin[i] < bins[j] + bin_width:
            Count_1=Count_1+1
        if Outcome_bin[i] < bins[j]:
            Count_2=Count_2+1
    # print(Count_1 - Count_2)
    Count[j] = Count_1 - Count_2
            
print('Count sum',sum(Count))
# print(Count)

# plt.plot(bins, Count, 'bo')
 
# plt.xlabel('Number of counts')
# plt.ylabel(' Bin number')
# plt.show()

#%%





CS137_scaling = 1.9672141

def MC_Simulation(Photon_en, Scattering_an, Gain, MCA_b_d, MCA_m_s, Detector_r):

    # Target
    Electron_r_mass = 511 #(keV)
    
    
    Scattering_Photon_e = Photon_en/(1 + (Photon_en/Electron_r_mass)*(1-math.cos(math.pi*Scattering_an/180)))
    Analog_signal = Scattering_Photon_e * Gain
    
    sc = Detector_r * Analog_signal / 2
    
    rand = []
    noise = []
    Outcome_bin = []
    N = 8000 # number of itterations
    
    for i in range(N):
        rand.append(random())
        noise.append(norm.ppf(rand[i], loc=0, scale=sc))
        Outcome_bin.append(math.floor((2**MCA_b_d)*(Analog_signal + noise[i])/MCA_m_s))

    bin_width = 5
    bins=[]
    range_bin = int(800/CS137_scaling)
    Count = np.zeros(int(range_bin/bin_width))
    
    for j in range(int(range_bin/bin_width)):
        bins.append(j*bin_width)
        Count_1=0
        Count_2=0
        for i in range(len(Outcome_bin)):
    
            if Outcome_bin[i] < bins[j] + bin_width:
                Count_1=Count_1+1
            if Outcome_bin[i] < bins[j]:
                Count_2=Count_2+1
        # print(Count_1 - Count_2)
        Count[j] = Count_1 - Count_2
    # print('Count sum',sum(Count))
    print('Outcome_bin sum',sum(Outcome_bin))
    for i in range(len(bins)):
        bins[i]= bins[i]*CS137_scaling 
    return bins, Count




b0,c0 = MC_Simulation(662, 0, 0.005, 9, 5, 0.075)
b20,c20 = MC_Simulation(662, 20, 0.005, 9, 5, 0.075)
b30,c30 = MC_Simulation(662, 30, 0.005, 9, 5, 0.075)
b45,c45 = MC_Simulation(662, 45, 0.005, 9, 5, 0.075)

plt.plot(b0, c0, '-', label='Scattering 0 deg')
plt.plot(b20, c20, '-', label='Scattering 20 deg')
plt.plot(b30, c30, '-', label='Scattering 30 deg')
plt.plot(b45, c45, '-', label='Scattering 45 deg')
 
plt.xlabel('Bin number')
plt.legend()
plt.ylabel(' Number of counts')
plt.show()





#%%
'Gaussian Fit Section'


# Fitting function
def gaussian(x, amp, cen, wid):
    return amp * np.exp(-(x-cen)**2 / wid)

bin_width = 5

'First we trim range Bin_No and for a specific source'
def gausian_fit(left,right,Bin_No,Count,initialGuess):

    bin_range  =Bin_No[left:right]
    count_range = Count[left:right]
    
    plt.plot(bin_range,count_range,'go',label='MC_simulatiom')

    popt, pcov = curve_fit(gaussian, bin_range , count_range, initialGuess)

    print(popt)
    #x values for the fitted function
    xFit = np.arange(bin_range[0],bin_range[-1] , 0.01)
    print('xfit',xFit)
     #Plot the fitted function
    plt.plot(xFit, gaussian(xFit, *popt), 'y', label='fit params: amp=%5.3f, cen=%5.3f , wid=%5.3f' % tuple(popt))
    # plt.ylim([0, 11000])
    plt.legend(loc='upper right')
    plt.show()


# left = int(400/CS137_scaling)
# right = int(700/CS137_scaling)

gausian_fit(55,100,b0,c0,[3000,650,50])


