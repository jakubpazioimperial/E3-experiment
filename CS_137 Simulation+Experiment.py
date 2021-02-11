#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 31 00:34:36 2021

@author: jakubpazio
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit
from random import random
from scipy.stats import norm
# import linear_interpolation as LI
#%% 
'Importing Data Section'

alldata = pd.read_csv('Experiment_Data.csv',header=6)
# print(alldata)
# type(alldata)

ALLDATA = np.array(alldata)
Bin_No = []
No_Source = []
Cs_137 = []
Cs_137_deg20_S = []
Cs_137_deg20_B = []
Cs_137_deg30_S = []
Cs_137_deg30_B = []
Cs_137_deg45_S = []
Cs_137_deg45_B = []

calibration_factor=3.399

for i in range(len(ALLDATA)):
    Bin_No.append(ALLDATA[i][0]*calibration_factor)
    No_Source.append(ALLDATA[i][1])
    Cs_137.append(ALLDATA[i][4])
    Cs_137_deg20_S.append(ALLDATA[i][6])
    Cs_137_deg20_B.append(ALLDATA[i][7])
    Cs_137_deg30_S.append(ALLDATA[i][8])
    Cs_137_deg30_B.append(ALLDATA[i][9])
    Cs_137_deg45_S.append(ALLDATA[i][10])
    Cs_137_deg45_B.append(ALLDATA[i][11])


# print(Bin_No)

    
#%%
'Gaussian Fit Section'

'First we trim range Bin_No and for a specific source'
left= int(400/calibration_factor) # from 0 KeV
right = int(800/calibration_factor) # to 800 KeV
bin_range  =Bin_No[left:right]
Cs_137_range = Cs_137[left:right]
Cs_137_deg20_S_range = Cs_137_deg20_S[left:right]

# plt.plot(bin_range,Cs_137_range,'bo',label='Caesium 137') # 0 deg Cesium experiment
# plt.plot(bin_range,Cs_137_deg20_S_range,'go',label='Caesium 137 20 deg') # 20 deg Cesium experiment


# Fitting function
def gaussian(x, amp, cen, wid):
    return amp * np.exp(-(x-cen)**2 / wid)




def interpolate0deg():

    initialGuess = [1200,661,15]   
    popt, pcov = curve_fit(gaussian, bin_range , Cs_137_range, initialGuess)
    print(popt)
     
    #x values for the fitted function
    xFit = np.arange(500, 800 , 0.01)     
    #Plot the fitted function
    plt.plot(xFit, gaussian(xFit, *popt), 'r', label='fit params: amp=%5.3f, cen=%5.3f , wid=%5.3f' % tuple(popt))
    # plt.ylim([0, 2000])


# interpolate0deg()


def interpolate20deg():

    initialGuess = [1200,661,15]   
    popt, pcov = curve_fit(gaussian, bin_range , Cs_137_deg20_S_range, initialGuess)
    print(popt)
     
    #x values for the fitted function
    xFit = np.arange(500, 800 , 0.01)     
    #Plot the fitted function
    plt.plot(xFit, gaussian(xFit, *popt), 'y', label='fit params: amp=%5.3f, cen=%5.3f , wid=%5.3f' % tuple(popt))
    # plt.ylim([0, 2000])

# interpolate20deg()



#%%
'Simulation Section'


plt.xlabel('Energy [keV]',fontsize="large",fontweight='demi')
plt.ylabel('Count',fontsize="large",fontweight='demi')
plt.title("Monte Carlo Simulation",fontweight='demi')
plt.legend()
plt.ylim(0,3000)

plt.show()

# Bin_Number = ALLDATA[0]
# print(Bin_Number)

# CS137_scaling = 1.92451452 
CS137_scaling = 1.9675
bin_width = 5

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
# b20,c20 = MC_Simulation(662, 20, 0.005, 9, 5, 0.075)
# b30,c30 = MC_Simulation(662, 30, 0.005, 9, 5, 0.075)
# b45,c45 = MC_Simulation(662, 45, 0.005, 9, 5, 0.075)



'Plotting section for simulation data'

plt.plot(b0, c0, 'go', label='Simulation Scattering 0 deg')
# plt.plot(b20, c20, '-', label='Scattering 20 deg')
# plt.plot(b30, c30, '-', label='Scattering 30 deg')
# plt.plot(b45, c45, '-', label='Scattering 45 deg')
 
# plt.xlabel('Bin number')
# plt.legend()
# plt.ylabel(' Number of counts')
# plt.show()


#%%
'Linear noise section addition'


def Noise_add(bins):
    y_noise = []
    y_intercept=800
    slope = -10
    for i in range(len(bins)):
        y_noise.append(y_intercept + slope * i)
    return y_noise
    
y_noise = Noise_add(b0)

plt.plot(b0, y_noise, 'r-', label='Noise')
plt.plot(b0,y_noise+c0, 'yo', label='Simulation Scattering 0 deg + Noise')




#%%


'First we trim range Bin_No and for a specific source'
'Section interpolates gausian fit into data'

def gausian_fit(left,right,Bin_No,Count,initialGuess):

    bin_range  =Bin_No[left:right]
    count_range = Count[left:right]
    

    popt, pcov = curve_fit(gaussian, bin_range , count_range, initialGuess)

    # print(popt)
    #x values for the fitted function
    xFit = np.arange(bin_range[0],bin_range[-1] , 0.01)
    # print('xfit',xFit)
     #Plot the fitted function
    print(pcov)

    error = math.sqrt(pcov[1][1])
    return xFit, gaussian(xFit, *popt), popt[1] ,error
    


l=60
r=75

# a,b,c=gausian_fit(l,r,b0,c0,[3000,300*2,50])
# gausian_fit(l,r,b20,c20,[3000,300*2,50])
# gausian_fit(l,r,b30,c30,[3000,300*2,50])
# gausian_fit(l,r,b45,c45,[3000,250*2,50])


# plt.plot(a, b, 'r', label='fit params: amp=%5.3f, cen=%5.3f , wid=%5.3f' % c)
# plt.plot(b0[l:r],c0[l:r],'yo',label='MC_simulatiom')


a1,b1,c1, error1=gausian_fit(l,r,b0,c0,[3000,300*2,50])
a2,b2,c2, error2 =gausian_fit(l,r,b0,y_noise+c0,[3000,300*2,50])
print('ERROR 1', error1)

plt.plot(a1, b1, 'g', label='fit params: cen=%5.3f $\pm$ %5.3f ' % (c1,error1) )
plt.plot(a2, b2, 'y', label='fit params: cen=%5.3f $\pm$ %5.3f,' % (c2,error2))
# plt.plot(b0[l:r],c0[l:r],'yo',label='MC_simulatiom')

plt.legend(loc='upper right')
plt.show()




#%%
# 'Parameters, change these to optimize analysis'

# h = 5   #This is the smoothing parameter, it should be an odd number
# grad = LI.popt[0]
# grad_err = LI.pcov[0][0]

# Range_20 = [500,670]            #Trims the energy range for each Gaussian 
# Range_30 = [387,680]
# Range_45 = [330,600]

# 'Estimate the Initial Parameters of each Gaussian'
# Params_20 = [7.5,570,70]        #[Amplitude,Mean,Width]
# Params_30 = [6,540,120]
# Params_45 = [14.5,500,100]

# 'This is the Energy of Cs_137 Gamma Photons without deflector'
# E0 = 192*grad

# Upper_Angle = 60                #Max value of angle plotted to for Compton Scattering
# Upper_Energy = 1000             #Max value of energy plotted to for Compton Scattering
# #%% 
# 'Importing Data Section'
# alldata = pd.read_csv('Experiment_Data.csv',header=6)
# # print(alldata)
# # type(alldata)
# ALLDATA = np.array(alldata)
# Bin_No = []

# Cs_137_deg20_S = []
# Cs_137_deg20_B = []
# Cs_137_deg30_S = []
# Cs_137_deg30_B = []
# Cs_137_deg45_S = []
# Cs_137_deg45_B = []

# for i in range(len(ALLDATA)):
#     Bin_No.append(ALLDATA[i][0])
#     Cs_137_deg20_S.append(ALLDATA[i][6])
#     Cs_137_deg20_B.append(ALLDATA[i][7])
#     Cs_137_deg30_S.append(ALLDATA[i][8])
#     Cs_137_deg30_B.append(ALLDATA[i][9])
#     Cs_137_deg45_S.append(ALLDATA[i][10])
#     Cs_137_deg45_B.append(ALLDATA[i][11])
    
# Energy = [i*grad for i in Bin_No]                            #A list of the energy of each bin so we can plot Count against Energy




