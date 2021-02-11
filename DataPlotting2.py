#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 12:41:39 2021

@author: robertoking
"""
'Callable Functions: PlotData PlotAll PlotGaussiansOnly PlotAngles1 PlotAngles2 PlotCompton'

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
# import linear_interpolation as LI
# from linear_interpolation import pcov,popt

#%%
'Parameters, change these to optimize analysis'

h = 3   #This is the smoothing parameter, it should be an odd number
grad = 3.403
# grad_err = LI.pcov[0][0]

Range_0 = [550,700]
Range_20 = [500,670]            #Trims the energy range for each Gaussian 
Range_30 = [387,680]
Range_45 = [330,600]

'Estimate the Initial Parameters of each Gaussian'
Params_0 = [7.5,670,70]  
Params_20 = [7.5,570,70]        #[Amplitude,Mean,Width]
Params_30 = [6,540,120]
Params_45 = [14.5,500,100]

'This is the Energy of Cs_137 Gamma Photons without deflector'
E0 = 662

Upper_Angle = 60                #Max value of angle plotted to for Compton Scattering
Upper_Energy = 1000             #Max value of energy plotted to for Compton Scattering
#%% 
'Importing Data Section'
alldata = pd.read_csv('Experiment_Data.csv',header=6)
# print(alldata)
# type(alldata)
ALLDATA = np.array(alldata)
Bin_No = []
Cs_137 = []
Cs_137_deg20_S = []
Cs_137_deg20_B = []
Cs_137_deg30_S = []
Cs_137_deg30_B = []
Cs_137_deg45_S = []
Cs_137_deg45_B = []

for i in range(len(ALLDATA)):
    Bin_No.append(ALLDATA[i][0])
    Cs_137.append(ALLDATA[i][4])
    Cs_137_deg20_S.append(ALLDATA[i][6])
    Cs_137_deg20_B.append(ALLDATA[i][7])
    Cs_137_deg30_S.append(ALLDATA[i][8])
    Cs_137_deg30_B.append(ALLDATA[i][9])
    Cs_137_deg45_S.append(ALLDATA[i][10])
    Cs_137_deg45_B.append(ALLDATA[i][11])
    
Energy = [i*grad for i in Bin_No]                            #A list of the energy of each bin so we can plot Count against Energy

#%%
'Smoothing Section'
'We take the average over the nearest h bins, eg 7 bins means 3 bins to the left and 3 to the right'

def smoothing(h,Data):
    Smoothed_Data = []                          #Smoothed_Data will be (h-1) shorter than Data
    for i in range(h//2,len(Data)-h//2):
        Smoothed_Data.append(sum(Data[i-h//2:i+h//2])/len(Data[i-h//2:i+h//2]))   #Could also replace len() with simply h
    Final_Smoothing = [0]*(h//2) + Smoothed_Data + [0]*(h//2)         #Need to add 3 zeros at the start of the list and 3 at the end
    return Data                                       #Final_Smoothing will be same length as Original Data

#%%
'Now time to smooth out the data and set it to new Variables'

deg_0_S = np.array(smoothing(h,Cs_137))
deg_20_S = np.array(smoothing(h,Cs_137_deg20_S))
deg_20_B = np.array(smoothing(h,Cs_137_deg20_B))
deg_30_S = np.array(smoothing(h,Cs_137_deg30_S))
deg_30_B = np.array(smoothing(h,Cs_137_deg30_B))
deg_45_S = np.array(smoothing(h,Cs_137_deg45_S))
deg_45_B = np.array(smoothing(h,Cs_137_deg45_B))

'Take Background Data away from Scattering data'

deg_0 = deg_0_S
deg_20 = deg_20_S - deg_20_B
deg_30 = deg_30_S - deg_30_B
deg_45 = deg_45_S - deg_45_B

#%%
'Fitting Gaussians Section'
def gaussian(x, amp, cen, wid):
    return amp * np.exp(-(x-cen)**2 / wid)

'Trim each data set for x & y to the range the Gaussian will be fit over'
Energy_Range_0 = Energy[int(Range_0[0]/grad) : int(Range_0[1]/grad)]
Energy_Range_20 = Energy[int(Range_20[0]/grad) : int(Range_20[1]/grad)]
Energy_Range_30 = Energy[int(Range_30[0]/grad) : int(Range_30[1]/grad)]
Energy_Range_45 = Energy[int(Range_45[0]/grad) : int(Range_45[1]/grad)]

Range_deg_0 = deg_0[int(Range_0[0]/grad) : int(Range_0[1]/grad)]
Range_deg_20 = deg_20[int(Range_20[0]/grad) : int(Range_20[1]/grad)]
Range_deg_30 = deg_30[int(Range_30[0]/grad) : int(Range_30[1]/grad)]
Range_deg_45 = deg_45[int(Range_45[0]/grad) : int(Range_45[1]/grad)]

'Optimise the initial estimated parameters from cell 1'
'[Amplitude, Mean, Width]'
Popt_0 , Pcov_0 = curve_fit(gaussian,Energy_Range_0,Range_deg_0,p0=Params_0)
Popt_20 , Pcov_20 = curve_fit(gaussian,Energy_Range_20,Range_deg_20,p0=Params_20)
Popt_30 , Pcov_30 = curve_fit(gaussian,Energy_Range_30,Range_deg_30,p0=Params_30)
Popt_45 , Pcov_45 = curve_fit(gaussian,Energy_Range_45,Range_deg_45,p0=Params_45)

'Create ranges in Energy with small steps to plot smooth Gaussians'
Energy_Fit_0 = np.arange(Range_0[0],Range_0[1],step=0.01)
Energy_Fit_20 = np.arange(Range_20[0],Range_20[1],step=0.01)
Energy_Fit_30 = np.arange(Range_30[0],Range_30[1],step=0.01)
Energy_Fit_45 = np.arange(Range_45[0],Range_45[1],step=0.01)

'Plot the Fitted Functions'

#%%
'Define Functions to plot Gaussians or Full Data'

'This function plots only the smoothed functions for each set of Scattering Data'
def PlotData():
    plt.plot(Energy,deg_20,label='20 degrees')
    plt.plot(Energy,deg_30,label='30 degrees')
    plt.plot(Energy,deg_45,label='45 degrees')

    plt.title('Cs137 Gamma Ray Count vs Energy', weight = 'bold')
    plt.xlabel('Energy / keV', weight = 'semibold')
    plt.ylabel('Gamma Ray Count', weight = 'semibold')
    plt.legend(title = 'Scattering Angle',title_fontsize=15)
    
    plt.show()
    
    
    
    
    

'Plots the full data as above and plots Fitted Gaussians over the top'
def PlotAll():
    plt.plot(Energy,deg_20,linewidth=0.5,color='blue')
    plt.plot(Energy,deg_30,linewidth=0.5,color='orange')
    plt.plot(Energy,deg_45,linewidth=0.5,color='green')
    
    plt.plot(Energy_Fit_20,gaussian(Energy_Fit_20,*Popt_20),color = 'blue',label='20 degrees, Mean = %5.1f keV' %Popt_20[1],linewidth=2)
    plt.plot(Energy_Fit_30,gaussian(Energy_Fit_30,*Popt_30),color = 'orange',label='30 degrees, Mean = %5.1f keV' %Popt_30[1],linewidth=2)
    plt.plot(Energy_Fit_45,gaussian(Energy_Fit_45,*Popt_45),color = 'green',label='45 degrees, Mean = %5.1f keV' %Popt_45[1],linewidth=2)
    
    plt.title('Cs137 Gamma Ray Count vs Energy', weight = 'bold')
    plt.xlabel('Energy / keV', weight = 'semibold')
    plt.ylabel('Gamma Ray Count', weight = 'semibold')
    plt.legend(title = 'Scattering Angle',title_fontsize=15)
    plt.show()
    
    
    
    
    

'Shows the Gaussians over the Trimmed Range'
def PlotGaussiansOnly():
    plt.plot(Energy_Fit_0,gaussian(Energy_Fit_0,*Popt_0),color = 'blue',label='0 degrees, Mean = %5.1f ± 10.2 keV ' %Popt_0[1])
    plt.plot(Energy_Range_0,Range_deg_0,'o',color='blue')   

    plt.plot(Energy_Fit_20,gaussian(Energy_Fit_20,*Popt_20),color = 'blue',label='20 degrees, Mean = %5.1f ± 9.32 keV' %Popt_20[1])
    plt.plot(Energy_Range_20,Range_deg_20,'o',color='blue')
    
    plt.plot(Energy_Fit_30,gaussian(Energy_Fit_30,*Popt_30),color = 'orange',label='30 degrees, Mean = %5.1f ± 8.53 keV' %Popt_30[1])
    plt.plot(Energy_Range_30,Range_deg_30,'o',color = 'orange')

    plt.plot(Energy_Fit_45,gaussian(Energy_Fit_45,*Popt_45),color = 'green',label='45 degrees, Mean = %5.1f ± 7.61 keV' %Popt_45[1])
    plt.plot(Energy_Range_45,Range_deg_45,'o',color = 'green')
    
    plt.title('Fitted Gaussians for different Angles of Scattering',weight='bold')
    plt.xlabel('Energy / keV', weight = 'semibold')
    plt.ylabel('Gamma Ray Count', weight = 'semibold')
    plt.legend(title='Scattering Angle',title_fontsize=15,fontsize='small')
    plt.show()


#%%
'This last section looks at the relationship between scattering angle and mean energy'

Mean_Energy = [Popt_0[1],Popt_20[1],Popt_30[1],Popt_45[1]]
Angles = [20,30,45]

def PlotAngles1():
    plt.plot(Angles,Mean_Energy,'o',color='red',label='Experimental Data')
    plt.title('Gamma Ray Energy vs Scattering Angle', weight='bold')
    plt.legend()
    plt.xlim(-5,Upper_Angle)
    plt.ylim(0,Upper_Energy)
    plt.xlabel('Scattering Angle/Degrees')
    plt.ylabel('Energy of Deflected Photons/keV')

def PlotAngles2():
    plt.errorbar(0, Mean_Energy[0], yerr=10.2, xerr=2, capsize=3, marker='o', linestyle="",capthick=1,label='0°, Energy = %5.1f ± 10.2 keV' %Mean_Energy[0])
    plt.errorbar(Angles[0],Mean_Energy[1], yerr=9.3, xerr=2, capsize=3, marker='o', linestyle="",capthick=1,color='blue',label='20°, Energy = %5.1f ± 9.3 keV' %Mean_Energy[1])
    plt.errorbar(Angles[1],Mean_Energy[2], yerr=8.53, xerr=2, capsize=3, marker='o', linestyle="",capthick=1,color='orange',label='30°, Energy = %5.1f ± 8.5 keV' %Mean_Energy[2])
    plt.errorbar(Angles[2],Mean_Energy[3], yerr=7.61, xerr=2, capsize=3, marker='o', linestyle="",capthick=1,color='green',label='45°, Energy = %5.1f ± 7.6 keV' %Mean_Energy[3])
    plt.title('Gamma Ray Energy vs Scattering Angle', weight='bold')
    plt.legend()
    plt.xlim(-5,Upper_Angle)
    plt.ylim(400,800)
    plt.xlabel('Scattering Angle [Degrees]')
    plt.ylabel('Energy of Scattered Photons [keV]')

'Define Compton Scattering Formula'
def Scatter(E0,θ):
    return E0/(1+(E0/511)*(1-np.cos(θ*np.pi/180)))

def PlotCompton():
    xvals = np.arange(0,Upper_Angle,step = 0.01)
    yvals = Scatter(E0,xvals)
    plt.plot(xvals,yvals,color='blue',label = 'Final-state Photon Energy')
    plt.title('Gamma Ray Energy vs Scattering Angle', weight='bold')
    plt.legend()
    plt.xlim(-5,Upper_Angle)
    plt.ylim(400,760)

def Scatter2(E0,θ):
    return 0.9948*E0/(1+(E0/511)*(1-np.cos(θ*np.pi/180)))-1.8
def Scatter3(E0,θ):
    return 0.9936*E0/(1+(E0/511)*(1-np.cos(θ*np.pi/180)))-4.6
def Scatter4(E0,θ):
    return 0.996*E0/(1+(E0/511)*(1-np.cos(θ*np.pi/180)))+0.1


def PlotCompton2():
    xvals = np.arange(0,Upper_Angle,step = 0.01)
    yvals2 = Scatter2(E0,xvals)
    yvals3 = Scatter3(E0,xvals)
    yvals4 = Scatter4(E0,xvals)
    plt.plot(xvals,yvals2,'-',color='orange',label = 'Final-state Photon Energy with correction')
    plt.fill_between(xvals, yvals3, yvals4, alpha=0.3, color='orange')
    # plt.plot(xvals,yvals,color='blue',label = 'Final-state Photon Energy')
    plt.title('Gamma Ray Energy vs Scattering Angle', weight='bold')
    plt.legend()
    plt.xlim(-5,Upper_Angle)
    plt.ylim(440,740)
    plt.show()




PlotAngles2()
PlotCompton()
PlotCompton2()







