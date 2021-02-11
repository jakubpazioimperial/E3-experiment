#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 15:13:41 2021

@author: robertoking
"""

import DataPlotting2 as DP
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import linear_interpolation as LI

n = 5    #This is the bin Amalgamation parameter, n bins summed into one
grad = LI.popt[0]
grad_err = LI.pcov[0][0]

Range_20 = [500,670]            #Trims the energy range for each Gaussian 
Range_30 = [387,680]
Range_45 = [330,600]

'Estimate the Initial Parameters of each Gaussian'
Params_20 = [7.5,570,70]        #[Amplitude,Mean,Width]
Params_30 = [6,540,120]
Params_45 = [14.5,500,100]

'This is the Energy of Cs_137 Gamma Photons without deflector'
E0 = 192*grad

Upper_Angle = 60                #Max value of angle plotted to for Compton Scattering
Upper_Energy = 1000             #Max value of energy plotted to for Compton Scattering
#%%
'Bin Consolidation Section'
'We combine n bins into one bigger bin'

def Consolidate(n,Data):                                #This function only works on y data
    Shorter_Data = []
    for i in range(len(Data)//n):
        Shorter_Data.append(sum(Data[3*i:3*i+2]))
    return Shorter_Data

'Need to select the right x values too'
def SelectX(n,Data=DP.Energy):
    Selectx = []
    for i in range(len(Data)//n):
        Selectx.append(Data[3*i+1])
    return Selectx
       

#%%
'Now time to Consolidate the Counts'
deg_20_S_ = np.array(Consolidate(n,DP.Cs_137_deg20_S))
deg_20_B_ = np.array(Consolidate(n,DP.Cs_137_deg20_B))
deg_30_S_ = np.array(Consolidate(n,DP.Cs_137_deg30_S))
deg_30_B_ = np.array(Consolidate(n,DP.Cs_137_deg30_B))
deg_45_S_ = np.array(Consolidate(n,DP.Cs_137_deg45_S))
deg_45_B_ = np.array(Consolidate(n,DP.Cs_137_deg45_B))

'Take Background Data away from Scattering data'

deg_20_ = deg_20_S_ - deg_20_B_
deg_30_ = deg_30_S_ - deg_30_B_
deg_45_ = deg_45_S_ - deg_45_B_

New_Energy = SelectX(n,Data = DP.Energy)


#%%
'Fitting Gaussians Section'
def gaussian(x, amp, cen, wid):
    return amp * np.exp(-(x-cen)**2 / wid)

'Trim each data set for x & y to the range the Gaussian will be fit over'
Energy_Range_20_ = New_Energy[int(Range_20[0]/(n*grad)) : int(Range_20[1]/(n*grad))]
Energy_Range_30_ = New_Energy[int(Range_30[0]/(n*grad)) : int(Range_30[1]/(n*grad))]
Energy_Range_45_ = New_Energy[int(Range_45[0]/(n*grad)) : int(Range_45[1]/(n*grad))]

Range_deg_20 = deg_20_[int(Range_20[0]/(n*grad)) : int(Range_20[1]/(n*grad))]
Range_deg_30 = deg_30_[int(Range_30[0]/(n*grad)) : int(Range_30[1]/(n*grad))]
Range_deg_45 = deg_45_[int(Range_45[0]/(n*grad)) : int(Range_45[1]/(n*grad))]

'Optimise the initial estimated parameters from cell 1'
'[Amplitude, Mean, Width]'
Popt_20 , Pcov_20 = curve_fit(gaussian,Energy_Range_20_,Range_deg_20,p0=Params_20)
Popt_30 , Pcov_30 = curve_fit(gaussian,Energy_Range_30_,Range_deg_30,p0=Params_30)
Popt_45 , Pcov_45 = curve_fit(gaussian,Energy_Range_45_,Range_deg_45,p0=Params_45)

'Create ranges in Energy with small steps to plot smooth Gaussians'
Energy_Fit_20 = np.arange(Range_20[0],Range_20[1],step=0.01)
Energy_Fit_30 = np.arange(Range_30[0],Range_30[1],step=0.01)
Energy_Fit_45 = np.arange(Range_45[0],Range_45[1],step=0.01)

#%%
'This Section looks at the errors on each Gaussian'
'Use Standard Deviation of the mean to work out each error: err = standdev / root(N)'

sigma_20,sigma_30,sigma_45 = np.sqrt(Popt_20[2]/2),np.sqrt(Popt_30[2]/2),np.sqrt(Popt_45[2]/2)
N_20,N_30,N_45 = len(Energy_Range_20_),len(Energy_Range_30_),len(Energy_Range_45_)
er_20,er_30,er_45 = sigma_20/np.sqrt(N_20),sigma_30/np.sqrt(N_30),sigma_45/np.sqrt(N_45)
errors = [er_20,er_30,er_45]


#%%
'Plotting Section'

'This function plots only the smoothed functions for each set of Scattering Data'
def PlotData():
    plt.plot(New_Energy,deg_20_,label='20 degrees')
    plt.plot(New_Energy,deg_30_,label='30 degrees')
    plt.plot(New_Energy,deg_45_,label='45 degrees')

    plt.title('Cs137 Gamma Ray Count vs Energy', weight = 'bold')
    plt.xlabel('Energy / keV', weight = 'semibold')
    plt.ylabel('Gamma Ray Count', weight = 'semibold')
    plt.legend(title = 'Scattering Angle',title_fontsize=15)
    
    plt.show()
    
'Shows the Gaussians over the Trimmed Range'
def PlotGaussiansOnly():

    plt.plot(Energy_Fit_20,gaussian(Energy_Fit_20,*Popt_20),color = 'blue',label='20 degrees, Mean = %5.1f keV' %Popt_20[1])
    plt.plot(Energy_Range_20_,Range_deg_20,'o',color='blue')
    
    plt.plot(Energy_Fit_30,gaussian(Energy_Fit_30,*Popt_30),color = 'orange',label='30 degrees, Mean = %5.1f keV' %Popt_30[1])
    plt.plot(Energy_Range_30_,Range_deg_30,'o',color = 'orange')

    plt.plot(Energy_Fit_45,gaussian(Energy_Fit_45,*Popt_45),color = 'green',label='45 degrees, Mean = %5.1f keV' %Popt_45[1])
    plt.plot(Energy_Range_45_,Range_deg_45,'o',color = 'green')
    
    plt.title('Fitted Gaussians for different Angles of Scattering',weight='bold')
    plt.xlabel('Energy / keV', weight = 'semibold')
    plt.ylabel('Gamma Ray Count', weight = 'semibold')
    plt.legend(title='Scattering Angle',title_fontsize=15,fontsize='small')
    plt.show()

'Plots the full data as above and plots Fitted Gaussians over the top'
def PlotAll():
    plt.plot(New_Energy,deg_20_,linewidth=0.5,color='blue')
    plt.plot(New_Energy,deg_30_,linewidth=0.5,color='orange')
    plt.plot(New_Energy,deg_45_,linewidth=0.5,color='green')
    
    plt.plot(Energy_Fit_20,gaussian(Energy_Fit_20,*Popt_20),color = 'blue',label='20 degrees, Mean = %5.1f keV' %Popt_20[1],linewidth=2)
    plt.plot(Energy_Fit_30,gaussian(Energy_Fit_30,*Popt_30),color = 'orange',label='30 degrees, Mean = %5.1f keV' %Popt_30[1],linewidth=2)
    plt.plot(Energy_Fit_45,gaussian(Energy_Fit_45,*Popt_45),color = 'green',label='45 degrees, Mean = %5.1f keV' %Popt_45[1],linewidth=2)
    
    plt.title('Cs137 Gamma Ray Count vs Energy', weight = 'bold')
    plt.xlabel('Energy / keV', weight = 'semibold')
    plt.ylabel('Gamma Ray Count', weight = 'semibold')
    plt.legend(title = 'Scattering Angle',title_fontsize=15)
    plt.show()


#%%
'This last section looks at the relationship between scattering angle and mean energy'

Mean_Energy = [Popt_20[1],Popt_30[1],Popt_45[1]]
Angles = [20,30,45]

def PlotAngles1():
    plt.plot(Angles,Mean_Energy,'o',color='red',label='Experimental Data')
    plt.errorbar(Angles,Mean_Energy,yerr=errors,fmt='o',color='red',capsize = 5)
    plt.title('Gamma Ray Energy vs Scattering Angle', weight='bold')
    plt.legend()
    plt.xlim(0,Upper_Angle)
    plt.ylim(0,Upper_Energy)
    plt.xlabel('Scattering Angle/Degrees')
    plt.ylabel('Energy of Deflected Photons/keV')

def PlotAngles2():
    plt.plot(Angles[0],Mean_Energy[0],'o',color='blue',label='20°, Energy = %5.1f keV' %Mean_Energy[0])
    plt.plot(Angles[1],Mean_Energy[1],'o',color='orange',label='30°, Energy = %5.1f keV' %Mean_Energy[1])
    plt.plot(Angles[2],Mean_Energy[2],'o',color='green',label='45°, Energy = %5.1f keV' %Mean_Energy[2])
    plt.errorbar(Angles,Mean_Energy,yerr=errors,fmt='.',color='gray',capsize = 5)
    plt.title('Gamma Ray Energy vs Scattering Angle', weight='bold')
    plt.legend()
    plt.xlim(0,Upper_Angle)
    plt.ylim(0,Upper_Energy)
    plt.xlabel('Scattering Angle/Degrees')
    plt.ylabel('Energy of Deflected Photons/keV')

'Define Compton Scattering Formula'
def Scatter(E0,θ):
    return E0/(1+(E0/511)*(1-np.cos(θ*np.pi/180)))

def PlotCompton():
    xvals = np.arange(0,Upper_Angle,step = 0.01)
    yvals = Scatter(E0,xvals)
    plt.plot(xvals,yvals,color='gray',label = 'Relation between Photon Energy and Scattering Angle ')
    plt.title('Gamma Ray Energy vs Scattering Angle', weight='bold')
    plt.legend()
    plt.xlim(0,Upper_Angle)
    plt.ylim(0,Upper_Energy)















