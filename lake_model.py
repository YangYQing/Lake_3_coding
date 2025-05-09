"""
Main program that models the lake (thermal structure + solids concentration). 
TO FILL!
"""

"comments"

import os
import numpy as np
from lake_functions import *
from lake_cases import *
import datetime as dtt
import matplotlib.pyplot as plt


#%% ###########################################################################
# PARAMETERS
###############################################################################

#Time
nyears=2
tyear=np.tile(np.arange(1,366,1),[1,nyears]) # DOY
tyear=tyear[0];
t=np.arange(1,365*nyears+1,1) # Time in days from the first day

# Meteorological data:
meteodata_NASA =np.genfromtxt("Meteo_data_NASA.txt", usecols=tuple(np.arange(3,9,1)),skip_header=1,dtype=float)
meteodata_NASA=np.tile(meteodata_NASA,[nyears,1]) 
Ta=meteodata_NASA[:,4]  # [°C]
RH=meteodata_NASA[:,2]  # [#]
Wsp=meteodata_NASA[:,5]  # [m/s]
R=meteodata_NASA[:,1]*10**6/(24*3600)  # [W/m^2]
P=meteodata_NASA[:,3]*1000  # [Pa]

# Atmospheric variables:
rho_air=1.2  # [kg/m^3], air density
CD=0.001  # [-], drag coefficient

#Basin:
hmax=8; # [m], lake depth

# Heat fluxes:
Twinter=4
Adiff=0.066  # [-], albedo of diffuse shortwave radiation
Adir=0.1  # [-], albedo of direct shortwave radiation (approximate value for 50°N)
C=0.5  # Cloudiness
Fdir=(1-C)/((1-C)+0.5*C) 
Fdiff=0.5*C/((1-C)+0.5*C) 
Hsw0=-R*(Fdir*(1-Adir)+Fdiff*(1-Adiff))  # [W/m^2], shortwave radiation (<0) passing through the lake surface

# Gravity
g=9.81 # [m.s^(-2)]

# Water and particles properties:
Dp=1.6*10**(-6)  # [m]
rho_p=2650 # [kg/m^3]
nu=1E-6  # [m^2/s]
Vs=g*(rho_p-1000)/1000*Dp**2/(18*nu)  # [m/s], settling velocity
f_fft=10  # [%], solids content in the FFT
C_FFT=f_fft*10**4  # [g/m^3], solids concentration in the FFT


#%% ###########################################################################
# THERMOCLINE
###############################################################################

iceoff=(dtt.datetime(2015,4,26)-dtt.datetime(2015,1,1)).days  # iceoff DOY (based on observations)
iceon=max(t)   # initial value for the iceon DOY, it will be calculated by the simulation
h_epi=np.zeros(np.shape(t)); h_epi[0]=hmax 

#%% ###########################################################################
# INITIALIZATION

 # Temperature:
T_epi=np.zeros(np.shape(t)); T_epi[0]=Twinter   # [°C]
T_hypo=np.zeros(np.shape(t)); T_hypo[0]=Twinter   # [°C]

 # Solids concentration:
Css_epi=np.zeros(np.shape(t)); Css_epi[0]=100   # [mg/L]
Css_hypo=np.zeros(np.shape(t)); Css_hypo[0]=100   # [mg/L]


#%% ###########################################################################
# ITERATIONS
###############################################################################

for k in np.arange(len(t)-2): 

    VAR_LAKE={"T_epi":T_epi[k],"T_hypo":T_hypo[k],"Css_epi":Css_epi[k],"Css_hypo":Css_hypo[k],"h_epi":h_epi[k]}

    # -------- TO FILL --------

#%% ###########################################################################
# PLOTS
###############################################################################

# -------- TO FILL --------
