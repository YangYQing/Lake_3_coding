"""
Three cases used by the model: winter (case1), fall turnover (case 2), stratified season (case 3).
For each case, computation of the temperatures, solids concentration and thermocline
depth at the next time step. 

TO FILL!
"""

import numpy as np
import sys
from lake_functions import alphaT,rho_TCss, seasonal_therm
from statistics import mean

#%% ###########################################################################
def case1(VAR_LAKE,Twinter=4,hmax=8):
    """
    Winter case: ice-covered lake, no heat fluxes = constant temperatures 
    and concentrations. 
     
    INPUTS:
        VAR_LAKE: dictionary containing
            T_epi, T_hypo: epilimnion and hypolimnion temperatures [°C]
            Css_epi, Css_hypo: epilimnion and hypolimnion solids concentrations [mg/L]
            h_epi: thermocline depth [m]
        Twinter: temperature of maximum density [°C]
        hmax: lake depth [m]
        
    OUTPUTS:
        VAR_LAKE: dictionary containing the same variables than VAR_LAKE but at the next time step
    """

    # -------- TO FILL --------
    
    
#%% ###########################################################################
def case2(VAR_LAKE,Hsurf,Hsw0,tyear,iceon,Vs,C_FFT,hmax=8,g=9.81,Cpw=4200,Twinter=4,A0=7.8*10**6):
    """
    Fall turnover case: one box with heat fluxes.
     
    INPUTS:
        VAR_LAKE: dictionary containing
            T_epi, T_hypo: epilimnion and hypolimnion temperatures [°C]
            Css_epi, Css_hypo: epilimnion and hypolimnion solids concentrations [mg/L]
            h_epi: thermocline depth [m]
        Hsurf: surface heat flux [W/m^2]
        Hsw0: surface radiation flux [W/m^2]
        tyear: time as DOY
        iceon: iceon date as DOY
        Vs: settling velocity [m/s]
        C_FFT: solids concentration in the FFT [g/m^3]
        hmax: lake depth [m]
        g: gravitational acceleration [m/s^2]
        Cpw: heat capacity of water [J.kg^(-1).K^(-1)]
        Twinter: temperature of maximum density [°C]
        A0: lake surface area [m^2]
        
        
    OUTPUTS:
        VAR_LAKE: dictionary containing the same variables than VAR_LAKE but at the next time step
        iceon: iceon date as DOY (modified only if ice formation after initial iceon date)
    """
   # -------- TO FILL --------
    
    
#%% ###########################################################################
def case3(VAR_LAKE,Hsurf,Hsw0,tyear,iceon,Vs,ustar,hmax=8,g=9.81,Cpw=4200,Twinter=4,A0=7.8*10**6,Dth=1.4E-7):
    """
    Fall turnover case: one box with heat fluxes.
     
    INPUTS:
        VAR_LAKE: dictionary containing
            T_epi, T_hypo: epilimnion and hypolimnion temperatures [°C]
            Css_epi, Css_hypo: epilimnion and hypolimnion solids concentrations [mg/L]
            h_epi: thermocline depth [m]
        Hsurf: surface heat flux [W/m^2]
        Hsw0: surface radiation flux [W/m^2]
        tyear: time as DOY
        iceon: iceon date as DOY
        Vs: settling velocity [m/s]
        ustar: friction velocity [m/s]
        hmax: lake depth [m]
        g: gravitational acceleration [m/s^2]
        Cpw: heat capacity of water [J.kg^(-1).K^(-1)]
        Twinter: temperature of maximum density [°C]
        A0: lake surface area [m^2]
        Dth: thermal molecular diffusivity [m^2/s]
        
  
    OUTPUTS:
        VAR_LAKE: dictionary containing the same variables than VAR_LAKE but at the next time step
        iceon: iceon date as DOY (modified only if ice formation after initial iceon date)
    """
    # -------- TO FILL --------

