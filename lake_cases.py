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

    # -------- CODE --------
    # Temperature
    T_epi1 = Twinter
    T_hypo1 = Twinter
    # solids concnetrations
    Css_epi1 = 100
    Css_hypo1 = 100

    h_epi1 = hmax
    
    VAR_LAKE={"T_epi":T_epi1,"T_hypo":T_hypo1,"Css_epi":Css_epi1,"Css_hypo":Css_hypo1,"h_epi":h_epi1}

    return VAR_LAKE
    
    
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
    # -------- CODE --------
    # ----Initial setup----
    h_epi2 = hmax  # Thermocline depth
    T2 = VAR_LAKE["T_epi"]  # Temperature, constant
    Css2 = VAR_LAKE["Cssepi"] # Initial solid concentration constant

    rho_w2 = rho_TCss(T_epi2, Css_epi2) # total density
    m_w2 = rho_w2*A0*hmax # water mass
    # ----Special process in fall----
    # Heat flux and temperature change
    Hsurf = surfheat(Tw=T2,Ta,Wsp,RH,P,C)
    Qnet2 = -(Hsurf + Hsw0)*A0
    dT2 = (Qnet2/(m_w2*Cpw))/

    # turbulent mixing
    
    # iceon (not sure)
    if T_epi < Twinter:
        iceon = tyear

    VAR_LAKE={"T_epi":T2,"T_hypo":T2,"Css_epi":Css2,"Css_hypo":Css2,"h_epi":h_epi2}
    return VAR_LAKE, iceon
    
    
#%% ###########################################################################
def case3(VAR_LAKE,Hsurf,Hsw0,tyear,iceon,Vs,ustar,hmax=8,g=9.81,Cpw=4200,Twinter=4,A0=7.8*10**6,Dth=1.4E-7):
    """
    spring-summer, 2 boxes.
     
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

