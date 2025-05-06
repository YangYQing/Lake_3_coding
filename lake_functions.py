"""
List of functions used by the model .

"""

import numpy as np
import sys

#%% ###########################################################################
def surfheat(Tw,Ta,Wsp,RH,P,C):
    """
    Compute the surface heat flux Hsurf [W/m^2] at a given time (Hsurf>0 for cooling).
     
    INPUTS:
        Tw: lake temperature [°C]
        Ta: air temperature [°C]
        Wsp: wind speed [m/s]
        RH: relative humidity [%]
        P: atmospheric pressure [Pa]
        C: couldiness [-]
        
    OUTPUTS:
        Hsurf: dsurface heat flux [W/m^2]
    """
    
    # 0. Parameters
    AL=0.03 # [-]
    Ew=0.972 # [-], emissivity of water
    sigma=5.67E-8 # [W/m^2/K^4], Boltzmann constant 
    Cpa=1005 # [J/kg/K], specific heat of air
    Lv=2470E3 # [J/kg]
    a=1.09
    P=P/100 # [hPa]
    def e_sat(T): # T: numerical array with temperature values
        return 611.2*np.exp(17.62*T/(243.12+T)) # [Pa] with T in [°C]  
    
    # 1. Compute the heat fluxes 
         
    ea=e_sat(Ta)*RH/100  # [Pa]
         
    # a. Sensible heat
    f=4.8+1.98*Wsp+0.28*(Tw-Ta)  # [W/m^2/mbar]
    Hsens=-1*(-Cpa*P/(0.622*Lv)*f*(Tw-Ta))  # [W/m^2], P in [hPa]
     
    # b. Latent heat
    Hlat=-1*(-f*(e_sat(Tw)-ea)/100)  # [W/m^2]
     
    # c. Longwave radiation
    Ea=a*(1+0.17*C**2)*1.24*(ea/100/(Ta+273.15))**(1/7)  # [-]
    HlwIN=-1*((1-AL)*sigma*Ea*(Ta+273.15)**4)  # [W/m^2]
    HlwOUT=(Ew*sigma*(Tw+273.15)**4) 
     
    # d. Total
    Hsurf=HlwOUT+HlwIN+Hsens+Hlat 
    return Hsurf
        
#%% ###########################################################################
def seasonal_therm(h,dt,B,N2):
    """
    Calculate the thermocline depth hnew [m] at the next time step (t+dt) from 
    the current depth (time t).
    
    INPUTS : 
        h: depth of the termocline at time t [m]
        dt: time step [s]
        B: buoyancy flux at time t [W/kg]
        N2: buoyancy frequency at time t [s^(-2)]
    
    
    OUTPUTS:
        hnew: depth of the termocline at time t+dt [m]
        
    """
    
    # Parameters
    A=0.2  # [-], entrainment coefficient
    
    if N2<=0:
        sys.exit("ERROR: The water column is not stratified")
    
    if B>0:  # Unstable conditions
        for t in np.arange(1,dt/3600+1,1): # Hourly time step
            h=h+3600*(1+2*A)*B/(N2*h) 
        hnew=h
        
    else:  # Stable conditions
       hnew=6
   
    return hnew

###############################################################################
def rho_TCss(T,Css,rho_p=2650):
    """
    Calculate the density from temperature and solids concentration.
    
    INPUTS : 
        T: water temperature [°C]
        Css: solids concentration [mg/L]
        rho_p: density of particles [kg/m^3]
    
    
    OUTPUTS:
        rho_w: water density [kg/m^3]
        
    """
    def rho_T(T):
        return 999.84298 + (65.4891*T-8.56272*T**2+0.059385*T**3)*10**(-3) # [kg/m^3]
    
    rho_w=rho_T(T)+Css*10**(-3)*(1-rho_T(T)/rho_p) # [kg/m^3]
    return rho_w

###############################################################################
def alphaT(T):
    """
    Calculate the thermal expansivity of water.
    
    INPUTS : 
        T: water temperature [°C]
        
    OUTPUTS:
        alphaval: thermal expansivity [K^-1]
        
    """
    alphaval=10**(-6)*(-65.4891+17.12544*T-0.178155*T**2)
    return alphaval







