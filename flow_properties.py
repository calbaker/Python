# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 15:01:59 2010

@author: calbaker
"""
import scipy as sp 

# Function for calculating the Reynolds number for a hydraulic diameter with 
# known temperature (K), pressure (kPa)
def Re_D_air(T,P,velocity,D):
    k_B = 1.38e-23 # Boltzmann's constant (J/K)
    Nhat = 6.022e26 # Avogadro's # (molecules/kmol)
    Rhat = k_B*1e-3*Nhat # Universal gas constant (kJ/kmol*K)
    Mhat_air = 28.964 # molar mass of air from Bird, Stewart, Lightfoot Table E.1 (kg/kmol)
    R_air = Rhat / Mhat_air # gas constant for air (kJ/kg*K)
    d_air = 3.617e-10 # collision diameter of "air molecule" from Bird, Stewart, Lightfoot Table E.1 (m)
    m_air = Mhat_air / Nhat # molecular mass (kg/molecule)
    rho_air = P/(R_air*T)
    mu = 5/16. * sp.sqrt(sp.pi*m_air*k_B*T)/(sp.pi*d_air**2) # viscosity of
    # exhaust flow from Bird, Stewart, Lightfoot Eq. 1.4-14 (Pa*s)
    Re_D_air = rho_air*velocity*D/mu
    return Re_D_air
   
# Turbulent pipe flow friction factor
def f(Re_D): 
    return 0.078*Re_D**(-1/4) # Adrian Bejan, Convection Heat Transfer, 3rd ed. Equation 8.13

# Nusselt number for turbulent pipe flow
def Nu_D(Re_D,Pr): 
    return 0.023*Re_D**(4/5.)*Pr**(1/3.) # Adrian Bejan, Convection Heat Transfer, 3rd ed. Equation 8.30

# Function for providing the mass based gas constant for air
def R_air():
    k_B = 1.38e-23 # Boltzmann's constant (J/K)
    Nhat = 6.022e26 # Avogadro's # (molecules/kmol)
    Rhat = k_B*1e-3*Nhat # Universal gas constant (kJ/kmol*K)
    Mhat_air = 28.964 # molar mass of air from Bird, Stewart, Lightfoot Table E.1 (kg/kmol)
    R_air = Rhat / Mhat_air # gas constant for air (kJ/kg*K)

    return R_air

# Function for providing molecular mass of air 
def m_air():
    Nhat = 6.022e26 # Avogadro's # (molecules/kmol)
    Mhat_air = 28.964 # molar mass of air from Bird, Stewart, Lightfoot Table E.1 (kg/kmol)
    m_air = Mhat_air / Nhat # molecular mass (kg/molecule)

    return m_air
    
# Function for providing c_p (kJ/kg-K) of air
def c_p_air(T):
    # Moran and Shapiro, Table A-21 constants for calculating specific heat of air
    alpha = 3.653 
    beta = -1.337e-3 
    gamma = 3.294e-6
    delta = 1.913e-9
    epsilon = 0.2763e-12
    coeff = sp.array([epsilon,delta,gamma,beta,alpha])
    c_p_air = sp.polyval(coeff,T)*R_air()
    
    return c_p_air
