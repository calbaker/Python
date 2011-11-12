# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 15:01:59 2010

@author: calbaker
"""
import scipy as sp

class flow:
    # Function for calculating the Reynolds number for a hydraulic diameter
    def get_Re_D(self):
        Re_D = self.rho*self.velocity*self.D/self.mu
        return Re_D
    # Turbulent pipe flow friction factor
    def get_f(self): 
        return 0.078*self.Re_D**(-1/4) # Adrian Bejan, Convection Heat Transfer, 3rd ed. Equation 8.1
    # Nusselt number for turbulent pipe flow
    def get_Nu_D(self): 
        return 0.023*self.Re_D**(4/5.)*self.Pr**(1/3.) # Adrian Bejan, Convection Heat Transfer, 3rd ed. Equation 8.30

class ideal_gas():
    def __init__(self,**kwargs):
        if 'Mhat' in kwargs:
            self.Mhat = kwargs['Mhat']
        else:
            self.Mhat = 28.964 # molar mass of air from Bird, Stewart, Lightfoot Table E.1 (kg/kmol)
        if 'd' in kwargs:
            self.d = kwargs['d']
        else:
            self.d = 3.617e-10 # collision diameter of "air molecule" from Bird, Stewart, Lightfoot Table E.1 (m)
    # Constant attributes for all gases
    k_B = 1.38e-23 # Boltzmann's constant (J/K)
    Nhat = 6.022e26 # Avogadro's # (molecules/kmol)
    Rhat = k_B*1e-3*Nhat # Universal gas constant (kJ/kmol*K)
    # Calculated attributes
    R = Rhat / Mhat # gas constant (kJ/kg*K)
    m = Mhat / Nhat # molecular mass (kg/molecule)
    def get_rho(self):
        rho = self.P/(self.R*self.T)
        return rho
    def get_mu(self):
        mu = 5./16. * sp.sqrt(sp.pi*self.m*self.k_B*self.T)/(sp.pi*self.d**2) # viscosity of
        # exhaust flow from Bird, Stewart, Lightfoot Eq. 1.4-14 (Pa*s)
    # Function for providing c_p (kJ/kg-K) of air
    def c_p_air(self):
        # Moran and Shapiro, Table A-21 constants for calculating specific heat of air
        alpha = 3.653 
        beta = -1.337e-3 
        gamma = 3.294e-6
        delta = 1.913e-9
        epsilon = 0.2763e-12
        coeff = sp.array([epsilon,delta,gamma,beta,alpha])
        c_p_air = sp.polyval(coeff,self.T)*R_air()
        return c_p_air
