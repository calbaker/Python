"""Classes for characterizing flows and determining properties of ideal
gases"""

import numpy as np
from scipy.integrate import quad

class flow(object):
    """Class for dealing with things that are common to all flows."""
    def set_Re(self):
        """Function for calculating the Reynolds number for a hydraulic
        diameter""" 
        self.Re_D = self.velocity * self.D / self.nu

class ideal_gas(flow):
    """Class for modeling ideal gases.  Inherits properties of
    flow""" 

    def __init__(self,**kwargs):
        """ Sets a bunch of constants common to all ideal gases, and
        some which are specific to air."""

        # Unless otherwise noted, all gas properties are coming from
        # Bird, Stewart, Lightfoot Transport Phenomena Table E.1
        if 'species' in kwargs:
            self.species = kwargs['species']
        else:
            self.species = 'air'

        if self.species == 'air':
            self.Mhat = 28.964 
            # molar mass of air from Bird, Stewart, Lightfoot Table
            # E.1 (kg/kmol)  
            self.d = 3.617e-10
            # collision diameter of "air molecule" from Bird, Stewart,
            # Lightfoot Table E.1 (m)  
            self.Pr = 0.74 
            # Pr of air from Bird, Stewart, Lightfoot Table 9.3-1              
        
        elif self.species == 'propane' or 'C3H8':
            self.Mhat = 44.10
            self.d = 4.934e-10

        if 'Mhat' in kwargs:
            self.Mhat = kwargs['Mhat']
        if 'd' in kwargs:
            self.d = kwargs['d']
        if 'Pr' in kwargs:
            self.Pr = kwargs['Pr']
            
        # Constant attributes for all gases
        self.k_B = 1.38e-26
        # Boltzmann's constant (kJ/K)
        self.Nhat = 6.022e26
        # Avogadro's # (molecules/kmol)
        self.Rhat = self.k_B * self.Nhat 
        # Universal gas constant (kJ/kmol*K) 
        # Calculated attributes
        self.R = self.Rhat / self.Mhat # gas constant (kJ/kg*K)
        self.m = self.Mhat / self.Nhat # molecular mass (kg/molecule)

    def get_entropy(self,T):
        """Returns entropy with respect to 0 K at 1 bar."""
        def get_integrand(T):
            integrand = self.get_c_p_air(T) / T
            return integrand
        entropy = (quad(get_integrand, 0.5, T)[0])
        return entropy

    def get_enthalpy(self,T):
        """Returns enthalpy."""
        enthalpy = (quad(self.get_c_p_air, 0., T)[0])
        return enthalpy

    def get_rho(self, T, P):
        """Returns density (kg/m^3)
        arguments: P(kPa), T(K)
        returns rho(kg/m**3)"""  
        rho = P / (self.R * T) # density (kg/m**3)
        return rho
    
    def get_n(self, rho):
        """Returns number density (#/m^3)
        argument: rho(kg/m**3)
        returns n(#/m**3)"""  
        n = rho / self.m # number density (#/m^3)

    def get_mu(self,T):
        """Returns viscosity (Pa*s) of ideal gas attached from Bird,
        Stewart, Lightfoot Eq. 1.4-14.  This expression works ok for
        nonpolar gases, even ones with multiple molecules.""" 
        mu = (5. / 16. * (np.pi * m * k_B * T)**0.5 / (np.pi * d**2))   
        return mu

    def get_c_p_air(self,T):
        """c_p (kJ/kg-K) of air calculated using Moran and Shapiro,
        Table A-21 constants for polynomial for specific heat of air"""  
        self.polyrep = np.poly1d([0.2763e-12, 1.913e-9, 3.294e-6,
        -1.337e-3, 3.653]) 
        c_p_air = self.polyrep(T) * self.R 
        return c_p_air

    def set_rho(self):
        """Sets density as rho (kg/m^3) and number density as n
        (#/m^3)"""  
        self.rho = self.get_rho(self.T, self.P)
        self.n = self.get_n(self.rho)

    def set_Temp_dependents(self):
        """Sets viscosity (Pa*s) of general ideal gas and specific
        heat (kJ/kg*K) of air.  For other gases, use a different
        specific heat correlation."""  
        self.mu = self.get_mu(self.T)       
        self.c_p_air = self.get_c_p_air(self.T)
        # constant pressure specific heat of air (kJ/kg*K)  
        self.entropy = self.get_entropy(self.T)
        self.enthalpy = self.get_enthalpy(self.T)

    def set_TempPres_dependents(self):
        """Sets temp dependent properties and then sets properties
        that depend on both pressure and temperature. Properties
        include density (kg/m^3), kinematic viscosity (m^2/s), thermal
        diffusivity (m^2/s), and thermal conductivity (kW/m-K)"""
        self.set_Temp_dependents()
        self.set_rho()
        self.nu = self.mu / self.rho # kinematic viscosity (m^2/s)

    def set_thermal_props(self):
        """Sets temp and press dependents, then sets thermal
        diffusivity (m^2/s) if Pr is known, and then sets thermal 
        conductivity (kW/m*K).""" 
        self.set_TempPres_dependents()
        self.alpha = self.nu / self.Pr # thermal diffusivity (m^2/s)
        self.k_air = (self.alpha * self.rho * self.c_p_air) # thermal
            # conductivity(kW/m-K) of air
    
