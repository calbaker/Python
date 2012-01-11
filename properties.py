"""Classes for characterizing flows and determining properties of ideal
gases"""

import numpy as np
from scipy.integrate import quad
from scimath.units import * 

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
            self.Mhat = 28.964 * mass.kg / substance.kmol 
            # molar mass of air from Bird, Stewart, Lightfoot Table
            # E.1 (kg/kmol)  
            self.d = 3.617e-10 * length.m 
            # collision diameter of "air molecule" from Bird, Stewart,
            # Lightfoot Table E.1 (m)  
            self.Pr = 0.74 
            # Pr of air from Bird, Stewart, Lightfoot Table 9.3-1              
        
        elif self.species == 'propane' or 'C3H8':
            self.Mhat = 44.10 * mass.kg / substance.kmol
            self.d = 4.934e-10 * length.m

        if 'Mhat' in kwargs:
            self.Mhat = kwargs['Mhat']
        if 'd' in kwargs:
            self.d = kwargs['d']
        if 'Pr' in kwargs:
            self.Pr = kwargs['Pr']
            
        # Constant attributes for all gases
        self.k_B = 1.38e-23 * energy.J / temperature.K 
        # Boltzmann's constant (J/K)
        self.Nhat = 6.022e26 / substance.kmol
        # Avogadro's # (molecules/kmol)
        self.Rhat = ( self.k_B * 1.e-3 * self.Nhat * energy.kJ /
        (substance.kmol * temperature.K) ) 
        # Universal gas constant (kJ/kmol*K) 
        # Calculated attributes
        self.R = self.Rhat / self.Mhat # gas constant (kJ/kg*K)
        self.m = self.Mhat / self.Nhat # molecular mass (kg/molecule)

    def set_rho(self):
        """Sets density as rho (kg/m^3) and number density as n
        (#/m^3)"""  
        self.rho = self.P / (self.R * self.T) # density (kg/m**3)
        self.n = self.rho / self.m # number density (#/m^3)

    def get_entropy(self,T):
        """Returns entropy with respect to 0 K at 1 bar."""
        T = T.value
        print "T =",T
        print type(T)
        def get_integrand(T):
            print "T =",T
            print type(T)
            integrand = self.get_c_p_air(T).value / T
            return integrand
        entropy = ( quad(get_integrand, 0.5, T)[0] * energy.kJ /
        (mass.kg * temperature.K) )   
        return entropy

    def get_enthalpy(self,T):
        """Returns enthalpy."""
        def get_integrand(T):
            integrand = self.get_c_p_air(T)
            return integrand
        enthalpy = ( quad(get_integrand, 0., T.value)[0] * energy.kJ /
        mass.kg )  
        return enthalpy

    def get_c_p_air(self,T):
        self.polyrep = np.poly1d([0.2763e-12, 1.913e-9, 3.294e-6,
        -1.337e-3, 3.653]) 
        # c_p (kJ/kg-K) of air
        # Moran and Shapiro, Table A-21 constants for calculating
        # specific heat of air 
        c_p_air = self.polyrep(T) * self.R 

    def set_Temp_dependents(self):
        """Sets viscosity (Pa*s) of general ideal gas and specific
        heat (kJ/kg*K) of air.  For other gases, use a different
        specific heat correlation."""  
        self.mu = (5./16. * (np.pi*self.m*self.k_B*self.T)**0.5 /
        (np.pi*self.d**2)) # viscosity (Pa*s) of ideal gas from Bird,
            # Stewart, Lightfoot Eq. 1.4-14.  This
            # expression works ok for nonpolar gases,
            # even ones with multiple molecules.  
        self.c_p_air = self.get_c_p_air(self.T.value)
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
        self.nu = self.mu/self.rho # kinematic viscosity (m^2/s)

    def set_alpha(self):
        """Sets temp and press dependents, then sets thermal
        diffusivity (m^2/s) if Pr is known, and then sets thermal 
        conductivity (kW/m*K).""" 
        self.set_TempPres_dependents()
        self.alpha = self.nu/self.Pr # thermal diffusivity (m^2/s)
        self.k_air = (self.alpha * self.rho * self.c_p_air) # thermal
            # conductivity(kW/m-K) of air
        
