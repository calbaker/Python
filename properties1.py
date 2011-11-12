# Classes for characterizing flows and determining properties of ideal
# gases 

import scipy as sp

class flow:
    # Function for calculating the Reynolds number for a hydraulic
    # diameter 
    def set_Re_dependents(self):
        self.Re_D = self.velocity*self.D/self.nu
        if self.Re_D < 2300:
            if self.geometry == 'parallel plates':
                self.f = 24*self.Re_D**-1 # friction factor for laminar flow
                # between parallel plates
            else:
                print "You don't have the right geometry"
        else:
            self.f = 0.078*self.Re_D**(-1./4.) # friction factor for
                                        # turbulent flow
                                        
    # Nusselt number for turbulent pipe flow
    def set_Nu_D(self): # check to see what geometry this is for
        self.set_Re_dependents()
        if self.Re_D < 2300:
            if self.geometry == 'parallel plates':
                self.Nu_D = 8.235 # Bejan Table 3.2, both sides have
                                  # const q''
        else:
            self.Nu_D = 0.023*self.Re_D**(4/5.)*self.Pr**(1/3.) # Adrian
                                        # Bejan, Convection Heat
                                        # Transfer, 3rd ed. Equation
                                        # 8.30  

class ideal_gas(flow):
    def __init__(self,**kwargs):
        if 'Mhat' in kwargs:
            self.Mhat = kwargs['Mhat']
        else:
            self.Mhat = 28.964 # molar mass of air from Bird, Stewart,
                               # Lightfoot Table E.1 (kg/kmol) 
        if 'd' in kwargs:
            self.d = kwargs['d']
        else:
            self.d = 3.617e-10 # collision diameter of "air molecule"
                               # from Bird, Stewart, Lightfoot Table
                               # E.1 (m) 
        if 'Pr' in kwargs:
            self.Pr = kwargs['Pr']
        else:
            self.Pr = 0.74 # Pr of air from Bird, Stewart, Lightfoot
                           # Table 9.3-1 # 
        # Constant attributes for all gases
        self.k_B = 1.38e-23 # Boltzmann's constant (J/K)
        self.Nhat = 6.022e26 # Avogadro's # (molecules/kmol)
        self.Rhat = self.k_B*1e-3*self.Nhat # Universal gas constant
                                        # (kJ/kmol*K) 
        # Calculated attributes
        self.R = self.Rhat / self.Mhat # gas constant (kJ/kg*K)
        self.m = self.Mhat / self.Nhat # molecular mass (kg/molecule)

    def set_Temp_dependents(self):
        self.mu = (5./16. * sp.sqrt(sp.pi*self.m*self.k_B*self.T) /
        (sp.pi*self.d**2)) # viscosity (Pa*s) of ideal gas from Bird,
                           # Stewart, Lightfoot Eq. 1.4-14 (might be
                           # for diatomic only)
        # c_p (kJ/kg-K) of air
        # Moran and Shapiro, Table A-21 constants for calculating
        # specific heat of air 
        coeff = sp.array([0.2763e-12,1.913e-9,3.294e-6,-1.337e-3
        ,3.653]) 
        self.c_p_air = sp.polyval(coeff,self.T)*self.R # constant
                                        # pressure specific heat of
                                        # air (kJ/kg*K) 

    def set_TempPres_dependents(self):
        self.set_Temp_dependents()
        self.rho = self.P/(self.R*self.T) # density (kg/m**3)
        self.nu = self.mu/self.rho # kinematic viscosity (m^2/s)
        self.alpha = self.nu/self.Pr # thermal diffusivity (m^2/s)
        self.k_air = (self.alpha * self.rho * self.c_p_air) # thermal
                                        # conductivity(kW/m-K) of air
