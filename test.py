# -*- coding: utf-8 -*-
"""
Created on Thu Dec 16 11:16:23 2010

@author: calbaker
"""

# Moran and Shapiro, Table A-21 constants for calculating specific heat of air
T = 300.
alpha = 3.653 
beta = -1.337e-3 
gamma = 3.294e-6
delta = 1.913e-9
epsilon = 0.2763e-12
coeff = sp.array([epsilon,delta,gamma,beta,alpha])
c_p_air = sp.polyval(coeff,T)

print c_p_air

c_p_dumb = alpha + beta*T + gamma*T**2 + delta*T**3 + epsilon*T**4

print c_p_dumb