# -*- coding: utf-8 -*-
"""
Created on Mon Dec  6 15:45:20 2010

@author: calbaker
"""

# Ensure variable is defined
def is_variable(x):
    try:
        x
    except NameError:
        x = None
    
    # Test whether variable is defined to be None
    if x is None:
        print "Variable did not exist"
    else:
        print "Variable was defined as "+str(x)
