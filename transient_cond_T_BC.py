# Created by Chad Baker on 10 Nov 2011

import numpy as np
import matplotlib.pyplot as plt

x_step = 0.001
x = np.arange(0, 0.1, x_step) # x-coordinate in slab (m)

t_step = 0.01
t = np.arange(0., 500., t_step)

T_cold = 500.

T_hot = 500. + 100. * np.sin(np.linspace(0,2.*np.pi,t.size/10)) 

alpha = 8.e-6
Fo = alpha * t_step / x_step**2

t_crit = x_step**2 / (2. * alpha)
margin = (t_crit - t_step) / t_crit * 100.
if t_step < t_crit:
    print "time step is", margin, "percent lower than necessary."
else:
    print "time step too large by", margin, "percent"

T = np.zeros([x.size, t.size])

T[:,0] = np.linspace(T_hot[0], T_cold, x.size)
T[-1,:] = T_cold

# creating and populating the coefficient matrix
coeff_mat = np.zeros([T.shape[0], T.shape[0]])
coeff_mat[0,0] = 1.
coeff_mat[-1,-1] = 1.
for pop in range(coeff_mat.shape[0]-2):
    coeff_mat[pop+1, pop] = Fo
    coeff_mat[pop+1, pop+1] = 1. - 2. * Fo
    coeff_mat[pop+1, pop+2] = Fo

# solving
for i in range(1,t.size):
    if i%10 == 0:
        T[0,i-1] = T_hot[i/10]
    T[:,i] = np.dot(coeff_mat, T[:,i-1])

markers = ['-', '--', '-.', ':']

plt.close('all')

fig1 = plt.figure()
for i in range(t.size):
    if i%250 == 0:
        plt.plot(x*100., T[:,i],marker=markers[i%3])

plt.xlabel('X coordinate (cm)')
plt.ylabel('Temperature (K)')
plt.grid()

plt.show()
