# Script to model 3D diffusion of Fe-Mg from olivine and olivine-hosted melt inclusions using imported mesh

from dolfin import *
import matplotlib.pyplot as plt
import math as m

import pandas as pd
#import vcalc_cy as vc
import numpy as np
from scipy.interpolate import interp1d
import random
import math as m

import os
import sys
import time
import datetime

#parameters["refinement_algorithm"] = "plaza_with_parent_facets" #refinement
set_log_active(False) #Option to supress outputs


class Model_calc(object):
    
    def __init__(self, rcoords, mod, err, obs):
        self.rcoords = rcoords
        self.mod = mod
        self.obs = obs
        self.err = err
        
    def modc(self):
        mod_i = self.mod[::-1]
        intp = interp1d(self.rcoords[:,0], mod_i,bounds_error=False, fill_value= "extrapolate") #mod_i[0])
        modx = intp(dist)
        return modx
    
    def unc(self):
        un = np.ones(len(self.obs))*self.err
        return un
        
    def chisqu(self):
        diff = self.modc()-self.obs
        diff_2 = (diff**2)/((self.unc())**2)
        tt = np.sum(diff_2)
        return tt
    
                 
def modc(model, fcoords, dist):
    # Interpolate model points onto observation distances
    mod_i = model[::-1]
    intp = interp1d(fcoords[:,0], mod_i,bounds_error=False, fill_value= mod_i[0])
    modx = intp(dist)
    return modx   
################################################################################
# DIFFUSION COEFFICIENTS
    
def D(C, T_K, P_Pa, fO2_Pa):
    x1 = 10**(-9.21)
    x2 = ((fO2_Pa/(1e-7))**(1.0/6.0))
    x3 = (10.0**(3.0*(0.9-C)))
    x4 = np.exp(-(201000.0 + ((7e-6)*(P_Pa - 1e5)))/(8.314*T_K))
    
    D_001 = (x1*x2*x3*x4*(1e12))
    D_100 = D_001/6.0
    D_010 = D_001/6.0
    
    D = (D_100*Constant(m.cos(psi)**2.0) + D_010*Constant(m.cos(phi)**2.0) + D_001*Constant(m.cos(gamma)**2.0))
    return D 


def mod_diff(time, t1, T, fO2_Pa, P_Pa, femg_rim, femg_core, maxstep):

    t1 *= 86400.0 # convert time into seconds
    dt1 = t1/maxstep # only 500 time steps
    
    print(dt1)
    
    dt_c = Constant(dt1)
    
    T_K = T + 273.0

    min_dist = min(dist) #Minimum observed distance along profile
    max_dist = max(dist) # Maximum observed distance along profile
    n = 999 # Inital number of mesh points

    mesh = IntervalMesh(n, min_dist, max_dist)
    xs = np.linspace(min_dist, max_dist, n+1)

    # Define finite element function space
    Q = FunctionSpace(mesh, "CG", 1)

    #################################################################################
    # BOUNDARIES
    # Define boundaries
    def inner_boundary(x):
        return near(x[0], max_dist)

    def edge_boundary(x):
        return near(x[0], min_dist)
        
        
    ####################################################################################
    # CREATE FUNCTION SPACE AND BOUNDARY CONDITIONS

    # Create additional functions

    C0 = Function(Q)
    C1 = TrialFunction(Q)
    S = TestFunction(Q)

    # Theta method: theta = 0.0 is forward Euler, theta = 0.5 is Crank-Nicholson, theta = 1.0 is backward Euler
    theta = Constant(0.5)
    C_mid = theta*C1 + (Constant(1.0)+Constant(-1.0)*theta)*C0
    
    D1 = D(fo_average, T_K, P_Pa, fO2_Pa)
    print(D1)

    F = S*(C1-C0)*dx + dt_c*(inner(D1*grad(S), grad(C_mid)))*dx

    a, L = lhs(F), rhs(F)

    # Initial conditions for olivine 

    C0.vector()[:] = Constant(femg_core)

    Cbc0 = DirichletBC(Q, Constant(femg_rim), edge_boundary)  # fixed boundary conditions on exterior of mesh
    Cbc1 = DirichletBC(Q, Constant(femg_rim), inner_boundary)  # fixed boundary conditions on exterior of mesh
    Cbcs = [Cbc0, Cbc1]
    
    ################################################################################
    # Set up misfit function
    t_s = np.empty([0])
    
    rcoords = mesh.coordinates()
    
       
    ################################################################################    
    # TIMESTEPPING

    i = 0
    t = 0
    while i < maxstep:
        if t/(86400.0) == time:
            mod = modc(C0.vector().get_local(), rcoords, dist)    
            
        t_s = np.append(t_s, t/(86400.0))
        # Calculate chi square value 
    
        solve(a==L, C0, bcs=Cbcs)
            
        t += dt1
        i += 1


    print('Model Complete!!!')

    return mod, t_s





#############################################################################
# PARAMETERS FROM TERMINAL

f_dat = sys.argv[1] # Mesh file
f_angles = sys.argv[2]

# set parameters here
T = float(sys.argv[3]) # Temperature (C)
f_fo2 = float(sys.argv[4]) # fO2 
f_P = float(sys.argv[5]) # Pressure Pa?
fo_core = float(sys.argv[6]) # Core composition
fo_rim = float(sys.argv[7]) # Rim composition
femg_core = float(sys.argv[8]) # Core composition
femg_rim = float(sys.argv[9]) # Rim composition
f_time = float(sys.argv[10]) # total time (days)
f_nts = int(float(sys.argv[11])) # Number of time steps

fo_average = (fo_rim + fo_core)/2.0

###############################################################################
#FILES
# Files to be imported and written

# Import data files
obs = pd.read_csv(f_dat) #This is 3D model data

# Import angles file
obs_angle = pd.read_csv(f_angles) # Angles

dat_1 = obs['FeMg_1yr'].values
dat_2 = obs['FeMg_2yr'].values 
dat_5 = obs['FeMg_5yr'].values
dat_10 = obs['FeMg_10yr'].values
dat_20 = obs['FeMg_20yr'].values

#err_femg = obs['FeMg_sd'].values

# Profile distances
dist = obs['Distance'].values

psi, phi, gamma = m.radians(obs_angle['angle100P']), m.radians(obs_angle['angle010P']), m.radians(obs_angle['angle001P'])

# Add plotted data


####################################################################################
# RUN MODEL

femg_mod_1, t_mod_1 = mod_diff(360.0, f_time, T, f_fo2, f_P, femg_rim, femg_core, f_nts)
femg_mod_2, t_mod_2 = mod_diff(740.0, f_time, T, f_fo2, f_P, femg_rim, femg_core, f_nts)
femg_mod_5, t_mod_5 = mod_diff(1820.0, f_time, T, f_fo2, f_P, femg_rim, femg_core, f_nts)
femg_mod_10, t_mod_10 = mod_diff(3660.0, f_time, T, f_fo2, f_P, femg_rim, femg_core, f_nts)
femg_mod_20, t_mod_20 = mod_diff(7300.0, f_time, T, f_fo2, f_P, femg_rim, femg_core, f_nts)

####################################################################################
# PLOTTING

# Plot up curves from 3D and 1D models

plt.plot(dist, dat_1[::-1], label = '3D model 1 yr', color = plt.cm.viridis(0.0), linewidth = 1.5, linestyle = '-')
plt.plot(dist, dat_2, label = '3D model 2 yr', color = plt.cm.viridis(0.25), linewidth = 1.5, linestyle = '-')
plt.plot(dist, dat_5, label = '3D model 5 yr', color = plt.cm.viridis(0.5), linewidth = 1.5, linestyle = '-')
plt.plot(dist, dat_10, label = '3D model 10 yr', color = plt.cm.viridis(0.75), linewidth = 1.5, linestyle = '-')
plt.plot(dist, dat_20, label = '3D model 20 yr', color = plt.cm.viridis(1.0), linewidth = 1.5, linestyle = '-')

plt.plot(dist, femg_mod_1, '--', label = '1D Model 1 yr', color = plt.cm.viridis(0.0), linewidth = 1.5)
plt.plot(dist, femg_mod_2, '--', label = '1D Model 2 yr', color = plt.cm.viridis(0.25), linewidth = 1.5) 
plt.plot(dist, femg_mod_5, '--', label = '1D Model 5 yr', color = plt.cm.viridis(0.5), linewidth = 1.5) 
plt.plot(dist, femg_mod_10, '--', label = '1D Model 10 yr', color = plt.cm.viridis(0.75), linewidth = 1.5) 
plt.plot(dist, femg_mod_20, '--', label = '1D Model 20 yr', color = plt.cm.viridis(1.0), linewidth = 1.5) 
 
plt.legend(fontsize = 'x-small', loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2)
plt.ylabel('Fe/Mg')
plt.xlabel(r'Distance ($\mu$m)')


plt.tight_layout()
plt.savefig("./model_plt.pdf")
plt.savefig("./model_plt.png")
plt.close()




####################################################################################
# DATA OUTPUT

# Output model fits 1D



#obs['Model_fit'] = femg_mod

#obs.to_csv('./model_output.csv', sep=',')

df_mf = pd.DataFrame({'Distance': dist, '1D_1yr': dat_1[::-1], '1D_2yr': dat_2, '1D_5yr': dat_5, '1D_5yr': dat_5, '1D_10yr': dat_10, '1D_20yr': dat_20})

df_mf.to_csv('./1D_model_outputs.csv', sep=',')


