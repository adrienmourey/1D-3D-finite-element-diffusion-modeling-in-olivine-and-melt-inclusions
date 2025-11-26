# Script to model 3D diffusion of water from olivine hosted melt inclusions using imported mesh

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


def mod_diff(t1, T, fO2_Pa, P_Pa, femg_rim, femg_core, maxstep):

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
    Cbcs = [Cbc0]
    
    ################################################################################
    # Set up misfit function
    dat_mf = np.empty([0])
    t_s = np.empty([0])
    
    rcoords = mesh.coordinates()
    
       
    ################################################################################    
    # TIMESTEPPING

    i = 0
    t = 0
    while i < maxstep:
        print(i)
        femg_m = Model_calc(rcoords, C0.vector().get_local(), err_femg, dat_femg)
        dat_mf = np.append(dat_mf, femg_m.chisqu())
        t_s = np.append(t_s, t/(86400.0)) 
        # Calculate chi square value 
    
        solve(a==L, C0, bcs=Cbcs)
            
        t += dt1
        i += 1

    # Re-run model and estimate minimum misfit curve

    C0.vector()[:] = Constant(femg_core)

    i = 0
    t = 0
    while i < maxstep:
        if t/(86400.0) == t_s[np.argmin(dat_mf)]:
            mod = modc(C0.vector().get_local(), rcoords, dist)

        # Calculate chi square value 
    
        solve(a==L, C0, bcs=Cbcs)
            
        t += dt1
        i += 1

    print('Model Complete!!!')

    return mod, t_s, dat_mf





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
obs = pd.read_csv(f_dat)

# Import angles file
obs_angle = pd.read_csv(f_angles) # Angles

dat_femg = obs['Fe/Mg'].values  
err_femg = obs['FeMg_sd'].values

# Profile distances
dist = obs['Distance'].values

psi, phi, gamma = m.radians(obs_angle['angle100P']), m.radians(obs_angle['angle010P']), m.radians(obs_angle['angle001P'])


####################################################################################
# RUN MODEL

femg_mod, t_mod, mf_mod = mod_diff(f_time, T, f_fo2, f_P, femg_rim, femg_core, f_nts)


####################################################################################
# PLOTTING

fig, axes = plt.subplots(nrows=2,ncols=1)

axes[0].errorbar(dist, dat_femg, yerr = err_femg, fmt = 'o', color = 'steelblue', ecolor = '0.7', label='Data', zorder = -3)
axes[0].plot(dist, femg_mod, '-', color = 'firebrick', label = 'Model')
axes[0].set_ylabel('Fe/Mg')
axes[0].set_xlabel(r'Distance ($\mu$m)')
axes[0].legend(loc = 0, fontsize = 'x-small')

axes[1].plot(t_mod, mf_mod, '-', color = 'firebrick')
axes[1].text(min(t_mod) + 0.05*t_mod[np.argmin(mf_mod)], max(mf_mod) - 0.05*max(mf_mod), 't = {0:.2f} days'.format(t_mod[np.argmin(mf_mod)]))

axes[1].set_xlim(0, t_mod[np.argmin(mf_mod)]*2.0)

axes[1].set_ylabel('Chi-sq Misfit')
axes[1].set_xlabel('time (days)')

plt.tight_layout()
plt.savefig("./model_plt.pdf")
plt.savefig("./model_plt.png")
plt.close()


####################################################################################
# DATA OUTPUT

obs['Model_fit'] = femg_mod

obs.to_csv('./model_output.csv', sep=',')

df_mf = pd.DataFrame({'time (days)': t_mod, 'Chi-sq_misfit': mf_mod})

df_mf.to_csv('./model_misfit.csv', sep=',')


