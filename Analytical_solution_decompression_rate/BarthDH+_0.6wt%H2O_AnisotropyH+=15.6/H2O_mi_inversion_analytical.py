from dolfin import *
import matplotlib.pyplot as plt
import math as m
#import water_melt as wm #water diffusivity function
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import random

from scipy.special import elliprf
from scipy.integrate import cumtrapz


import os
import sys
import time
import datetime

# MC
import pymultinest
import json

# plotting
import matplotlib.pyplot as plt
#import seaborn as sb

# multi thread
#from mpi4py import MPI
#from petsc4py import PETSc


# custom
import pmc
#import KC_fO2 as kc # Calculates oxygen fugacity from melt compositions
#import vcalc_cy as vc

# SETUP MPI variables-------------------------------------------------
size = MPI.comm_world.Get_size()
rank = MPI.comm_world.Get_rank()
comm = MPI.comm_self


###############################################################################


def Aniso_analytical(mi_r, f_Kd, f_aniso, f_dpdt, f_degas, p_degas, T_K):     
    """
    Function for calculating analytical solution for anisotropic equilibration of a melt inclusion. Returns the composition of the melt inclusion following ascent upon degassing.

    ...

    Parameters
    ----------
    mi_r : float
        Melt inclusion radius (microns)
    f_Kd : float
        olivine-melt partition coefficient
    f_aniso : float
        Diffusive anisotropy in the olivine
    f_dpdt : float
        Decompression rate (MPa s-1)
    f_degas : array
        Water content along the degassing pathway (wt%)
    p_degas : array
        Pressure along the degassing pathway (MPa)

    Returns
    -------
    cmi : array
        Water contents in the centre of the melt inclusion (wt%)
    """
    
    c = f_degas
    p = p_degas

    t = (p[0] - p) / f_dpdt
    
    # DIFFUSION COEFFICIENTS
    
    # Fo 84/85 

    # Olivine diffusion coefficient - H+ with anisotropy - values from Barth et al. (2019) # Fo - 80
    D_100 = ((9.6e-6)*m.exp(-125000.0/(8.314*T_K)))*(1e12)
    D_010 = (1.0/f_aniso)*D_100 # Assume only 10x anisotropy for now. 40 
    D_001 = (1.0/f_aniso)*D_100
    
    ### Ferris 2018 Diffusion coefficients # 85 - Fo 90
    #D_100 = (10.0**(-5.4))*m.exp(-130000.0/(8.314*T_K))*(1e12)
    #D_010 = (10.0**(-6.9))*m.exp(-130000.0/(8.314*T_K))*(1e12)
    #D_001 = (10.0**(-6.6))*m.exp(-130000.0/(8.314*T_K))*(1e12)

    dr_olm = 1.2 # density ratio of olivine/melt


    D_eff = np.sqrt(D_100 * D_010 * D_001) / elliprf(D_100**-1.0, D_010**-1.0, D_001**-1.0)

    tau = (1.0 / (dr_olm * f_Kd)) * mi_r**2 / (3.0 * D_eff)
    
      
    cmi = np.exp(-t / tau) * (cumtrapz(c * np.exp(t / tau) / tau, x=t, initial=0.0) + c[0]) 
    
    return cmi



def mod_diff(cube, nparams):
    dpdt = cube[0] 
    T = cube[1]
    kd = cube[2]
    mi_r = cube[3]
        
    T_K = T + 273.0
    aniso = 15.6
    
    C_end = Aniso_analytical(mi_r, kd, aniso, dpdt, H2Onew, Pnew, T_K)
    
    # Test for numerical stability and equilibrium
    # If unstable it means it is likely equilibrated - set to final value
    if m.isnan(C_end[-1]) == True:
        gl_mod = H2Onew[-1]
    else:    
        gl_mod = C_end[-1]
            
    return gl_mod


########################################################################################################################################
    
# loglikelihood function
    
# ---------------- log likelihood function
def loglike(model, data, data_err):
    return np.sum( -0.5*((model - data)/data_err)**2 )
    

#------------------------------------------------------------------------
# MAIN
#------------------------------------------------------------------------    
if __name__ == "__main__":
    
    
    ## DATA-------------------------------------------------
    
    # import data from files
    f_dat = sys.argv[1] # Melt inclusion datafile
    f_dg = sys.argv[2] # Degassing file
    f_inclusion = sys.argv[3]
    f_out = sys.argv[4]

    # Import eruption parameters

    df_p = pd.read_csv(f_dg)
    
    
    # Import decompression pathway

    dfm = pd.read_csv(f_dg)
    
    h2o_degas = dfm['H2O_liq'].values
    P_degas = dfm['P_MPa'].values  
    
    # Use 14 order polynomial to fit the degassing curve 
    
    #p14 = np.poly1d(np.polyfit(P_degas, h2o_degas, 14, rcond=None, full=False, w=None, cov=False))

    #maxstep = 10000

    #Pnew = np.linspace(P_degas[0], P_degas[-1], maxstep)
    #H2Onew = p14(Pnew)
    
    Pnew = P_degas
    H2Onew = h2o_degas   

    #print(h2o_degas)
    #print(P_degas)


    # Import data set
    
    df_dat = pd.read_csv(f_dat)
    
    # Filter dataset based on melt inclusion name

    df_datf = df_dat[df_dat['Label'] == f_inclusion]

    # Create global arrays for observations and errors
    gl_obs = df_datf['H2O'].values[0]
    gl_err = df_datf['H2O_err'].values[0] #df_datf['weighted_err'].values[0]

    cov_s = np.empty([0])

    mi_r = df_datf['MI_radius'].values[0]
    mi_err = df_datf['MI_radius_err'].values[0]
    
    # check for output
    #dir = '/'.join(f_out.split('/')[:-1])
    #if not os.path.exists(dir):
    #    print("Making directory for output:", dir)
    #    os.makedirs(dir)
        
    # PARAMETER SETUP-------------------------------------------------
    # number of weight parameters (melt region sections = N(weights) + 1)
    
    parameters = ["dpdt", "T", "Kd", "MI_radius"]
    pti = np.empty([len(parameters)])
    pti[:] = np.nan
    #pti[4:10], pti[10:16], pti[16:] = 0, 1 , 2

        
   #==============================================================================
   # for i in Ds.index.values:
    #     parameters.append(i)
    #==============================================================================
    n_dim = len(parameters)
    n_params = n_dim
        
    # MCMC-------------------------------------------------
    # setup prior cube
    #type of distribution LU = ln uniform, U = uniform, MG= multivariate gaussian
    ptype = ["LU"] * n_dim
    ptype[1:] = ["G"]*3
    #ptype = ["U"] * n_dim  
    #ptype[4:] = ["G"]*1
        
    pcube = np.full((n_dim,2), np.nan)
    pcube[0,:] = [np.log10(0.0000001), np.log10(10.0)]  # dpdt - MPa/s
    pcube[1,:] = [1210.0, 13.0]  # Temperature C 
    pcube[2,:] = [0.0009, 0.0001]  # Kd - Towbin et al., (2023)    
    pcube[3,:] = [mi_r, mi_err]  # Kd - Towbin et al., (2023) 
    
    invMC = pmc.Pmc(n_dim, ptype, pcube, cov_s, pti, loglike, mod_diff,
                    gl_obs, gl_err,
                    evidence_tolerance = 0.5, sampling_efficiency = 0.8,
                    n_live_points = 400, n_params= n_params)
    
    json.dump(parameters, open(f_out + 'params.json', 'w')) # save parameter names
    
    invMC.run_mc(f_out)
    result = invMC.result(f_out)
       
    if rank == 1 :
            # PLOT-------------------------------------------------
    # Prevents having to generate 4 plots when weaving on 4 separate processors
    #==============================================================================
        bf_params = result.get_best_fit()['parameters']
        C_bf = mod_diff(bf_params, len(bf_params))
                
        C_mod = Aniso_analytical(bf_params[3], bf_params[2], 15.6, bf_params[0], H2Onew, Pnew, bf_params[1] + 273.0)
        
        pd.DataFrame({'P_MPa': Pnew, 'H2O_liq': H2Onew, 'H2O_MI': C_mod}).to_csv(f_out + 'mod_cv.csv', sep=',')
            
        # data fit
        #fig, axes = plt.subplots(3,1)
        
        plt.errorbar(Pnew[-1], gl_obs, yerr = gl_err,fmt='o', color='grey', label='data')
        #plt.plot(dist, inicon['Mgnum'].values, color = 'black')
        plt.plot(Pnew, H2Onew, label = 'Liq degassing curve')
        plt.plot(Pnew, C_mod, label = 'MI fit')
        plt.ylabel('H2O (wt%)')
        plt.xlabel('P (MPa)')
        plt.legend(loc = 0, fontsize = 'x-small')
        #plt.errorbar(dist, obs['XFo'].values, yerr = obs['Fo_std'].values,fmt='o', color='red', label='data')
        #plt.plot(dist, C_Fo_bf, color = 'blue')
        
        
        plt.savefig(f_out + 'fit.png')
        plt.savefig(f_out + 'fit.pdf')
        plt.close()
       
       
       
       
       
       
       
        


