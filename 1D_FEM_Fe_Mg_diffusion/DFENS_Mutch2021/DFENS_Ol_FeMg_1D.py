# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 10:56:36 2017

@author: ejfm2
"""
#FENICS
# from dolfin import *

# May need to rewrite so that intensive parameters are required for classes?
# read command line arguments
import os
import sys
import time

# numerical functions
import numpy as np
import scipy.integrate as integrate
import math as m
from scipy import special as sp
from scipy.interpolate import interp1d
from dolfin import *


# data
import pandas as pd

# MC
import pymultinest
import json

# plotting
import matplotlib.pyplot as plt
#import seaborn as sb

# custom
import pmc
#import KC_fO2 as kc

set_log_active(False)

# SETUP MPI variables-------------------------------------------------
size = MPI.comm_world.Get_size()
rank = MPI.comm_world.Get_rank()
comm = MPI.comm_self

# FUNCTIONS : DIFFUSION-------------------------------------------------
# Diffusion Functions
       
#==============================================================================
# def D(T_i, lnfO2, lnD0, clnfO2, cT_i):
#     lnD = lnD0 + clnfO2*lnfO2 + cT_i*T_i
#     
#     D = exp(lnD)*Constant(1e12)
#     
#     return D 
#==============================================================================

def D(C, T_i, P, lnfO2, lnD0, clnfO2, cXFo, cT_i, cP, cT_iP, aniso):
    # Calculate diffusion coefficient given parameters generated in each realisation
    lnD = lnD0 + clnfO2*lnfO2 + cXFo*C + cT_i*T_i + cP*P + cT_iP*T_i*P
    
    D_001 = exp(lnD)
    D_100 = (Constant(1.0)/Constant(aniso))*D_001
    D_010 = (Constant(1.0)/Constant(aniso))*D_001

    D = (D_100*Constant(m.cos(psi)**2.0) + D_010*Constant(m.cos(phi)**2.0) + D_001*Constant(m.cos(gamma)**2.0))*Constant(1e12)
    return D        
               
        
def modc(model, fcoords, dist):
    mod_i = model[::-1]
    intp = interp1d(fcoords[:,0], mod_i,bounds_error=False, fill_value= "extrapolate")
    modx = intp(dist)
    return modx
        

def ICs_import(dist_init, IC, xs, rc, Q):
    i_int = interp1d(dist_init, IC ,bounds_error=False, fill_value=rc)
    ic = i_int(xs)
    Cinit = ic[dof_to_vertex_map(Q)] 
    return Cinit
        
# Model function

def mod_diff(cube, nparams):
    
    t1, T1, lnfO2, P, comp0, comp1, intf1, DFo = cube[0], cube[1], cube[2], cube[3], cube[4], cube[5], cube[6], cube[7:nparams]

    t1 *= (86400.0) # convert time into seconds

    dt1 = t1/500.0
    dT1 = Constant(dt1)

    #P *= 1.0e8
    T1 += 273.0

    T1_i = 1.0/T1

    P *= 1.0e8 # Convert pressure to Pa 
    #lnfO2 = kc.fO2calc_eq7(melt_comp, T1, P, fe3)  # fO2 in bars from Kress and Carmichael

    min_dist = min(dist)
    max_dist = max(dist)
    n = 999
    R = Constant(0.008314)

    #Plagioclase Diffusion Models
    mesh = IntervalMesh(comm, n, min_dist, max_dist)
    xs = np.linspace(min_dist, max_dist, n+1)
    
    Q = FunctionSpace(mesh, "CG", 1)
    # Define boundaries
    def inner_boundary(x):
        return near(x[0], max_dist)

    #def edge_boundary1(x):
        #return near(x[0], intf2)

    def edge_boundary(x):
        return near(x[0], min_dist)
        
    cell_dim = mesh.topology().dim()
    cmarker = MeshFunction("size_t", mesh, cell_dim, 0)
    fmarker = MeshFunction("size_t", mesh, cell_dim - 1, 0)

    # Create cell marker
    for c in cells(mesh):
        if (c.midpoint().x() < intf1):
            cmarker[c] = 0
        else:
            cmarker[c] = 1

    # Mark facets between the two regions
    for f in facets(mesh):
        cc = f.entities(cell_dim)
        if len(cc) == 2 and cmarker[cc[0]] != cmarker[cc[1]]:
            fmarker[f] = 1

    dofmap = Q.dofmap()
    d1_dofs = []
    d0_dofs = []

    for c in cells(mesh): # compute dofs in the domains
        if cmarker[c] == 0:
            d0_dofs.extend(dofmap.cell_dofs(c.index()))
        elif cmarker[c] == 1:
            d1_dofs.extend(dofmap.cell_dofs(c.index()))

    # unique
    d1_dofs = list(set(d1_dofs))
    d0_dofs = list(set(d0_dofs))

    # An expression which maps the MeshFunction cmarker to cells
    class CellMark(UserExpression):
        def eval_cell(self, values, x, cell):
            values[0] = cmarker[cell.index]
        def value_shape(self):
            return ()

    cm = CellMark()

    # Construct weak form
    C0_fo = Function(Q)       # Composition at current time step Forsterite
    C1_fo = Function(Q)  # Composition at next time step
    S_fo = TestFunction(Q)
    
    T1_i = Constant(T1_i)

    P = Constant(P)
    lnfO2 = Constant(lnfO2)

    theta = Constant(0.5)
    C_mid_fo = theta*C1_fo + Constant((1.0-theta))*C0_fo


    D1_fo = D(C_mid_fo, T1_i, P, lnfO2, Constant(DFo[0]), Constant(DFo[1]), Constant(DFo[2]), Constant(DFo[3]), Constant(DFo[4]), Constant(DFo[5]), 6.0)
    F1_fo = S_fo*(C1_fo-C0_fo)*dx + dT1*(inner(D1_fo*grad(S_fo), grad(C_mid_fo)))*dx


    # Initial conditions - create simple step pattern in each domain
    C0_fo.vector()[d0_dofs] = Constant(comp0)
    C0_fo.vector()[d1_dofs] = Constant(comp1)

    inicon_1 = C0_fo.vector().get_local()

    # Boundary conditions
    Cbc0 = DirichletBC(Q, Constant(comp0), edge_boundary) 
    Cbcs = [Cbc0]
    
    #Cinit_1 = ICs_import(dist_init, inicon['Mgnum'].values, xs, rim_comp, Q)
    #C0.vector()[:] = Cinit_1
        
    # Run first diffusion model for t1 and T1   

    # Timestepping
    i = 0
    
    while i < t1:
        solve(F1_fo==0, C1_fo, Cbcs)
        C0_fo.assign(C1_fo)

        i += dt1 # Need to decide on time increment
                
    rcoords = mesh.coordinates()
        
    gl_mod = modc(C0_fo.vector().get_local(), rcoords, dist)

    dat_t1 = C0_fo.vector().get_local()

    #plt.errorbar(dist, obs['XFo'].values, yerr = obs['Fo_stdev'].values ,fmt='o', color='black', label='data')
    #plt.plot(rcoords, dat_t1[::-1], color = 'blue')
    #plt.plot(rcoords, inicon_1[::-1], color = 'blue', linestyle = '--')
    #plt.ylabel('XFo')
    #plt.xlabel('Distance (um)')
    #plt.savefig('ol_test_model_{0}_{1}_{2}_{3}.png'.format(t1/(86400.0), t2/(86400.0), T1 -273, T2-273))
    #plt.close()

           
    return gl_mod

    
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
    
    
    f_dat = sys.argv[1]
    #f_melt = sys.argv[2]
    f_modpar = sys.argv[2]
    f_angles = sys.argv[3]
    #f_inicon = sys.argv[4]
    #f_mk = sys.argv[5]
    
    f_out = sys.argv[4]
    
    obs = pd.read_csv(f_dat)
    
    #mc = np.genfromtxt(f_melt,delimiter=',',dtype=None, names=True, encoding = None) #Melt composition
    obs_mpar = pd.read_csv(f_modpar) # Model parameters
    obs_angle = pd.read_csv(f_angles)

    #inicon = pd.read_csv(f_inicon)
    #df_mk = pd.read_csv(f_mk)
    
#    if len(sys.argv) == 8:
#        f_boundary = sys.argv[7]
#        Fixed = False
#    else:
#        Fixed = True
    
    # Import vcov files
    # Mg in plagioclase
    fo_all_vcov = pd.read_csv('./D_vcov/FeMg_all_001_vcov.csv') 
    fo_TAMED_vcov = pd.read_csv('./D_vcov/FeMg_TAMED_001_vcov.csv')
    
    
    # combine different elements into a single global array
    
    gl_obs = obs['XFo'].values  
    gl_err = obs['Fo_stdev'].values
    
    dist = obs['Distance'].values
    #dist_init = inicon['Distance'].values

    #melt_comp = mc['Composition']

    psi, phi, gamma = m.radians(obs_angle['angle100P']), m.radians(obs_angle['angle010P']), m.radians(obs_angle['angle001P'])
    
    #rim_comp, core_comp = df_mk['Mgnum_markers'][0], df_mk['Mgnum_markers'][1]
    
    #n, L = 299, max(dist)
    #n, L = 299, max(dist)
    
    # check for output
    #dir = '/'.join(f_out.split('/')[:-1])
    #if not os.path.exists(dir):
    #    print("Making directory for output:", dir)
    #    os.makedirs(dir)
        
    # PARAMETER SETUP-------------------------------------------------
    # number of weight parameters (melt region sections = N(weights) + 1)
    parameters = ["t", "T", "lnfO2", "P", "C0", "C1", "b1", "lnD0_Fo", "clnfO2_Fo", "cXFo_Fo", "cT_i_Fo", "cP_Fo", "cT_iP_Fo"]

    pti = np.empty([len(parameters)])
    pti[:] = np.nan
    pti[7:] = 0
   
    cov_fo = fo_TAMED_vcov.to_numpy()
    
    cov_s = np.array([cov_fo])
    
    #==============================================================================
    # for i in Ds.index.values:
    #     parameters.append(i)
    #==============================================================================
    n_dim = len(parameters)
    n_params = n_dim
    
    # MCMC-------------------------------------------------
    # setup prior cube
    ptype = ["MG"] * n_dim
    ptype[0:1] = ["LU"]*1   #type of distribution LU = ln uniform, U = uniform
    ptype[1:4] = ["G"]*3
    ptype[4:7] = ["U"]*3

    #ptype[23:24] = ["G"]
        
    pcube = np.full((n_dim,2), np.nan)
    pcube[0,:] = [np.log10(obs_mpar['time'][0]), np.log10(obs_mpar['time'][1])]     # time (days)
    pcube[1,:] = [obs_mpar['T'][0], obs_mpar['T'][1]] # temperature (C)
    pcube[2,:] = [obs_mpar['lnfO2_bar'][0], obs_mpar['lnfO2_bar'][1]]     #log fo2_bars
    pcube[3,:] = [obs_mpar['P_kbar'][0], obs_mpar['P_kbar'][1]]     #P kbar
    pcube[4,:] = [obs_mpar['C0'][0], obs_mpar['C0'][1]]     # Composition zone 0
    pcube[5,:] = [obs_mpar['C1'][0], obs_mpar['C1'][1]]     # Composition zone 1
    pcube[6,:] = [obs_mpar['step'][0], obs_mpar['step'][1]]     # Boundary position 1
    pcube[7:,0] = np.array([-6.755, 2.244e-1, -7.181, -2.674e4, -5.213e-10, -1.028e-7]) # DFo
    
    invMC = pmc.Pmc(n_dim, ptype, pcube, cov_s, pti, loglike, mod_diff,
                    gl_obs, gl_err,
                    evidence_tolerance = 0.5, sampling_efficiency = 0.8,
                    n_live_points = 400, n_params= n_params)
    
    json.dump(parameters, open(f_out + 'params.json', 'w')) # save parameter names
    
    invMC.run_mc(f_out)
    result = invMC.result(f_out)
        
        # fiddle around with this to gaussian with variance and covariance, manually change prior cube
        
    if rank == 1 :
            # PLOT-------------------------------------------------
    # Prevents having to generate 4 plots when weaving on 4 separate processors
    #==============================================================================
        bf_params = result.get_best_fit()['parameters']
        C_bf = mod_diff(bf_params, len(bf_params))
        pd.DataFrame({'XFo': C_bf}).to_csv(f_out + 'mod_cv.csv', sep=',')
            
        # data fit
        #fig, axes = plt.subplots(3,1)
        
        plt.errorbar(dist, obs['XFo'].values, yerr = obs['Fo_stdev'].values,fmt='o', color='red', label='data')
        #plt.plot(dist, inicon['Mgnum'].values, color = 'black')
        plt.plot(dist, C_bf, color = 'blue')
        plt.ylabel('XFo')
        plt.xlabel('Distance ($\mu$m)')
        #plt.errorbar(dist, obs['XFo'].values, yerr = obs['Fo_std'].values,fmt='o', color='red', label='data')
        #plt.plot(dist, C_Fo_bf, color = 'blue')
        
        
        plt.savefig(f_out + 'fit.png')
        plt.close()
        
        



