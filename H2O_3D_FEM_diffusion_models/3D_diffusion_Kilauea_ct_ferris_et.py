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
import water_melt as wm

import os
import sys
import time
import datetime

#parameters["refinement_algorithm"] = "plaza_with_parent_facets" #refinement
set_log_active(False) #Option to supress outputs

#############################################################################
# PARAMETERS FROM TERMINAL

f_mesh = sys.argv[1] # Mesh file
f_melt = sys.argv[2] # Mesh coordinates information
f_degas = sys.argv[3] # Degassing curve curve file

# set parameters here
f_aniso = float(sys.argv[4]) # Diffusive anisotropy in olivine
f_nts = int(sys.argv[5]) # Number of time steps
f_Kd = float(sys.argv[6])  # Olivine-melt partition coefficient for water
f_dpdt = float(sys.argv[7]) # Magma decompression rate MPa/s
#f_h2oi = float(sys.argv[8]) # Initial water content in the melt 
#f_h2oe = float(sys.argv[9]) # water content at the end of decompression path based on embayment wt%
T = float(sys.argv[8]) # Temperature (C)

T_K = T + 273.0 # Convert temperature into Kelvin


###############################################################################
#FILES
# Files to be imported and written

# Import mesh here
mesh = Mesh()
hdf = HDF5File(mesh.mpi_comm(), './{0}.h5'.format(f_mesh), "r")
hdf.read(mesh, "/mesh", False)
domains = MeshFunction('size_t', mesh, mesh.topology().dim())
hdf.read(domains, "/domains")
boundaries = MeshFunction('size_t', mesh, mesh.topology().dim() - 1)
hdf.read(boundaries, "/boundaries")

# Import dataframe of coordinates and concentrations

#df_coords = pd.read_csv('./{0}/{1}.csv'.format(f_mesh, f_coords))

#X = df_coords['MI_cen_x'].values # For data collection coordinates
#Y = df_coords['MI_cen_y'].values
#Z = df_coords['MI_cen_z'].values

#mi_c = np.array([X[0], Y[0], Z[0]]) # Define centre of melt inclusion based on coordinates
#mm = Point(mi_c[0], mi_c[1], mi_c[2])
#mm_cell, mm_distance = mesh.bounding_box_tree().compute_closest_entity(mm)

xdmf =  XDMFFile("./{0}_output_ferris_{1}_{2}.xdmf".format(f_mesh, f_dpdt, T))
xdmf2 =  XDMFFile("./{0}_output_ferris_final_{1}_{2}.xdmf".format(f_mesh, f_dpdt, T))

# Define finite element function space
Q = FunctionSpace(mesh, "CG", 1)

###############################################################################
# TIMESTEP AND DEGASSING CALCULATIONS

maxstep = f_nts #100 # maximum number of time steps

# Import degassing curve from file
degas = pd.read_csv(f_degas)
df_melt = pd.read_csv(f_melt)

vc_P = degas['P_MPa'].values # Pressure 
vc_h2o = degas['H2O'].values # Water content 

dsteps = len(vc_h2o)

dgci = interp1d(vc_h2o, vc_P)
Ps, Pe = vc_P[0], vc_P[-1]
dp = (Ps - Pe)/maxstep # change in pressure at each step

# Calculate timestep based on decompression rate and change in pressure

dpdt = f_dpdt #MPa/s 0.01 0.16
dt = dp/dpdt
dt_c = Constant(dt) 

# interpolate degassing curve onto new array of pressures if degassing curve != maxstep

if dsteps != maxstep:
    dgc = interp1d(vc_P, vc_h2o)
    Pnew = np.linspace(Ps, Pe, maxstep)
    H2Onew = dgc(Pnew)   # use interpolation function returned by `interp1d`

else:
    H2Onew = vc_h2o

################################################################################
# DIFFUSION COEFFICIENTS

# Olivine diffusion coefficient - H+ with anisotropy - values from Barth et al. (2019)
D_100 = (10.0**(-5.4))*m.exp(-130000.0/(8.314*T_K))*(1e12)
D_010 = (10.0**(-6.9))*m.exp(-130000.0/(8.314*T_K))*(1e12)
D_001 = (10.0**(-6.6))*m.exp(-130000.0/(8.314*T_K))*(1e12)

D0 = as_matrix(((Constant(D_100), Constant(0.0), Constant(0.0)), (Constant(0.0), Constant(D_010), Constant(0.0)), (Constant(0.0), Constant(0.0), Constant(D_001))))

# Calculate a constant diffusivity using Ni and Zhang equation

D_H2O =  wm.Dwtri(vc_h2o[0], df_melt, T, Ps/1000.0, "H2Ot") # T_C, P_GPa

D1 = as_matrix(((Constant(D_H2O), Constant(0.0), Constant(0.0)), (Constant(0.0), Constant(D_H2O), Constant(0.0)), (Constant(0.0), Constant(0.0), Constant(D_H2O))))

#######################################################################################
# MESH TAGGING

cell_dim = mesh.topology().dim()
cmarker = MeshFunction("size_t", mesh, cell_dim, 0)
fmarker = MeshFunction("size_t", mesh, cell_dim - 1, 0)

vmarker = MeshFunction("size_t", mesh, cell_dim - 2, 0)
tmarker = MeshFunction("size_t", mesh, cell_dim - 2, 0)

bmarker = MeshFunction("bool", mesh, cell_dim, False)

# Section here to map tags from gmsh to tags used below

MI_flag = 2
Ol_flag = 4

eb = 3
ib = 1

for c in cells(mesh):
    if domains[c] == MI_flag:
        cmarker[c] = 1
    else:
        cmarker[c] = 0

for f in facets(mesh):
    if boundaries[f] == ib:
        fmarker[f] = 1
    elif boundaries[f] == eb:
        fmarker[f] = 2
    else:
        fmarker[f] = 0


for c in cells(mesh):
    for f in facets(c):
        if fmarker[f] == 2:
            bmarker[c] = True
        elif fmarker[f] == 1:
            bmarker[c] = True

####################################################################################
# CREATE FUNCTION SPACE AND BOUNDARY CONDITIONS

dr_olm = 1.2 # density ratio of olivine/melt

dr_mol = 1.0/1.2 # density ratio of melt/olivine

# Define function spaces
C0 = Function(Q)         
C1 = TrialFunction(Q)
S = TestFunction(Q)

CK0 = Function(Q)

# An expression which maps the MeshFunction cmarker to cells
class CellMark(UserExpression):
    def eval_cell(self, values, x, cell):
        values[0] = cmarker[cell.index]
    def value_shape(self):
        return ()
cm = CellMark()

# Create a DG0 Function which has the same values as cmarker
DG0 = FunctionSpace(mesh, "DG", 0)
d_mark = interpolate(cm, DG0)

# Theta method:
# theta = 0.0 is forward Euler
# theta = 0.5 is Crank-Nicholson
# theta = 1.0 is backward Euler
theta = Constant(0.5)
C_mid = theta*C1 + (Constant(1.0)+Constant(-1.0)*theta)*C0

# Set the partition coefficient
Kd = f_Kd  # Olivine/melt partition coefficient from Barth et al., (2019) # 0.001

dofmap = Q.dofmap()
dof_first, dof_last = dofmap.ownership_range() 

d1_dofs = []
d0_dofs = []

for c in cells(mesh): # compute dofs in the domains
    if cmarker[c] == 1:
        d1_dofs.extend(dofmap.cell_dofs(c.index()))
    else:
        d0_dofs.extend(dofmap.cell_dofs(c.index()))

# unique
d1_dofs = list(set(d1_dofs))
d0_dofs = list(set(d0_dofs))

# Calculate dofs for melt inclusion

unowned = dofmap.local_to_global_unowned()

dofs = list(filter(lambda dof: dofmap.local_to_global_index(dof) not in unowned, range(dof_last-dof_first)))

d1_dofsx = list(filter(lambda dof: dof in d1_dofs, dofs))

# Initial conditions for olivine 

# Use degassing curve start point as initial condition

ol_h2o = (H2Onew*10000.0)*Kd
ol_init = ol_h2o[0]
C0.vector()[:] = Constant(ol_init)

# Combined diffusion coefficients for different parts of the mesh
D = D1*d_mark + D0*(1-d_mark)
rho = ((dr_mol)/Kd)*d_mark + 1.0*(1-d_mark)

F = rho*S*(C1-C0)*dx + dt_c*(dot(rho*D*grad(C_mid), grad(S)))*dx

a, L = lhs(F), rhs(F)
#J = derivative(F, C1)

# boundary conditions to start with
Cbc0 = DirichletBC(Q, ol_h2o[0], fmarker, 2)  # fixed boundary conditions on exterior of mesh
Cbcs = [Cbc0]

# Calculate composition in one domain
#CK0.vector()[:] = C0.vector()[:]

# Calculate the composition in the melt inclusion domain with the corresponding partition coefficient
#CK0.vector()[d1_dofs] = (C0.vector()[d1_dofs])/Kd

################################################################################
# TIMESTEPPING

#out_file = File("output_seguam/conc.pvd")


# Create file signifying time steps and pressures etc.
arr_j = np.empty([0])
arr_i = np.empty([0])
arr_P = np.empty([0])
arr_ol = np.empty([0])
arr_h2o = np.empty([0])
arr_t = np.empty([0])  

# Time-stepping
i = 0
t = 0
j = 0
while i < maxstep:

    CK0.vector()[:] = C0.vector()[:]
    CK0.vector()[d1_dofsx] = (C0.vector()[d1_dofsx])/Kd
    
    
    if i % 100 == 0:
        arr_j = np.append(arr_j, j)
        arr_i = np.append(arr_i, i)
        arr_t = np.append(arr_t, t)
        arr_P = np.append(arr_P, Pnew[i])
        arr_ol = np.append(arr_ol, ol_h2o[i])
        arr_h2o = np.append(arr_h2o, H2Onew[i])  
        xdmf.write(CK0, j, XDMFFile.Encoding.HDF5)
        j += 1

    # changing boundary conditions
    Cbc0 = DirichletBC(Q, Constant(ol_h2o[i]), fmarker, 2)  # fixed boundary conditions on exterior of mesh - changing with degassing curve.
    Cbcs = [Cbc0]
    
    solve(a==L, C0, bcs=Cbcs, solver_parameters={'linear_solver': 'cg', 'preconditioner': 'hypre_amg'})
        
    if mesh.mpi_comm().rank == 0:
        print('step {0}/{1}'.format(i, maxstep))

    t += dt
    i += 1


CK0.vector()[:] = C0.vector()[:]
CK0.vector()[d1_dofsx] = (C0.vector()[d1_dofsx])/Kd
    
xdmf2.write(CK0, i, XDMFFile.Encoding.HDF5)

print('Model Complete!!!')

# output compiled model files

df_out = pd.DataFrame({'Index': arr_j, 'Model step': arr_i, 'Time (s)': arr_t, 'Pressure (MPa)': arr_P, 'H2O_ol boundary (ppm)': arr_ol, 'H2O_melt (wt%)': arr_h2o})
df_out.to_csv("./{0}_output_ferris_{1}_{2}_stats.csv".format(f_mesh, f_dpdt, T), sep=',')



