# Script to model 3D diffusion of Mg from olivine-hosted melt inclusions using imported mesh

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

#############################################################################
# PARAMETERS FROM TERMINAL

f_mesh = sys.argv[1] # Mesh file

# set parameters here
Ts = float(sys.argv[2]) # Temperature (C) at the start
Te = float(sys.argv[3]) # Temperature at the end
f_fo2 = float(sys.argv[4]) # fO2 
f_P = float(sys.argv[5]) # Pressure in Pa
fo_average = float(sys.argv[6])
f_Kd = float(sys.argv[7])  # Olivine-melt partition coefficient for water
f_cooling = float(sys.argv[8]) # cooling rate - C/s?
f_nts = int(float(sys.argv[9])) # Number of time steps

Ts_K = Ts + 273.0 # Convert temperature into Kelvin
Te_K = Te + 273.0
fO2_Pa = f_fo2
P_Pa = f_P
#fo_average = (fo_rim + fo_core)/2.0
#ol_core = f_core
#ol_rim = f_rim


dT = (Ts_K - Te_K)/f_nts # change in pressure at each step
dTdt = f_cooling #MPa/s 0.01 0.16
dt1 = dT/dTdt
dt_c = Constant(dt1)

###########
# Need to consider relationship between Fe/Mg and T - geothermometer or empirical?

#def temp_mg(T_K, Kd):
#    T_C = T_K - 273.0
#    femg_melt = (-3.6926313343e-8)*(T_C**3.0) + (0.000141343755)*(T_C**2.0) + (-0.1807396823276)*(T_C) + 77.33299628856
#    femg_ol = femg_melt/Kd
#    return femg_ol
    
def temp_mg(T_K, Kd):
    # Rearranged Shea et al., 2022 thermometer - MgO in glass
    T_C = T_K -273.0
    mgo_melt = (T_C - 1017.0)/21.2
    mgo_ol = mgo_melt*Kd
    return mgo_ol

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


#xdmf_fo =  XDMFFile("./{0}_output_Fo.xdmf".format(f_mesh))

xdmf_femg =  XDMFFile("./{0}_{1}_output_FeMg_cooling.xdmf".format(f_mesh, f_cooling))


# Define finite element function space
Q = FunctionSpace(mesh, "CG", 1)

###############################################################################
# TIMESTEP AND DEGASSING CALCULATIONS

maxstep = f_nts #100 # maximum number of time steps


################################################################################
# DIFFUSION COEFFICIENTS

################################################################################
# Mutch version of the olivine diffusion coefficent
#def D0_mutch(C, T_i, P, lnfO2, lnD0, clnfO2, cXFo, cT_i, cP, cT_iP, aniso):
#    # Calculate diffusion coefficient given parameters generated in each realisation
#    lnD = lnD0 + clnfO2*lnfO2 + cXFo*C + cT_i*T_i + cP*P + cT_iP*T_i*P
    
#    D_001 = exp(lnD)*Constant(1e12)
#    D_100 = (Constant(1.0)/Constant(aniso))*D_001
#    D_010 = (Constant(1.0)/Constant(aniso))*D_001

#    D0 = as_matrix(((Constant(D_100), Constant(0.0), Constant(0.0)), (Constant(0.0), Constant(D_010), Constant(0.0)), (Constant(0.0), Constant(0.0), Constant(D_001))))
#    return D0
    
    
def D0(C, T_K, P_Pa, fO2_Pa):
    x1 = 10**(-9.21)
    x2 = ((fO2_Pa/(1e-7))**(1.0/6.0))
    x3 = (10.0**(3.0*(0.9-C)))
    x4 = np.exp(-(201000.0 + ((7e-6)*(P_Pa - 1e5)))/(8.314*T_K))
    
    D_001 = (x1*x2*x3*x4*(1e12))/6.0
    D_100 = D_001
    D_010 = D_001
    
    D0 = as_matrix(((Constant(D_100), Constant(0.0), Constant(0.0)), (Constant(0.0), Constant(D_010), Constant(0.0)), (Constant(0.0), Constant(0.0), Constant(D_001))))
    return D0 

def D1(T_K):
    D_m = exp(Constant(-7.895) - (Constant(26257.0/T_K)))*Constant((1e12))
    
    D1 = as_matrix(((Constant(D_m), Constant(0.0), Constant(0.0)), (Constant(0.0), Constant(D_m), Constant(0.0)), (Constant(0.0), Constant(0.0), Constant(D_m))))
    
    return D1

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

# Define function spaces for Fe-Mg and forsterite
#C0_fo = Function(Q)         
#C1_fo = Function(Q)
#S_fo = TestFunction(Q)

# Create additional functions

C0_femg = Function(Q)
C1_femg = TrialFunction(Q)
S_femg = TestFunction(Q)


CK0_femg = Function(Q)

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
#C_fo_mid = theta*C1_fo + (Constant(1.0)+Constant(-1.0)*theta)*C0_fo
C_femg_mid = theta*C1_femg + (Constant(1.0)+Constant(-1.0)*theta)*C0_femg

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

#C0_fo.vector()[:] = Constant(fo_core)


#a_fo, L_fo = lhs(F_fo), rhs(F_fo)
#J_fo = derivative(F_fo, C1_fo)

# boundary conditions to start with
#Cbc0_fo = DirichletBC(Q, fo_rim, fmarker, 2)  # fixed boundary conditions on exterior of mesh
#Cbcs_fo = [Cbc0_fo]



# Calculate composition in one domain
#CK0.vector()[:] = C0.vector()[:]

# Calculate the composition in the melt inclusion domain with the corresponding partition coefficient
#CK0.vector()[d1_dofs] = (C0.vector()[d1_dofs])/Kd

################################################################################
# TIMESTEPPING


#out_file = File("output_seguam/conc.pvd")

# Time-stepping
i = 0
t = 0
T_K = Ts_K

# Generate initial conditions

femg_i = temp_mg(T_K, Kd)

C0_femg.vector()[:] = Constant(femg_i)


while i < maxstep:

    # Combined diffusion coefficients for different parts of the mesh
    D = D1(T_K)*d_mark + D0(fo_average, T_K, P_Pa, fO2_Pa)*(1-d_mark)
    rho = ((dr_mol)/Kd)*d_mark + 1.0*(1-d_mark)
    
    #F_fo = rho*S_fo*(C1_fo-C0_fo)*dx + dt_c*(dot(rho*D*grad(C_fo_mid), grad(S_fo)))*dx
    F_femg = rho*S_femg*(C1_femg-C0_femg)*dx + dt_c*(dot(rho*D*grad(C_femg_mid), grad(S_femg)))*dx
    
    a_femg, L_femg = lhs(F_femg), rhs(F_femg)
    #J_femg = derivative(F_femg, C1_femg)
    
    # set boundary conditions
    
    femg_b = temp_mg(T_K, Kd) # Function that relates temperature and Mg/Fe

    Cbc0_femg = DirichletBC(Q, femg_b, fmarker, 1)  # fixed boundary conditions on exterior of mesh    
    Cbc1_femg = DirichletBC(Q, femg_b, fmarker, 2)  # fixed boundary conditions on exterior of mesh
    Cbcs_femg = [Cbc0_femg, Cbc1_femg]

    CK0_femg.vector()[:] = C0_femg.vector()[:]
    CK0_femg.vector()[d1_dofsx] = (C0_femg.vector()[d1_dofsx])/Kd
    
    #xdmf_fo.write(C0_fo, i, XDMFFile.Encoding.HDF5)
    
    solve(a_femg==L_femg, C0_femg, bcs=Cbcs_femg, solver_parameters={'linear_solver': 'cg', 'preconditioner': 'hypre_amg'}) # Solve with CG and algebraic multigrid
    #solve(F_fo==0, C1_fo, bcs=Cbcs_fo, J=J_fo, solver_parameters={'newton_solver':{'linear_solver': 'cg', 'preconditioner': 'hypre_amg'}}) # Solve with CG and algebraic multigrid
    #out_file << C1
    #C0_fo.assign(C1_fo)

        
    if mesh.mpi_comm().rank == 0:
        print('step {0}/{1}'.format(i, maxstep))
        print(T_K)
        print(femg_b)
    T_K -= dT
    t += dt1
    i += 1

xdmf_femg.write(CK0_femg, i, XDMFFile.Encoding.HDF5)

print('Model Complete!!!')
 









