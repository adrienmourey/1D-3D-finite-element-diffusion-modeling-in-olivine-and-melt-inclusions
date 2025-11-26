# 1D-3D-finite-element-diffusion-modeling-in-olivine-and-melt-inclusions
This repository contain all the python scripts used in the manuscript "Refining magmatic timescales and volatile transfer with 3D diffusion modeling of complex crystals" 

Adrien J. Mourey, adrien.mourey@ntu.edu.sg
Euan J.F. Mutch

CONTENTS
- Folder 1 "1D_3D_FEM_comparison_k1820_ol5B": Script modeling 1D diffusion of Fe–Mg between the core and rim of a melt inclusion hosted in an olivine crystal, to compare with a precomputed 3D diffusion dataset. The code uses the finite element method (FEM) through the FEniCS package (dolfin) to solve the diffusion equation along a linear transect, then compares the result to imported 3D model outputs at equivalent times.

- Folder 2 "1D_FEM_Fe_Mg_diffusion": Script modeling the 1D diffusion of Fe–Mg in olivine using DFENS or Dohmen & Chakraborty (2007) diffusion coefficients

- Folder 3 "3D_FEM_Fe_Mg_diffusion": Script modeling the 3D diffusion of Fe-Mg in melt inclusions and olivine using a pre-imported mesh. It applies the finite element method (FEM) with FEniCS (dolfin) to solve the diffusion partial differential equation, taking into account anisotropic, temperature-, pressure-, and oxygen fugacity-dependent diffusion coefficients.

- Folder 4 "3D_FEM_MgO_cooling_rate": Script modeling the 3D diffusion of Mg from olivine-hosted melt inclusions using an imported mesh and the finite element method with FEniCS. 

- Folder 5 "Analytical_solution_decompression_rate": Script using the analytical solution from Mutch et al. (2024) that model the diffusive equilibration of a spherical melt inclusion in an anisotropic host and Bayesian inference via Markov chain Monte Carlo (MCMC) to estimate the decompression rates for individual melt inclusions.

- Folder 6 "H2O_3D_FEM_diffusion_models": Script modeling the 3D diffusion of H2O in melt inclusions and olivine using a pre-imported mesh.
