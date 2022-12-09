#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plasmonic Eigenvalue Problem (PEP) for a corner in Euler coordinates (z,theta)
with a PML region.

This script:
    - assembles the PEP by separating real and imaginary parts
    - solves the generalized eigenvalue problem
    - compares the computed spectrum to the exact essential spectrum
"""

import numpy as np
import matplotlib.pyplot as plt
import fenics as fe
import gmsh_utils
import gmsh_utils_fenics
import fenics_utils
import SLEPc_utils
from slepc4py import SLEPc
import itertools
import meshio_utils
import os
    # PEP-specific
import PlasmonicEigenvalueProblem_utils_fenics as PEP_utils
from PlasmonicEigenvalueProblem_utils_fenics import DIR_MESH
#%% Geometrical parameters
geofile=os.path.join(DIR_MESH,"Corner-Euler-PML.geo")
gmshfile="Corner-Euler-PML.msh"
R=1; R_PML = 1e-2; R_TR =1e-10; # Radiuses
phi = np.pi/2; # Camember angle (rad)
Nr=10; # No. of element on (ln(R_TR),ln(R))
Ntheta=10; # No. of elements in (-pi,pi)
case_param={'R':R,'R_PML':R_PML,'R_TR':R_TR,'phi':phi,
            'Nz':Nr,'Ntheta':Ntheta}
    # Exact spectrum
psi = lambda eta : np.sinh(eta*np.pi)/np.sinh(eta*(np.pi-phi));
k_even = lambda eta,alpha : (1-psi(eta*alpha))/(1+psi(eta*alpha));
k_odd = lambda eta,alpha : (1+psi(eta*alpha))/(1-psi(eta*alpha));
#%% Load and plot mesh
gmsh_utils.generate_mesh_cli(geofile,gmshfile,2,refinement=0,binary=True,
                             parameters=case_param)
meshio_utils.print_diagnostic_gmsh(gmshfile)
xdmfiles = meshio_utils.convert_gmsh_to_XDMF(gmshfile,prune_z=True)    
dmesh = fenics_utils.DolfinMesh.init_from_xdmf(xdmfiles['triangle'],xdmfiles['line'])
gmsh_phys = meshio_utils.read_phys_gmsh(xdmfiles['gmsh_physical_entities'])
PhysName_1D_tag = gmsh_phys[0]; PhysName_2D_tag = gmsh_phys[1]
fe.plot(dmesh.mesh,title=f"{gmshfile} \n {dmesh.mesh.num_vertices()} vertices")
#%% Assemble standard formulation
FE = fe.FiniteElement("P", dmesh.mesh.ufl_cell(), 1) # 'CR' (best) or 'P'
    # Periodic boundary condition on theta, bottom -> top
per = PEP_utils.get_pbc_corner(0)
#per = fenics_utils.PeriodicBoundary(); per.init() # no periodic boundary cond.
    # Assembly command when all physical entities are defined
(A_m_l,A_d_l,idx_ur,idx_ui, V_r,V_i) = PEP_utils.assemble_PEP_PML(FE,dmesh,per,
                        PhysName_2D_tag['omega-m'], 
                        PhysName_2D_tag['omega-d'], 
                        PhysName_1D_tag['corner-bnd-right'],
                        PhysName_2D_tag['omega-m-pml'], 
                        PhysName_2D_tag['omega-d-pml'],
                        diag_Am=1e2, diag_Ad=1e-2)
idx_ur=np.array(idx_ur,dtype='int32'); idx_ui=np.array(idx_ui,dtype='int32');
#%% Solve sparse eigenvalue problem: using SLEPc
alpha = np.exp(1j*np.pi/8); # PML angle
A_m = A_m_l[0] + np.real(alpha)*A_m_l[1] + np.imag(alpha)*A_m_l[2] + \
      np.real(1/alpha)*A_m_l[3] + np.imag(1/alpha)*A_m_l[4];
A_d = A_d_l[0] + np.real(alpha)*A_d_l[1] + np.imag(alpha)*A_d_l[2] + \
      np.real(1/alpha)*A_d_l[3] + np.imag(1/alpha)*A_d_l[4];
Am_petsc = A_m.mat(); Ad_petsc = A_d.mat()
EPS = SLEPc_utils.solve_GEP_shiftinvert(Am_petsc,Ad_petsc,
                          problem_type=SLEPc.EPS.ProblemType.GNHEP,
                          solver=SLEPc.EPS.Type.KRYLOVSCHUR,
                          nev=20,tol=1e-6,max_it=500,
                          target=0.6,shift=0.6)
(eigval,eigvec)=SLEPc_utils.EPS_get_spectrum_ReImFormulation(EPS,
                              idx_ur.astype(np.int32),idx_ui.astype(np.int32))
    # filter lambda=1
mask=[((np.abs(x)+1)>1e-4) for x in eigval]
eigval = list(itertools.compress(eigval, mask))
eigvec = list(itertools.compress(eigvec, mask))
#%% Plot spectrum
# Add number of nodes
eta = np.linspace(1e-16,8,200); # in (0,inf)
f = plt.figure()
ax=f.add_subplot(1,1,1)
ax.plot(np.real(k_even(eta,alpha)),np.imag(k_even(eta,alpha)),label=r'$I_{c,\alpha}^{\rm{even}}$',linestyle='-',marker='None',mfc='none')
ax.plot(np.real(k_odd(eta,alpha)),np.imag(k_odd(eta,alpha)),label=r'$I_{c,\alpha}^{\rm{odd}}$',linestyle='-',marker='None',mfc='none')
ax.plot(-np.real(eigval),-np.imag(eigval),label=r'FEM ($\mathbb{C}$ DoF'+f'={len(idx_ur):.2g})',linestyle='none',marker='x',mfc='none')
ax.set_xlabel(r"$\Re(\kappa)$")
ax.set_ylabel(r"$\Im(\kappa)$")
ax.set_title(r"Spectrum $arg(\alpha)=$"+f"{np.rad2deg(np.angle(alpha)):3.4g} deg")
ax.set_xlim([-3,0])
ax.set_ylim([-1,1])
ax.legend(loc="upper right",ncol=3)
ax.set_aspect('equal') # orthonormal axis
#%% Plot eigenfunction (using a dolfin function)
i_plot = 6
eigval_r = np.real(eigval[i_plot]); eigval_i = np.imag(eigval[i_plot])
eigvec_r = eigvec[i_plot]; eigvec_r = eigvec_r/abs(eigvec_r[idx_ur]).max()
u_fe=fe.Function(V_r,fe.PETScVector(eigvec_r)) # Re(u)
#u_fe=fe.Function(V_i,fe.PETScVector(eigvec_r)) # Im(u)
f = plt.figure()
ax=f.add_subplot(1,1,1)
coll=fe.plot(u_fe,title=f"Eigenvalue {eigval_r+1j*eigval_i:.4g}",cmap='jet')
ax.set_xlabel("x")
ax.set_ylabel("y")
plt.colorbar(coll)
