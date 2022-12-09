#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plasmonic Eigenvalue Problem with PML. The corner region is meshed
in Euler coordinates (z,theta).

"""

import numpy as np
import matplotlib.pyplot as plt
import fenics as fe
import gmsh_utils
import gmsh_utils_fenics
import fenics_utils
import SLEPc_utils
from slepc4py import SLEPc
import PlasmonicEigenvalueProblem_utils_fenics as PEP_utils
from PlasmonicEigenvalueProblem_utils_fenics import DIR_MESH
import os
import itertools
import matplotlib.pyplot as plt
import meshio_utils
#%% Geometrical parameters
geofile=os.path.join(DIR_MESH,"Ellipse-with-1-corner_Structured-v3.geo")
gmshfile=os.path.join(os.getcwd(),"Ellipse-with-1-corner_Structured-v3.msh")
    # Case name
name = "Perturbed ellipse"
    # Ellipse m
a_m =2.5; b_m = 1
c = np.sqrt(a_m**2-b_m**2)
    # Ellipse d
a_d =3.5; b_d = np.sqrt(a_d**2-c**2)
    # Corner
phi = 0.75*np.pi
R_PML = 0.5  # R_PML/R radius of PML region
R_TR = 1e-5 # R_TR/R radius of truncated region
    # Mesh: structured layer (outer/inner)
a_o = 1.03*a_m; b_o = np.sqrt(a_o**2-c**2)
a_i = c*np.cosh(2*np.arccosh(a_m/c) - np.arccosh(a_o/c))
b_i = np.sqrt(a_i**2-c**2)
    # Mesh: 
Nint=80 # nodes along sign-changing interface
GeomProgint=1.2 # geometric progression
Nmu = 1 # element in structured layer (1 recommended)
Nz = 10 # element in z-direction in corner region
    # Corner region offset
x_offset=-1.40*a_d
y_offset=-0
case_param={'a_i':a_i,'b_i':b_i,'a_m':a_m,'b_m':b_m,'a_o':a_o,'b_o':b_o,
            'a_d':a_d,'b_d':b_d, 'phi':phi,'R_PML':R_PML,'R_TR':R_TR,
            'x_ofst':x_offset, 'y_ofst': y_offset,
            'Nmu':Nmu,'Nint':Nint,'GeomProgint':GeomProgint,'Nz':Nz}
    # Periodic boundary conditions
per = PEP_utils.get_pbc_ellipse1corner(a_m,b_m,phi,R_TR,x_offset,y_offset,tol=1e-10)
    # Unperturbed eigenvalues
n = np.r_[1:10]
n_ev = np.concatenate(( [0], n )) # add 0 for even eigenvalues
mu_m = np.arccosh(a_m/c)
mu_d = np.arccosh(a_d/c)
kappa_ex_ev = -np.tanh(n_ev*(mu_d-mu_m))*np.tanh(n_ev*mu_m)
kappa_ex_odd = -np.tanh(n*(mu_d-mu_m))/np.tanh(n*mu_m)
#%% Load and plot mesh
gmsh_utils.generate_mesh_cli(geofile,gmshfile,2,refinement=0,log=1,\
                             parameters=case_param,order=1,binary=True)

meshio_utils.print_diagnostic_gmsh(gmshfile)
xdmfiles = meshio_utils.convert_gmsh_to_XDMF(gmshfile,prune_z=True)    
dmesh = fenics_utils.DolfinMesh.init_from_xdmf(xdmfiles['triangle'],xdmfiles['line'])
gmsh_phys = meshio_utils.read_phys_gmsh(xdmfiles['gmsh_physical_entities'])
PhysName_1D_tag = gmsh_phys[0]; PhysName_2D_tag = gmsh_phys[1]
#fe.plot(dmesh.mesh,title=f"{gmshfile} \n {dmesh.mesh.num_vertices()} vertices")
#%% Assemble standard formulation
FE = fe.FiniteElement("P", dmesh.mesh.ufl_cell(), 1) # 'CR' (best) or 'P'
(A_m_l,A_d_l,idx_ur,idx_ui, V_r, V_i) = PEP_utils.assemble_PEP_PML(FE, dmesh,per,
                        PhysName_2D_tag['omega-m']+PhysName_2D_tag['corner-omega-m'],
                        PhysName_2D_tag['omega-d']+PhysName_2D_tag['corner-omega-d'],
                        PhysName_1D_tag['gamma-d'],
                        PhysName_2D_tag['corner-omega-m-pml'],
                        PhysName_2D_tag['corner-omega-d-pml'],
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
    # Exact spectrum
psi = lambda eta : np.sinh(eta*np.pi)/np.sinh(eta*(np.pi-phi));
k_even = lambda eta,alpha : (1-psi(eta*alpha))/(1+psi(eta*alpha));
k_odd = lambda eta,alpha : (1+psi(eta*alpha))/(1-psi(eta*alpha));
eta = np.linspace(1e-16,8,200); # in (0,inf)
f = plt.figure()
ax=f.add_subplot(1,1,1)
ax.plot(np.real(k_even(eta,alpha)),np.imag(k_even(eta,alpha)),label=r'$I_{c,\alpha}^{\rm{even}}$',color='b',linestyle='-',marker='None',mfc='none')
ax.plot(np.real(k_odd(eta,alpha)),np.imag(k_odd(eta,alpha)),label=r'$I_{c,\alpha}^{\rm{odd}}$',color='r',linestyle='-',marker='None',mfc='none')
ax.plot(np.real(kappa_ex_ev),np.imag(kappa_ex_ev),label=r'$\kappa_{\rm{unperturbed}}^{\rm{even}}$',color='b',linestyle='none',marker='o',mfc='none')
ax.plot(np.real(kappa_ex_odd),np.imag(kappa_ex_odd),label=r'$\kappa_{\rm{unperturbed}}^{\rm{odd}}$',color='r',linestyle='none',marker='o',mfc='none')
ax.plot(-np.real(eigval),-np.imag(eigval),label=r'FEM ($\mathbb{C}$ DoF'+f'={len(idx_ur):.2g})',color='k',linestyle='none',marker='x',mfc='none')
ax.set_xlabel(r"$\Re(\kappa)$")
ax.set_ylabel(r"$\Im(\kappa)$")
title = f"{name} |  $\phi$={np.rad2deg(phi):3.4g}° " +r"$\frac{R_{TR}}{R}$"+f"={R_TR:3.1e}"
title += f"\nComplex scaling: "+r"$arg(\alpha)=$"+f"{np.rad2deg(np.angle(alpha)):3.4g}°"
ax.set_title(r"Spectrum | "+title)
ax.set_xlim([-3,0])
ax.set_ylim([-1,1])
ax.legend(loc="upper right",ncol=3)
ax.set_aspect('equal') # orthonormal axis
#%% Plot eigenfunction (using a dolfin function)
i_plot = 3
eigval_r = np.real(eigval[i_plot]); eigval_i = np.imag(eigval[i_plot])
eigvec_r = eigvec[i_plot]; eigvec_r = eigvec_r/abs(eigvec_r[idx_ur]).max()
u_fe=fe.Function(V_r,fe.PETScVector(eigvec_r)) # Re(u)
#u_fe=fe.Function(V_i,fe.PETScVector(eigvec_r)) # Im(u)
f = plt.figure()
ax=f.add_subplot(1,1,1)
coll=fe.plot(u_fe,title=title+f"\n"+r"$\kappa=$"+f"{eigval_r+1j*eigval_i:.4g}",cmap='jet')
ax.set_xlabel("x")
ax.set_ylabel("y")
plt.colorbar(coll)
ax.set_aspect('equal') # orthonormal axis
