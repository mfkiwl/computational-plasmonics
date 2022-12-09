#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plasmonic Eigenvalue Problem with two corner PMLs.


This script plots: (i) unperturbed spectrum, known analytically
                   (ii) spectrum computed using corner complex scaling
                   at each corner.

This demonstrates the existence of complex resonances and embedded eigenvalues 
for an ellipse perturbed by two corners. The result produced by this script 
have been shown at Mathematics of Wave Phenomena, Feb. 2022.
"""

import numpy as np
import matplotlib.pyplot as plt
import PlasmonicEigenvalueProblem_utils_fenicsx as PEP_utils
from PlasmonicEigenvalueProblem_utils_fenicsx import DIR_MESH
import PlasmonicEigenvalueProblem_utils_analytical as PEP_utils_ana
from petsc4py import PETSc
from slepc4py import SLEPc
import fenicsx_utils
import gmsh_utils
import os

def build_mesh(geofile,geo_param,gmsh_param):
    """ Build mesh file from geometry file. """
    gmshfile = os.path.splitext(geofile)[0]+'.msh'
    gmsh_utils.generate_mesh_cli(geofile,gmshfile,2,parameters=geo_param,**gmsh_param)
    return gmshfile

    # Common
SLEPc_params = {'nev': 50,
             'target': 0.65,
             'shift': 0.65,
             'problem_type': SLEPc.EPS.ProblemType.GNHEP,
             'solver': SLEPc.EPS.Type.KRYLOVSCHUR,
              'tol': 1e-8,
             'max_it': 1000}
OptDB = PETSc.Options()
OptDB["st_ksp_type"] = "preonly"
OptDB["st_pc_type"] = "lu"
OptDB["st_pc_factor_mat_solver_type"] = "mumps"

    # Quasi-circular case, with smaller perturbation
a_m =2.5; b_m = 2.45; c = np.sqrt(a_m**2-b_m**2)
phi_c = 0.75*np.pi # same angle for both corners
(x_c,y_c,x_j,y_j,R) = PEP_utils_ana.get_C1corner_ellipse(a_m,b_m,phi_c)
a_d =np.abs(x_c) + 4*R; b_d = np.sqrt(a_d**2-c**2)
     # Even major / even minor
kappa_ex_ev_ev_fun = lambda N : PEP_utils_ana.ellipse_eigenvalues(a_m, b_m, a_d, b_d, N)[0][0::2]
     # Even major / odd minor
kappa_ex_ev_odd_fun = lambda N : PEP_utils_ana.ellipse_eigenvalues(a_m, b_m, a_d, b_d, N)[0][1::2]
    # Odd major / even minor    
kappa_ex_odd_ev_fun = lambda N : PEP_utils_ana.ellipse_eigenvalues(a_m, b_m, a_d, b_d, N)[1][0::2]
    # Odd major / odd minor
kappa_ex_odd_odd_fun = lambda N : PEP_utils_ana.ellipse_eigenvalues(a_m, b_m, a_d, b_d, N)[1][1::2]
    # Unperturbed eigenvalues
k_ex_ev_ev = kappa_ex_ev_ev_fun(10)
k_ex_ev_odd = kappa_ex_ev_odd_fun(10)
k_ex_odd_odd = kappa_ex_odd_odd_fun(10)
k_ex_odd_ev = kappa_ex_odd_ev_fun(10)
    # Critical interval
PEP = PEP_utils.PEP_with_PML_fenicsx([],[])
eta = np.linspace(1e-16,8,200); # in (0,inf)
PEP.set_deformed_critical_interval([phi_c])
k_Ic_even = PEP.k_Ic_even[0](eta,1.0)
k_Ic_odd = PEP.k_Ic_odd[0](eta,1.0)

# Computation with 1 corner

# == Quasi-circular boundaries, small corner perturbation
    # -- Case definition
name = r"Quasi-circular ellipse with corner on major axis ($\phi$=0.75$\pi$)"
    # Corner
phi = [phi_c]
corner_pos = ["left"]
(x_c,y_c,x_m,y_m,R) = PEP_utils_ana.get_C1corner_ellipse(a_m,b_m,phi[0])
    # Mesh: structured layer (outer/inner)
fac = 0.06
a_o = (1-fac)*a_m+fac*a_d; b_o = np.sqrt(a_o**2-c**2)
a_i = c*np.cosh(2*np.arccosh(a_m/c) - np.arccosh(a_o/c))
b_i = np.sqrt(a_i**2-c**2)
    # Corner region offset
x_offset = [-a_d-0.7]
y_offset = [np.pi-b_d]
    # Complex resonances computed with COMSOL
kappa_ex_even_fun = lambda N : []
kappa_ex_odd_fun = lambda N : []
    # -- Mesh tags
tag_names = {
    'omega_m': ['omega-m'],
    'omega_d': ['omega-d'],
    'gamma_d': ['gamma-d'],
    'omega_m_pml': ['corner-omega-m-pml'],
    'omega_d_pml': ['corner-omega-d-pml'],
    'Euler_bottom_bnd': ['corner-bnd-bot'],
    'Euler_right_bnd': ['corner-bnd-right'],
}
    # -- Geometrical file: no buffer layer, i.e. the PML is all of the Euler
    # region (R_PML=R)
geofile=os.path.join(DIR_MESH,"Ellipse-with-1-corner_Structured-v1.geo")
gmshfile=os.path.join(os.getcwd(),"Ellipse-with-1-corner_Structured-v1.msh")
gmsh_param={
    'save_and_exit': True,
    'binary': True,
    'order' : 2,
    'meshing' : -1,
    'recombination' : -1,
    'flexible_transfinite' : True    
}
geo_param={
    'a_i':a_i,'b_i':b_i,
    'a_m':a_m,'b_m':b_m,
    'a_o':a_o,'b_o':b_o,
    'a_d':a_d,'b_d':b_d,
    'phi1':phi[0],
    'R1_TR': 1e-25, # R_TR/R radius of truncated region
    'x_ofst_1':x_offset[0], 'y_ofst_1': y_offset[0],
    'Nmu': 1, # element in structured layer (1 recommended)
    'Nint': 100, # nodes along sign-changing interface
    'GeomProgint': 1.04, # geometric progression (larger than 1)
    'Nz': 200, # element in z-direction in corner region
    'CharLengthMin_adim': 1, # Characteristic mesh size (ratio)
    'CharLengthMax_adim': 1,
    'GenerateQuadMesh': 1
}

# Mesh generation, assembly, and solve
gmshfile = build_mesh(geofile,geo_param,gmsh_param)
dmesh = fenicsx_utils.DolfinxMesh.init_from_gmsh(gmshfile,2)
phys_tag_1D = gmsh_utils.getPhysicalNames(gmshfile,1)
phys_tag_2D = gmsh_utils.getPhysicalNames(gmshfile,2)
    # Periodic boundary conditions
pbc = PEP_utils.get_pbc_ellipseNcorners_fenicsx(a_m,b_m,
    phi,corner_pos,x_offset,y_offset,
    [phys_tag_1D[s] for s in tag_names['Euler_bottom_bnd']],
    [phys_tag_1D[s] for s in tag_names['Euler_right_bnd']])
    # Physical entities
phys_tags = dict()
phys_tags['omega-m'] = [k for tag in tag_names['omega_m'] for k in phys_tag_2D[tag]]
phys_tags['omega-d'] = [k for tag in tag_names['omega_d'] for k in phys_tag_2D[tag]]
phys_tags['gamma-d'] = [k for tag in tag_names['gamma_d'] for k in phys_tag_1D[tag]]
phys_tags['omega-m-pml'] = [k for tag in tag_names['omega_m_pml'] for k in phys_tag_2D[tag]]
phys_tags['omega-d-pml'] = [k for tag in tag_names['omega_d_pml'] for k in phys_tag_2D[tag]]
    # Define PEP
PEP = PEP_utils.PEP_with_PML_fenicsx(dmesh,phys_tags,pbc=pbc,case_name=f"{name}",mesh_name=f"{os.path.basename(geofile)}")
    # Assemble
V_fem = ("CG", gmsh_param['order'])
PEP.assemble(V_fem)
    # Solve PEP
alpha = np.exp(1j*np.pi/5); # PML angle
PEP.solve(SLEPc_params,alpha=alpha,tolReIm=1e-6)
PEP.export_eigenfunction_pyvista(f"{os.path.basename(geofile)}")
k_FEM_1corner = [-l for l in PEP.eigval]
N_FEM_1corner = PEP.eigvec_r[0].size
alpha_1corner = PEP.alpha
# Computation with 2 corners
name = r"Quasi-circular ellipse with corners on major and minor axes ($\phi_1$=$\phi_2$=0.75$\pi$)"
    # Ellipse m
c = np.sqrt(a_m**2-b_m**2)
    # Corner
phi = [phi_c,phi_c] # major axis, minor axis, angle in (0,pi)
corner_pos = ["left","top"]
(x_c,y_c,x_m,y_m,R) = PEP_utils_ana.get_C1corner_ellipse(a_m,b_m,phi[0],pos=corner_pos[0])
R2 = PEP_utils_ana.get_C1corner_ellipse(a_m,b_m,phi[1],pos=corner_pos[1])[4]
    # Mesh: structured layer (outer/inner)
fac = 0.06
a_o = (1-fac)*a_m+fac*a_d; b_o = np.sqrt(a_o**2-c**2)
a_i = c*np.cosh(2*np.arccosh(a_m/c) - np.arccosh(a_o/c))
b_i = np.sqrt(a_i**2-c**2)
    # Corner region offset
x_offset = [-a_d-0.5]; y_offset = [np.pi-b_d]
x_offset.append(x_offset[0] + np.log(R/R2))
y_offset.append(y_offset[0] + 2*np.pi + 0.1)
    # Complex resonances computed with COMSOL
kappa_ex_even_fun = lambda N : []
kappa_ex_odd_fun = lambda N : []
    # -- Mesh tags
tag_names = {
    'omega_m': ['omega-m'],
    'omega_d': ['omega-d'],
    'gamma_d': ['gamma-d'],
    'omega_m_pml': ['corner-1-omega-m-pml','corner-2-omega-m-pml'],
    'omega_d_pml': ['corner-1-omega-d-pml','corner-2-omega-d-pml'],
    'Euler_bottom_bnd': ['corner-1-bnd-bot','corner-2-bnd-bot'],
    'Euler_right_bnd': ['corner-1-bnd-right','corner-2-bnd-right'],
}

    # -- Geometrical file: R_PML equal to R
geofile=os.path.join(DIR_MESH,"Ellipse-with-2-corner_Structured.geo")
gmshfile=os.path.join(os.getcwd(),"Ellipse-with-2-corner_Structured.msh")
gmsh_param={
    'save_and_exit': True,
    'binary': True,
    'order' : 2,
    'meshing' : -1,
    'recombination' : -1,
    'flexible_transfinite' : False    
}
geo_param={
    'a_i':a_i,'b_i':b_i,
    'a_m':a_m,'b_m':b_m,
    'a_o':a_o,'b_o':b_o,
    'a_d':a_d,'b_d':b_d, 
    'phi1': phi[0],
    'R1_TR': 1e-25, # R_TR/R radius of truncated region
    'x_ofst_1':x_offset[0], 'y_ofst_1': y_offset[0],
    'phi2':phi[1],
    'R2_TR': 1e-25, # R_TR/R radius of truncated region
    'x_ofst_2':x_offset[1], 'y_ofst_2': y_offset[1],
    'Nmu': 1, # element in structured layer (1 recommended)
    'Nint': 100, # nodes along sign-changing interface
    'GeomProgint': 1.04, # geometric progression (larger than 1)
    'Nz': 200, # element in z-direction in corner region
    'CharLengthMin_adim': 1, # Characteristic mesh size (ratio)
    'CharLengthMax_adim': 1,
    'GenerateQuadMesh': 1
}    
gmshfile = build_mesh(geofile,geo_param,gmsh_param)
dmesh = fenicsx_utils.DolfinxMesh.init_from_gmsh(gmshfile,2)
phys_tag_1D = gmsh_utils.getPhysicalNames(gmshfile,1)
phys_tag_2D = gmsh_utils.getPhysicalNames(gmshfile,2)
    # Periodic boundary conditions
pbc = PEP_utils.get_pbc_ellipseNcorners_fenicsx(a_m,b_m,
    phi,corner_pos,x_offset,y_offset,
    [phys_tag_1D[s] for s in tag_names['Euler_bottom_bnd']],
    [phys_tag_1D[s] for s in tag_names['Euler_right_bnd']])
    # Physical entities
phys_tags = dict()
phys_tags['omega-m'] = [k for tag in tag_names['omega_m'] for k in phys_tag_2D[tag]]
phys_tags['omega-d'] = [k for tag in tag_names['omega_d'] for k in phys_tag_2D[tag]]
phys_tags['gamma-d'] = [k for tag in tag_names['gamma_d'] for k in phys_tag_1D[tag]]
phys_tags['omega-m-pml'] = [k for tag in tag_names['omega_m_pml'] for k in phys_tag_2D[tag]]
phys_tags['omega-d-pml'] = [k for tag in tag_names['omega_d_pml'] for k in phys_tag_2D[tag]]
    # Define PEP
PEP = PEP_utils.PEP_with_PML_fenicsx(dmesh,phys_tags,pbc=pbc,case_name=f"{name}",mesh_name=f"{os.path.basename(geofile)}")
PEP.set_deformed_critical_interval(phi)
    # Unperturbed eigenvalues
def compute_unperturbed_eigenvalues(N,mu_m,mu_d):
    nev = np.r_[0:N]; nodd = np.r_[1:N+1];
    kappa_ex_ev = -np.tanh(nev*(mu_d-mu_m))*np.tanh(nev*mu_m)
    kappa_ex_odd = -np.tanh(nodd*(mu_d-mu_m))/np.tanh(nodd*mu_m)
    return (kappa_ex_ev,kappa_ex_odd)
PEP.set_unperturbed_eigenvalues(
    lambda N: compute_unperturbed_eigenvalues(N,np.arccosh(a_m/c),np.arccosh(a_d/c))[0],
    lambda N: compute_unperturbed_eigenvalues(N,np.arccosh(a_m/c),np.arccosh(a_d/c))[1])
PEP.set_exact_eigenvalues(kappa_ex_even_fun,kappa_ex_odd_fun)

    # Assemble
V_fem = ("CG", gmsh_param['order'])
PEP.assemble(V_fem)
    # Solve PEP
alpha = np.exp(1j*np.pi/5); # PML angle
SLEPc_params = {'nev': 50,
             'target': 0.65,
             'shift': 0.65,
             'problem_type': SLEPc.EPS.ProblemType.GNHEP,
             'solver': SLEPc.EPS.Type.KRYLOVSCHUR,
              'tol': 1e-8,
             'max_it': 1000}
OptDB = PETSc.Options()
OptDB["st_ksp_type"] = "preonly"
OptDB["st_pc_type"] = "lu"
OptDB["st_pc_factor_mat_solver_type"] = "mumps"
PEP.solve(SLEPc_params,alpha=alpha,tolReIm=1e-6)
PEP.export_eigenfunction_pyvista(f"{os.path.basename(geofile)}")
k_FEM_2corner = [-l for l in PEP.eigval]
N_FEM_2corner = PEP.eigvec_r[0].size
alpha_2corner = PEP.alpha

# Plot
title = "Elliptical particle perturbed by two corners\n"+f"$(a_m,b_m)$=({a_m:2.2g},{b_m:2.2g}), $(a_d,b_d)$=({a_d:2.2g},{b_d:2.2g}), "+r"$\phi_{major}=\phi_{minor}$"+f"={np.rad2deg(phi_c):3g}°"
f = plt.figure()
ax=f.add_subplot(1,1,1)
ax.plot(k_ex_ev_ev,0*k_ex_ev_ev,marker="o",linestyle='none',label=r'$\kappa_n$ (even/even)',color='blue',fillstyle='none')
ax.plot(k_ex_ev_odd,0*k_ex_ev_odd,marker="o",linestyle='none',label='$\kappa_n$ (even/odd)',color='blue')
ax.plot(k_ex_odd_odd,0*k_ex_odd_odd,marker="s",linestyle='none',label='$\kappa_n$ (odd/odd)',color='red',fillstyle='none')
ax.plot(k_ex_odd_ev,0*k_ex_odd_ev,marker="s",linestyle='none',label='$\kappa_n$ (odd/even)',color='red')
ax.plot(np.real(k_Ic_even),0*np.real(k_Ic_even),label=r'$I_c$ (even/even)',color='b',linestyle='--')
ax.plot(np.real(k_Ic_odd),0*np.real(k_Ic_odd),label=r'$I_c$ (odd/odd)',color='r',linestyle='--')
ax.plot(np.real(k_FEM_1corner),np.imag(k_FEM_1corner),label=r'FEM corner on major'+f'\n{N_FEM_1corner:d} DoF\narg('+r'$\alpha$'+f')={np.rad2deg(np.angle(alpha_1corner)):3g}°',linestyle='none',marker='x',mfc='none',color='k')
ax.plot(np.real(k_FEM_2corner),np.imag(k_FEM_2corner),label=r'FEM 2 corners'+f'\n{N_FEM_2corner:d} DoF\narg('+r'$\alpha$'+f')={np.rad2deg(np.angle(alpha_2corner)):3g}°',linestyle='none',marker='x',mfc='none',color='r')
ax.set_title(title)
ax.set_xlim([-1,-0.58])
ax.set_ylim([-0.2,0.05])
ax.set_xlabel(r"$\Re(\kappa)$")
ax.set_ylabel(r"$\Im(\kappa)$")
ax.set_aspect('equal') # orthonormal axis
#ax.legend(loc='upper left', bbox_to_anchor=(0, -0.4),
#          ncol=3, fancybox=True, shadow=False)

ax.legend(loc='upper left', bbox_to_anchor=(1, 1),
          ncol=1, fancybox=True, shadow=False)
# Create animation in Inkscape
f.savefig("Spectrum-two-corners.svg")