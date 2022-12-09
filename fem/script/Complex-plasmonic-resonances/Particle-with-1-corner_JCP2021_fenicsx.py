#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plasmonic Eigenvalue Problem with one corner PML. 

The corner region is meshed in Euler coordinates (z,theta). This script
reproduces the results (cases A and B) of 10.1016/j.jcp.2021.110433.

"""
import os
    # Matrix
from petsc4py import PETSc
import numpy as np
    # Eigenvalue problem
from slepc4py import SLEPc
import SLEPc_utils
    # Dolfinx
import dolfinx
import fenicsx_utils
    # Mesh
import gmsh_utils
import meshio_utils
    # Plot
import matplotlib.pyplot as plt
import pyvista as pv
    # PEP-specific
import PlasmonicEigenvalueProblem_utils_fenicsx as PEP_utils
from PlasmonicEigenvalueProblem_utils_fenicsx import DIR_MESH
import PlasmonicEigenvalueProblem_utils_analytical as PEP_ana

def build_mesh(geofile,geo_param,gmsh_param):
    """ Build mesh file from geometry file. """
    gmshfile = os.path.splitext(geofile)[0]+'.msh'
    gmsh_utils.generate_mesh_cli(geofile,gmshfile,2,parameters=geo_param,**gmsh_param)
    return gmshfile

def build_PEP(gmshfile,case_name,
              a_m,b_m,a_d,b_d,
              phi,corner_pos, x_offset, y_offset,
              tag_names):
    """ Instanciate a PEP object from mesh file. """
        # Mesh importation
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
    PEP = PEP_utils.PEP_with_PML_fenicsx(dmesh,phys_tags,pbc=pbc,
                                         case_name=f"{name}",
                                         mesh_name=f"{os.path.basename(gmshfile)}")
    PEP.set_deformed_critical_interval(phi)
        # Unperturbed eigenvalues
    PEP.set_unperturbed_eigenvalues(
        lambda N : PEP_ana.ellipse_eigenvalues(a_m, b_m, a_d, b_d, N)[0],
        lambda N : PEP_ana.ellipse_eigenvalues(a_m, b_m, a_d, b_d, N)[1])
    return PEP

# Common
SLEPc_params = {
    'nev': 50,
    'target': 0.65,
    'shift': 0.65,
    'problem_type': SLEPc.EPS.ProblemType.GNHEP,
    'solver': SLEPc.EPS.Type.KRYLOVSCHUR,
    'tol': 1e-8,
    'max_it': 1000
}
OptDB = PETSc.Options()
OptDB["st_ksp_type"] = "preonly"
OptDB["st_pc_type"] = "lu"
OptDB["st_pc_factor_mat_solver_type"] = "mumps"
#%% Case A from JCP
name = "JCP Case A"
    # -- Case definition
(a_m,b_m,a_d,b_d,x_c,y_c,x_m,y_m,phi,corner_pos) = PEP_ana.perturbed_ellipse_parameters("JCP2021_case_A")
c = np.sqrt(a_m**2-b_m**2)
    # Mesh: structured layer (outer/inner)
fac = 0.15
a_o = (1-fac)*a_m+fac*a_d; b_o = np.sqrt(a_o**2-c**2)
a_i = c*np.cosh(2*np.arccosh(a_m/c) - np.arccosh(a_o/c))
b_i = np.sqrt(a_i**2-c**2)
    # Corner region offset
x_offset = [-1.40*a_d]
y_offset = [np.pi-b_d]
    # Complex resonances computed with COMSOL
kappa_ex_even_fun = lambda N : [-0.6547-1j*0.033,-0.8086-1j*0.025,-0.8847-1j*0.015,-0.9297-1j*0.00499]
kappa_ex_odd_fun = lambda N : []
    # -- Geometrical file: no buffer layer, i.e. the PML is all of the Euler
    # region (R_PML=R)
geofile=os.path.join(DIR_MESH,"Ellipse-with-1-corner_Structured-v1.geo")
gmshfile=os.path.join(os.getcwd(),"Ellipse-with-1-corner_Structured-v1.msh")
tag_names = {
    'omega_m': ['omega-m'],
    'omega_d': ['omega-d'],
    'gamma_d': ['gamma-d'],
    'omega_m_pml': ['corner-omega-m-pml'],
    'omega_d_pml': ['corner-omega-d-pml'],
    'Euler_bottom_bnd': ['corner-bnd-bot'],
    'Euler_right_bnd': ['corner-bnd-right'],
}
gmsh_param = {
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
    'R1_TR': 1e-25, # R_TR/R radius of truncated region,
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
PEP = build_PEP(gmshfile,name,
              a_m,b_m,a_d,b_d,
              phi,corner_pos, x_offset, y_offset,
              tag_names)
PEP.set_exact_eigenvalues(kappa_ex_even_fun,kappa_ex_odd_fun)
V_fem = ("CG", gmsh_param['order'])
PEP.assemble(V_fem)
alpha = np.exp(1j*np.pi/5); # PML angle
PEP.solve(SLEPc_params,alpha=alpha,tolReIm=1e-6)
eta = np.linspace(1e-16,8,200); # in (0,inf)
f = plt.figure()
ax = PEP.plot_spectrum(f,eta,plot_mesh_name=True)
ax.set_xlim([-1,-0.2])
ax.set_ylim([-0.17,0.1])
ax.xaxis.set_major_locator(plt.MultipleLocator(0.1))
ax.xaxis.set_major_formatter('{x:.1f}')
ax.grid(b=True, axis='x',which='major', color='gray', linestyle='-')
ax.xaxis.set_minor_locator(plt.MultipleLocator(0.01))
ax.minorticks_on()
ax.grid(b=True, which='minor', color='gray', linestyle='-', alpha=0.4)
PEP.export_eigenfunction_pyvista(f"{name}-{os.path.basename(geofile)}")

    # -- Geometrical file: with buffer layer, i.e. the PML is part of the Euler
    # region (R_PML<R).
geofile=os.path.join(DIR_MESH,"Ellipse-with-1-corner_Structured-v3.geo")
gmshfile=os.path.join(os.getcwd(),"Ellipse-with-1-corner_Structured-v3.msh")
tag_names = {
    'omega_m': ['omega-m','corner-omega-m'],
    'omega_d': ['omega-d','corner-omega-d'],
    'gamma_d': ['gamma-d'],
    'omega_m_pml': ['corner-omega-m-pml'],
    'omega_d_pml': ['corner-omega-d-pml'],
    'Euler_bottom_bnd': ['corner-bnd-bot'],
    'Euler_right_bnd': ['corner-bnd-right'],
}
gmsh_param = {
    'save_and_exit': True,
    'binary': True,
    'order' : 2,
    'meshing' : -1,
    'recombination' : 3,
    'flexible_transfinite' : False
}
geo_param={
    'a_i':a_i,'b_i':b_i,
    'a_m':a_m,'b_m':b_m,
    'a_o':a_o,'b_o':b_o,
    'a_d':a_d,'b_d':b_d, 
    'phi':phi[0],
    'R_PML': 1e-1,  # R_PML/R radius of PML region
    'R_TR': 1e-20, # R_TR/R radius of truncated region
    'x_ofst':x_offset[0], 'y_ofst': y_offset[0],
    'Nmu': 1, # element in structured layer (1 recommended)
    'Nint': 100, # nodes along sign-changing interface
    'GeomProgint': 1.01, # geometric progression (larger than 1)
    'Nz': 50, # element in z-direction in corner region
    'CharLengthMin_adim': 1, # Characteristic mesh size (ratio)
    'CharLengthMax_adim': 5,
    'GenerateQuadMesh': 1
}
    # Mesh generation, assembly, and solve
gmshfile = build_mesh(geofile,geo_param,gmsh_param)
PEP = build_PEP(gmshfile,name,
              a_m,b_m,a_d,b_d,
              phi,corner_pos, x_offset, y_offset,
              tag_names)
PEP.set_exact_eigenvalues(kappa_ex_even_fun,kappa_ex_odd_fun)
V_fem = ("CG", gmsh_param['order'])
PEP.assemble(V_fem)
alpha = np.exp(1j*np.pi/5); # PML angle
PEP.solve(SLEPc_params,alpha=alpha,tolReIm=1e-6)
eta = np.linspace(1e-16,8,200); # in (0,inf)
f = plt.figure()
ax = PEP.plot_spectrum(f,eta,plot_mesh_name=True)
ax.set_xlim([-1,-0.2])
ax.set_ylim([-0.17,0.1])
ax.xaxis.set_major_locator(plt.MultipleLocator(0.1))
ax.xaxis.set_major_formatter('{x:.1f}')
ax.grid(b=True, axis='x',which='major', color='gray', linestyle='-')
ax.xaxis.set_minor_locator(plt.MultipleLocator(0.01))
ax.minorticks_on()
ax.grid(b=True, which='minor', color='gray', linestyle='-', alpha=0.4)
PEP.export_eigenfunction_pyvista(f"{name}-{os.path.basename(geofile)}")

    # -- Geometrical file: with buffer layer (R_PML<R), where the buffer layer
    # is discretized in Cartesian coordinates. This delivers poor results. 
geofile=os.path.join(DIR_MESH,"Ellipse-with-1-corner_Structured-v2.geo")
gmshfile=os.path.join(os.getcwd(),"Ellipse-with-1-corner_Structured-v2.msh")
tag_names = {
    'omega_m': ['omega-m'],
    'omega_d': ['omega-d'],
    'gamma_d': ['gamma-d'],
    'omega_m_pml': ['corner-omega-m-pml'],
    'omega_d_pml': ['corner-omega-d-pml'],
    'Euler_bottom_bnd': ['corner-bnd-bot'],
    'Euler_right_bnd': ['corner-bnd-right'],
}
gmsh_param = {
    'save_and_exit': True,
    'binary': True,
    'order': 2,
    'meshing': -1,
    'recombination': -1,
    'flexible_transfinite': True
}
geo_param={
    'a_i':a_i,'b_i':b_i,
    'a_m':a_m,'b_m':b_m,
    'a_o':a_o,'b_o':b_o,
    'a_d':a_d,'b_d':b_d, 
    'phi1':phi[0],
    'R1_PML': 5e-1,  # R_PML/R radius of PML region
    'R1_TR': 1e-20, # R_TR/R radius of truncated region
    'x_ofst':x_offset[0], 'y_ofst': y_offset[0],
    'Nmu': 1, # element in structured layer (1 recommended)
    'Nint': 100, # nodes along sign-changing interface
    'Nint_corner': 1,
    'GeomProgint': 1 - 0.1, # geometric progression (smaller than 1)
    'Nz': 20, # element in z-direction in corner region
    'CharLengthMin_adim': 1, # Characteristic mesh size (ratio)
    'CharLengthMax_adim': 2,
    'GenerateQuadMesh': 1
}
    # Mesh generation, assembly, and solve
gmshfile = build_mesh(geofile,geo_param,gmsh_param)
PEP = build_PEP(gmshfile,name,
              a_m,b_m,a_d,b_d,
              phi,corner_pos, x_offset, y_offset,
              tag_names)
PEP.set_exact_eigenvalues(kappa_ex_even_fun,kappa_ex_odd_fun)
V_fem = ("CG", gmsh_param['order'])
PEP.assemble(V_fem)
alpha = np.exp(1j*np.pi/5); # PML angle
PEP.solve(SLEPc_params,alpha=alpha,tolReIm=1e-6)
eta = np.linspace(1e-16,8,200); # in (0,inf)
f = plt.figure()
ax = PEP.plot_spectrum(f,eta,plot_mesh_name=True)
ax.set_xlim([-1,-0.2])
ax.set_ylim([-0.17,0.1])
ax.xaxis.set_major_locator(plt.MultipleLocator(0.1))
ax.xaxis.set_major_formatter('{x:.1f}')
ax.grid(b=True, axis='x',which='major', color='gray', linestyle='-')
ax.xaxis.set_minor_locator(plt.MultipleLocator(0.01))
ax.minorticks_on()
ax.grid(b=True, which='minor', color='gray', linestyle='-', alpha=0.4)
PEP.export_eigenfunction_pyvista(f"{name}-{os.path.basename(geofile)}")
#%% Case B from JCP
name = "JCP Case B"
    # -- Case definition
(a_m,b_m,a_d,b_d,x_c,y_c,x_m,y_m,phi,corner_pos) = PEP_ana.perturbed_ellipse_parameters("JCP2021_case_B")
c = np.sqrt(a_m**2-b_m**2)
    # Mesh: structured layer (outer/inner)
fac = 0.15
a_o = (1-fac)*a_m+fac*a_d; b_o = np.sqrt(a_o**2-c**2)
a_i = c*np.cosh(2*np.arccosh(a_m/c) - np.arccosh(a_o/c))
b_i = np.sqrt(a_i**2-c**2)
    # Corner region offset
x_offset = [-1.40*a_d]
y_offset = [np.pi-b_d]
    # Complex resonances computed with COMSOL
kappa_ex_even_fun = lambda N : [-0.4076788197626322+1j*3.146916301514629e-8, -0.7047442590370551-0.054283845652155514*1j,-0.8655484713006019-0.017652368661287627*1j]
kappa_ex_odd_fun = lambda N : []
    # -- Mesh tags
tag_names = {
    'omega_m': ['omega-m','corner-omega-m'],
    'omega_d': ['omega-d','corner-omega-d'],
    'gamma_d': ['gamma-d'],
    'omega_m_pml': ['corner-omega-m-pml'],
    'omega_d_pml': ['corner-omega-d-pml'],
    'Euler_bottom_bnd': ['corner-bnd-bot'],
    'Euler_right_bnd': ['corner-bnd-right'],
}
    # -- Geometrical file: with buffer layer, i.e. the PML is part of the Euler
    # region (R_PML<R).
geofile=os.path.join(DIR_MESH,"Ellipse-with-1-corner_Structured-v3.geo")
gmshfile=os.path.join(os.getcwd(),"Ellipse-with-1-corner_Structured-v3.msh")
gmsh_param = {
    'save_and_exit': True,
    'binary': True,
    'order': 2,
    'meshing': -1,
    'recombination': 3,
    'flexible_transfinite': False
}
geo_param={
    'a_i':a_i,'b_i':b_i,
    'a_m':a_m,'b_m':b_m,
    'a_o':a_o,'b_o':b_o,
    'a_d':a_d,'b_d':b_d, 
    'phi':phi[0],
    'R_PML': 1e-1,  # R_PML/R radius of PML region
    'R_TR': 1e-20, # R_TR/R radius of truncated region
    'x_ofst': x_offset[0], 'y_ofst': y_offset[0],
    'Nmu': 1, # element in structured layer (1 recommended)
    'Nint': 100, # nodes along sign-changing interface
    'GeomProgint': 1.01, # geometric progression (larger than 1)
    'Nz': 50, # element in z-direction in corner region
    'CharLengthMin_adim': 1, # Characteristic mesh size (ratio)
    'CharLengthMax_adim': 5,
    'GenerateQuadMesh': 1
}
    # Mesh generation, assembly, and solve
gmshfile = build_mesh(geofile,geo_param,gmsh_param)
PEP = build_PEP(gmshfile,name,
              a_m,b_m,a_d,b_d,
              phi,corner_pos, x_offset, y_offset,
              tag_names)
PEP.set_exact_eigenvalues(kappa_ex_even_fun,kappa_ex_odd_fun)
V_fem = ("CG", gmsh_param['order'])
PEP.assemble(V_fem)
alpha = np.exp(1j*np.pi/5); # PML angle
PEP.solve(SLEPc_params,alpha=alpha,tolReIm=1e-6)
eta = np.linspace(1e-16,8,200); # in (0,inf)
f = plt.figure()
ax = PEP.plot_spectrum(f,eta,plot_mesh_name=True)
ax.set_xlim([-1,-0.2])
ax.set_ylim([-0.17,0.1])
ax.xaxis.set_major_locator(plt.MultipleLocator(0.1))
ax.xaxis.set_major_formatter('{x:.1f}')
ax.grid(b=True, axis='x',which='major', color='gray', linestyle='-')
ax.xaxis.set_minor_locator(plt.MultipleLocator(0.01))
ax.minorticks_on()
ax.grid(b=True, which='minor', color='gray', linestyle='-', alpha=0.4)
PEP.export_eigenfunction_pyvista(f"{name}-{os.path.basename(geofile)}")