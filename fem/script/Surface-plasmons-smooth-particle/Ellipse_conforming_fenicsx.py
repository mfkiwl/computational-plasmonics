#!/usr/bin/env python
# coding: utf-8
"""
Plasmonic Eigenvalue Problem (PEP) for an elliptical particle.

This script illustrates the influence of the local mesh symmetry around the
sign-changing interface. The weaker the symmetry the worse the spectral
pollution.

"""

import numpy as np
from petsc4py import PETSc
from slepc4py import SLEPc
import fenicsx_utils
import gmsh_utils
import meshio_utils
import matplotlib.pyplot as plt
import os
import PlasmonicEigenvalueProblem_utils_fenicsx as PEP_utils
from PlasmonicEigenvalueProblem_utils_fenicsx import DIR_MESH
import PlasmonicEigenvalueProblem_utils_analytical as PEP_ana

def build_mesh(geofile,geo_param,gmsh_param):
    """ Build mesh file from geometry file. """
    gmshfile = os.path.splitext(geofile)[0]+'.msh'
    gmsh_utils.generate_mesh_cli(geofile,gmshfile,2,parameters=geo_param,**gmsh_param)
    return gmshfile

def build_PEP(gmshfile,a_m,b_m,a_d,b_d):
    """ Instanciate a PEP object from mesh file. """
        # Mesh importation
    dmesh = fenicsx_utils.DolfinxMesh.init_from_gmsh(gmshfile,2)
    phys_tag_1D = gmsh_utils.getPhysicalNames(gmshfile,1)
    phys_tag_2D = gmsh_utils.getPhysicalNames(gmshfile,2)
        # Physical entities
    phys_tags = dict()
    phys_tags['omega-m'] = phys_tag_2D['omega-m']
    phys_tags['omega-d'] = phys_tag_2D['omega-d']
    phys_tags['gamma-d'] = phys_tag_1D['gamma-d']
        # Define PEP
    PEP = PEP_utils.PEP_with_PML_fenicsx(dmesh,phys_tags)
        # Exact eigenvalues
    kappa_ex_ev_fun = lambda N : PEP_ana.ellipse_eigenvalues(a_m, b_m, a_d, b_d, N)[0]
    kappa_ex_odd_fun = lambda N : PEP_ana.ellipse_eigenvalues(a_m, b_m, a_d, b_d, N)[1]
    PEP.set_exact_eigenvalues(kappa_ex_ev_fun,kappa_ex_odd_fun)
    return PEP

    # Common
V_fem = ("CG", 2)
SLEPc_params = {
    'nev': 50,
    'target': 0.75,
    'shift': 0.75,
    'problem_type': SLEPc.EPS.ProblemType.GNHEP,
    'solver': SLEPc.EPS.Type.KRYLOVSCHUR,
    'tol': 1e-16,
    'max_it': 1000}
OptDB = PETSc.Options()
OptDB["st_ksp_type"] = "preonly"
OptDB["st_pc_type"] = "lu"
OptDB["st_pc_factor_mat_solver_type"] = "mumps"

    # Geometrical case
#(a_m,b_m,a_d,b_d) = PEP_ana.ellipse_parameters("JCP2021_figure_9")
#(a_m,b_m,a_d,b_d) = PEP_ana.ellipse_parameters("JCP2021_figure_10")
(a_m,b_m,a_d,b_d) = PEP_ana.ellipse_parameters("JCP2021_case_A")
c = np.sqrt(a_m**2-b_m**2)

    # Mesh with symmetric layer at sign-changing interface
a_o = 1.03*a_m; b_o = np.sqrt(a_o**2-c**2)
a_i = c*np.cosh(2*np.arccosh(a_m/c) - np.arccosh(a_o/c))
b_i = np.sqrt(a_i**2-c**2)
geofile=os.path.join(DIR_MESH,"Ellipse_Structured.geo")
gmsh_param = {
    'save_and_exit': True,
    'binary': True,
    'order' : 2,
    'meshing' : 9,
    'recombination' : 3,
    'flexible_transfinite' : False
}
geo_param={
    'a_i':a_i,'b_i':b_i, # inner ellipse
    'a_m':a_m,'b_m':b_m, # sign-changing interface
    'a_o':a_o,'b_o':b_o, # outer ellipse
    'a_d':a_d,'b_d':b_d, # Dirichlet ellipse
    'N_mu':1,
    'Nint':2*10,
    'GeomProgint':1, # larger than one
    'CharLengthMin_adim':1, # Characteristic mesh size (ratio)
    'CharLengthMax_adim':2,
    'GenerateQuadMesh': 1,
    'TransfiniteSurface': 1
}
gmshfile = build_mesh(geofile,geo_param,gmsh_param)
PEP = build_PEP(gmshfile,a_m,b_m,a_d,b_d)
PEP.assemble(V_fem)
PEP.solve(SLEPc_params)
PEP.export_eigenfunction_pyvista(f"{os.path.basename(geofile)}")
PEP.export_eigenvalues_to_csv(f"{os.path.basename(geofile)}.csv")
f = plt.figure()
ax = PEP.plot_spectrum(f,N=20)
ax.set_ylim([-5e-2,5e-2])
ax.set_xlim([-1,-0.65])
ax.minorticks_on()
ax.set_aspect('equal') # orthonormal axis
mesh_cells_type = meshio_utils.get_cells_type_str(gmshfile)
title = f"Elliptical particle, a_d={a_d:2.2g}, a_m={a_m:2.2g}\n{V_fem}, {mesh_cells_type}, {PEP.dmesh.mesh.geometry.x.shape[0]} nodes \n {os.path.basename(geofile)}"
ax.set_title(title)

    # Mesh without explicit symmetric layer
    # Ellipse split in four
geofile=os.path.join(DIR_MESH,"Ellipse_Unstructured-v1.geo")
gmsh_param = {
    'save_and_exit': True,
    'binary': True,
    'order' : 2,
    'meshing' : 9,
    'recombination' : 3,
    'flexible_transfinite' : False
}
geo_param={
    'a_m':a_m,'b_m':b_m, # sign-changing interface
    'a_d':a_d,'b_d':b_d, # Dirichlet ellipse
    'Nint':130,
    'GeomProgint':1.02, # larger than one
    'CharLengthMin_adim':1, # Characteristic mesh size (ratio)
    'CharLengthMax_adim':2,
    'GenerateQuadMesh': 1, # Recombination can be challenging here
}
gmshfile = build_mesh(geofile,geo_param,gmsh_param)
PEP = build_PEP(gmshfile,a_m,b_m,a_d,b_d)
PEP.assemble(V_fem)
PEP.solve(SLEPc_params)
PEP.export_eigenfunction_pyvista(f"{os.path.basename(geofile)}")
PEP.export_eigenvalues_to_csv(f"{os.path.basename(geofile)}.csv")
f = plt.figure()
ax = PEP.plot_spectrum(f,N=20)
ax.set_ylim([-5e-2,5e-2])
ax.set_xlim([-1,-0.65])
ax.minorticks_on()
ax.set_aspect('equal') # orthonormal axis
mesh_cells_type = meshio_utils.get_cells_type_str(gmshfile)
title = f"Elliptical particle, a_d={a_d:2.2g}, a_m={a_m:2.2g}\n{V_fem}, {mesh_cells_type}, {PEP.dmesh.mesh.geometry.x.shape[0]} nodes \n {os.path.basename(geofile)}"
ax.set_title(title)

    # Mesh without explicit symmetric layer
    # Ellipse split in half
geofile=os.path.join(DIR_MESH,"Ellipse_Unstructured-v2.geo")
gmsh_param = {
    'save_and_exit': True,
    'binary': True,
    'order' : 2,
    'meshing' : -1,
    'recombination' : -1,
    'flexible_transfinite' : False
}
geo_param={
    'a_m':a_m,'b_m':b_m, # sign-changing interface
    'a_d':a_d,'b_d':b_d, # Dirichlet ellipse
    'Nint':150,
    'GeomProgint':0.9, # smaller than one (refinment towards tip)
    'CharLengthMin_adim':1, # Characteristic mesh size (ratio)
    'CharLengthMax_adim':2,
    'GenerateQuadMesh': 1,
}
gmshfile = build_mesh(geofile,geo_param,gmsh_param)
PEP = build_PEP(gmshfile,a_m,b_m,a_d,b_d)
PEP.assemble(V_fem)
PEP.solve(SLEPc_params)
PEP.export_eigenfunction_pyvista(f"{os.path.basename(geofile)}")
PEP.export_eigenvalues_to_csv(f"{os.path.basename(geofile)}.csv")
f = plt.figure()
ax = PEP.plot_spectrum(f,N=20)
ax.set_ylim([-5e-2,5e-2])
ax.set_xlim([-1,-0.65])
ax.minorticks_on()
ax.set_aspect('equal') # orthonormal axis
mesh_cells_type = meshio_utils.get_cells_type_str(gmshfile)
title = f"Elliptical particle, a_d={a_d:2.2g}, a_m={a_m:2.2g}\n{V_fem}, {mesh_cells_type}, {PEP.dmesh.mesh.geometry.x.shape[0]} nodes \n {os.path.basename(geofile)}"
ax.set_title(title)
PEP.export_eigenfunction_pyvista(f"{os.path.basename(geofile)}")
PEP.export_eigenvalues_to_csv(f"{os.path.basename(geofile)}.csv")