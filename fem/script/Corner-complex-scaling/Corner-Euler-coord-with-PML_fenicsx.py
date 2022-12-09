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

    # Geometrical parameters
geofile=os.path.join(DIR_MESH,"Corner-Euler-PML.geo")
gmshfile="Corner-Euler-PML.msh"
R=1; R_PML = 1e-1; R_TR =1e-10; # Radiuses
phi = [np.pi/2]; # Camember angle (rad)
Nr=int(1e2); # No. of element on (ln(R_TR),ln(R))
Ntheta=int(1e1) # No. of elements in (-pi,pi)
gmsh_param={'refinement': 0,'binary': True}
geo_param={
    'R':R,
    'R_PML':R_PML,
    'R_TR':R_TR,
    'phi':phi[0],
    'Nz':Nr,
    'Ntheta':Ntheta
}
    # Load and plot mesh
gmsh_utils.generate_mesh_cli(geofile,gmshfile,2,parameters=geo_param,**gmsh_param)
mesh_cells_type = meshio_utils.get_cells_type_str(gmshfile)
dmesh = fenicsx_utils.DolfinxMesh.init_from_gmsh(gmshfile,2)
phys_tag_1D = gmsh_utils.getPhysicalNames(gmshfile,1)
phys_tag_2D = gmsh_utils.getPhysicalNames(gmshfile,2)
    # Periodic boundary condition
pbc = PEP_utils.get_pbc_corner_fenicsx(phys_tag_1D['corner-bnd-bot'])
    # Physical entities
phys_tags = dict()
phys_tags['omega-m'] = phys_tag_2D['omega-m']
phys_tags['omega-d'] = phys_tag_2D['omega-d']
phys_tags['gamma-d'] = phys_tag_1D['corner-bnd-right']
phys_tags['omega-m-pml'] = phys_tag_2D['omega-m-pml'] 
phys_tags['omega-d-pml'] = phys_tag_2D['omega-d-pml']
    # Define PEP
PEP = PEP_utils.PEP_with_PML_fenicsx(dmesh,phys_tags,pbc=pbc,split_real_imag=True)
PEP.set_deformed_critical_interval(phi)
    # Assemble
V_fem = ("CG", 2)
PEP.assemble(V_fem)
    # Solve PEP
alpha = np.exp(1j*np.pi/8); # PML angle
SLEPc_params = {'nev': 30,
             'target': 0.5,
             'shift': 0.5,
             'problem_type': SLEPc.EPS.ProblemType.GNHEP,
             'solver': SLEPc.EPS.Type.KRYLOVSCHUR,
             'tol': 1e-2,
             'max_it': 1000}
PEP.solve(SLEPc_params,alpha=alpha,tolReIm=1e-6)
    # Plot spectrum
eta = np.linspace(1e-16,8,200); # in (0,inf)
f = plt.figure()
ax = PEP.plot_spectrum(f,eta)
PEP.export_eigenfunction_pyvista(f"corner-Euler-{np.rad2deg(phi[0]):1.2g}deg")