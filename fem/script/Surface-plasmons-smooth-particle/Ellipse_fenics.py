#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plasmonic Eigenvalue Problem (PEP) for an elliptical particle.

Discretization using conformal and non-conformal approximation spaces.

"""

import numpy as np
import matplotlib.pyplot as plt

import fenics as fe
from petsc4py import PETSc
fe.SubSystemsManager.init_petsc()
    # Ensure errors are properly reported
PETSc.Sys.pushErrorHandler("python")
import SLEPc_utils
from slepc4py import SLEPc
import gmsh_utils
import itertools
import meshio_utils
import fenics_utils
import os
import PlasmonicEigenvalueProblem_utils_fenics as PEP_utils
from PlasmonicEigenvalueProblem_utils_fenics import DIR_MESH
import PlasmonicEigenvalueProblem_utils_analytical as PEP_ana

def build_mesh(geofile,geo_param,gmsh_param):
    """ Build mesh file from geometry file. """
    gmshfile = os.path.splitext(geofile)[0]+'.msh'
    gmsh_utils.generate_mesh_cli(geofile,gmshfile,2,parameters=geo_param,**gmsh_param)
    return gmshfile

    # Common
OptDB = PETSc.Options()
OptDB["st_ksp_type"] = "preonly"
OptDB["st_pc_type"] = "lu"
OptDB["st_pc_factor_mat_solver_type"] = "umfpack"
SLEPc_params = {
    'nev': 50,
    'target': 0.75,
    'shift': 0.75,
    'problem_type': SLEPc.EPS.ProblemType.GNHEP,
    'solver': SLEPc.EPS.Type.KRYLOVSCHUR,
    'tol': 1e-16,
    'max_it': 1000}
(a_m,b_m,a_d,b_d) = PEP_ana.ellipse_parameters("JCP2021_figure_10")
c = np.sqrt(a_m**2-b_m**2)
kappa_ex_ev = PEP_ana.ellipse_eigenvalues(a_m, b_m, a_d, b_d, 10)[0]
kappa_ex_odd = PEP_ana.ellipse_eigenvalues(a_m, b_m, a_d, b_d, 10)[1]

    # Mesh with symmetric layer at sign-changing interface
geofile=os.path.join(DIR_MESH,"Ellipse_Structured.geo")
a_o = 1.02*a_m; b_o = np.sqrt(a_o**2-c**2)
a_i = c*np.cosh(2*np.arccosh(a_m/c) - np.arccosh(a_o/c))
b_i = np.sqrt(a_i**2-c**2)
gmsh_param = {
    'save_and_exit': True,
    'binary': True,
    'order' : 1,
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
    'GenerateQuadMesh': 0,
    'TransfiniteSurface': 1
}

    # Mesh without explicit symmetric layer
    # Ellipse split in half
geofile=os.path.join(DIR_MESH,"Ellipse_Unstructured-v2.geo")
gmsh_param = {
    'save_and_exit': True,
    'binary': True,
    'order' : 1,
    'meshing' : -1,
    'recombination' : -1,
    'flexible_transfinite' : False
}
geo_param={
    'a_m':a_m,'b_m':b_m, # sign-changing interface
    'a_d':a_d,'b_d':b_d, # Dirichlet ellipse
    'Nint':500,
    'GeomProgint':0.9, # smaller than one (refinment towards tip)
    'CharLengthMin_adim':1, # Characteristic mesh size (ratio)
    'CharLengthMax_adim':2,
    'GenerateQuadMesh': 0,
}

    # Load and plot mesh (dolfin)
gmshfile = build_mesh(geofile,geo_param,gmsh_param)
xdmfiles = meshio_utils.convert_gmsh_to_XDMF(gmshfile,prune_z=True)
dmesh = fenics_utils.DolfinMesh.init_from_xdmf(xdmfiles['triangle'],xdmfiles['line'])
gmsh_phys = meshio_utils.read_phys_gmsh(xdmfiles['gmsh_physical_entities'])
PhysName_1D_tag = gmsh_phys[0]; PhysName_2D_tag = gmsh_phys[1]
fe.plot(dmesh.mesh,title=f"{os.path.basename(gmshfile)} \n {dmesh.mesh.num_vertices()} vertices")

    # Finite element assembly
# Assemble standard formulation
V = fe.FunctionSpace(dmesh.mesh, 'CR', 1) # 'CR' (best) or 'P'
(A_m,A_d) = PEP_utils.assemble_PEP_no_PML(V, dmesh, 
                        PhysName_2D_tag['omega-m'], 
                        PhysName_2D_tag['omega-d'], 
                        PhysName_1D_tag['gamma-d'],
                        diag_Am=1e2, diag_Ad=1e-2)
#% Assemble mixed formulation
BDM = fe.FiniteElement("BDM", dmesh.mesh.ufl_cell(), 1)
DG  = fe.FiniteElement("DG", dmesh.mesh.ufl_cell(), 1)  
# W = fe.FunctionSpace(dmesh.mesh, BDM * DG)
# (A_m,A_d) = PEP_utils.assemble_PEP_no_PML_mixed(W,dmesh, 
#                         PhysName_2D_tag['omega-m'], 
#                         PhysName_2D_tag['omega-d'], 
#                         PhysName_1D_tag['gamma-d'],
#                         diag_Am=1e2, diag_Ad=1e-2)
#% Assemble DG formulation based on BZ Primal
(V,A_m,A_d) = PEP_utils.assemble_PEP_no_PML_DG_BZ(1,dmesh,
                         PhysName_2D_tag['omega-m'], 
                         PhysName_2D_tag['omega-d'], 
                         PhysName_1D_tag['gamma-d'],
                         diag_Am=1e2, diag_Ad=1e-2)
#% Assemble DG formulation based on BZ mixed
# (W,A_m,A_d) = PEP_utils.assemble_PEP_no_PML_DG_BZM(1,dmesh,
#                         PhysName_2D_tag['omega-m'], 
#                         PhysName_2D_tag['omega-d'], 
#                         PhysName_1D_tag['gamma-d'],
#                         diag_Am=1e2, diag_Ad=1e-2)
#% Assemble DG formulation based on Local DG mixed
# (W,A_m,A_d) = PEP_utils.assemble_PEP_no_PML_DG_LDG(1,dmesh,
#                         PhysName_2D_tag['omega-m'], 
#                         PhysName_2D_tag['omega-d'], 
#                         PhysName_1D_tag['gamma-d'],
#                         diag_Am=1e2, diag_Ad=1e-2)
    # Solve sparse eigenvalue problem:  using SLEPc
Am_petsc = A_m; Ad_petsc = A_d
EPS = SLEPc_utils.solve_GEP_shiftinvert(Am_petsc,Ad_petsc,
                          problem_type=SLEPc.EPS.ProblemType.GNHEP,
                          solver=SLEPc.EPS.Type.KRYLOVSCHUR,
                          nev=40,tol=1e-3,max_it=300,
                          target=0.4,shift=0.4)
(eigval,eigvec_r,eigvec_i) = SLEPc_utils.EPS_get_spectrum(EPS)
    # filter lambda=0
mask=[(np.abs(x)>1e-2) for x in eigval]
eigval = list(itertools.compress(eigval, mask))
eigvec_r = list(itertools.compress(eigvec_r, mask))
eigvec_i = list(itertools.compress(eigvec_i, mask))
# Plot (using a dolfin function)
i_plot = 7

eigvec_rp = fe.PETScVector(eigvec_r[i_plot])
#eigvec_rp=eigvec_r[i_plot]

f = plt.figure()
ax=f.add_subplot(2,1,1)
fun=fe.Function(V,eigvec_rp)
#fun=fe.Function(W,eigvec_rp/eigvec_rp.norm('linf'))
#(sigma,fun) = fun.split()

c=fe.plot(fun,axes=ax,title=f"$\lambda=${eigval[i_plot]:1.4g}",
          mode='color')
ax.set_xlabel("x")
ax.set_ylabel("y")
plt.colorbar(c)
ax=f.add_subplot(2,1,2)
ax.plot(np.real(kappa_ex_ev),np.imag(kappa_ex_ev),label=r'$\kappa_{\rm{exact}}^{\rm{even}}$',linestyle='none',marker='o',mfc='none')
ax.plot(np.real(kappa_ex_odd),np.imag(kappa_ex_odd),label=r'$\kappa_{\rm{exact}}^{\rm{odd}}$',linestyle='none',marker='o',mfc='none')
ax.plot(-np.real(eigval),-np.imag(eigval),label='FEM',linestyle='none',marker='x',mfc='none')
ax.plot(-np.real(eigval[i_plot]),-np.imag(eigval[i_plot]),linestyle='none',marker='o')
ax.set_xlabel(r"$\Re(\kappa)$")
ax.set_ylabel(r"$\Im(\kappa)$")
# ax.set_title(f"Crouzeix-Raviart (CR1) -- {A_m.size[0]} DoF")
# ax.set_title(f"Lagrange (P1) ALT -- {A_m.size[0]} DoF")
ax.set_title(f"Mixed -- {A_m.size[0]} DoF")
ax.set_xlim([-1,-0.2])
ax.minorticks_on()
#ax.grid(b=True,which='both')
ax.legend(loc="upper left",ncol=3)#%% Weak formulation
#%% Alternative formulation to obtain S.P.D. B matrix
# Does not yields better spectrum that the direct formulation
V = fe.FunctionSpace(dmesh.mesh, 'CR', 1) # 'CR' (best) or 'P'
    # Dirichlet boundary condition
u_D = fe.Constant(0.0)
bc = fe.DirichletBC(V, u_D, dmesh.boundaries,PhysName_1D_tag['gamma-d'][0])
g = fe.Constant(1.0+1j*1e0)
    # Variational forms
u = fe.TrialFunction(V)
v = fe.TestFunction(V)
a_m = fe.dot(fe.grad(u), fe.grad(v))*dmesh.dx(PhysName_2D_tag['omega-m'][0])-fe.dot(fe.grad(u), fe.grad(v))*dmesh.dx(PhysName_2D_tag['omega-d'][0])
a_d = g*fe.dot(fe.grad(u), fe.grad(v))*dmesh.dx(PhysName_2D_tag['omega-m'][0])+fe.dot(fe.grad(u), fe.grad(v))*dmesh.dx(PhysName_2D_tag['omega-d'][0])
A_m = fe.PETScMatrix(); A_d = fe.PETScMatrix()
fenics_utils.assemble_GEP(A_m,a_m,A_d,a_d,[bc],diag_A=1e2,diag_B=1e-2)
#%% Solve sparse eigenvalue problem:  using SLEPc
EPS = fenics_utils.create_GEP(A_m,A_d,opts={
    'spectrum': "target real",
    'solver': "arpack",
    'tolerance': 1e-4,
    'maximum_iterations': 100,
#    'problem_type': 'gen_hermitian',
    'spectral_transform' : 'shift-and-invert',
    'spectral_shift' : (0.5-1.0)/(0.5+1)
    })
EPS.solve(30)
print(f"No. of converged eigenvalues: {EPS.get_number_converged()}")
eigval = np.zeros((EPS.get_number_converged(),))
eigval = [EPS.get_eigenvalue(i)[0]+1j*EPS.get_eigenvalue(i)[1]
          for i in range(len(eigval))]
eigval=np.array(eigval)
eigval=(1+eigval)/(1-eigval)

eigvec_r = list(); eigvec_i = list()
for i in range(len(eigval)):
    eigvec_r.append(EPS.get_eigenpair(i)[2])
    eigvec_i.append(EPS.get_eigenpair(i)[3])