#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plasmonic Eigenvalue Problem (PEP) for an elliptical particle.

Discretization using T-conformal approximation, implemented using two auxiliary
variables.
"""
#%%
import numpy as np
import matplotlib.pyplot as plt

import fenics as fe
from petsc4py import PETSc
fe.SubSystemsManager.init_petsc()
import PETSc_utils
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
import multiphenics as mpfe
import multiphenics_utils

def build_mesh(geofile,geo_param,gmsh_param):
    """ Build mesh file from geometry file. """
    gmshfile = os.path.splitext(geofile)[0]+'.msh'
    gmsh_utils.generate_mesh_cli(geofile,gmshfile,2,parameters=geo_param,**gmsh_param)
    return gmshfile

    # Common
(a_m,b_m,a_d,b_d) = PEP_ana.ellipse_parameters("JCP2021_figure_10")
c = np.sqrt(a_m**2-b_m**2)
kappa_ex_ev = PEP_ana.ellipse_eigenvalues(a_m, b_m, a_d, b_d, 10)[0]
kappa_ex_odd = PEP_ana.ellipse_eigenvalues(a_m, b_m, a_d, b_d, 10)[1]

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
meshtag = meshio_utils.read_phys_gmsh(xdmfiles['gmsh_physical_entities'])
meshtag_1D = meshtag[0]; meshtag_2D = meshtag[1]
fe.plot(dmesh.mesh,title=f"{os.path.basename(gmshfile)} \n {dmesh.mesh.num_vertices()} vertices")
restriction_gamma =  multiphenics_utils.build_MeshRestriction_from_tags(dmesh.mesh,
                                                dmesh.boundaries,
                                                meshtag_1D['gamma-m'])
restriction_omega_m =  multiphenics_utils.build_MeshRestriction_from_tags_subdomains(dmesh.mesh,
                                                dmesh.domains,
                                                meshtag_2D['omega-m'])
#%% Finite element assembly for direct problem
# Standard H^1_0(Omega) space
degree = gmsh_param['order']
V_lag = fe.FunctionSpace(dmesh.mesh, 'P', 1) # 'CR' (best) or 'P'
#f_expr = 'exp(-x[0]*x[0]-x[1]*x[1])'
f_expr = 'cos(4*x[0]+4*x[1])'
f = fe.Expression(f_expr, element=V_lag.ufl_element())
sigma_1 = 1.0
sigma_2 = -0.9
a,b = list(), list()
u = fe.TrialFunction(V_lag)
v = fe.TestFunction(V_lag)
for tag in meshtag_2D["omega-d"]:
    a.append(sigma_1*fe.inner(fe.grad(u),fe.grad(v))*dmesh.dx(tag))    
for tag in meshtag_2D["omega-m"]:
    a.append(sigma_2*fe.inner(fe.grad(u),fe.grad(v))*dmesh.dx(tag))
b.append(fe.inner(f,v)*dmesh.dx)
u_fe = fe.Function(V_lag)
bcs = []
for tag in meshtag_1D["gamma-d"]:
    bcs.append(fe.DirichletBC(V_lag, fe.Constant(0.0), dmesh.boundaries,tag))
fe.solve(sum(a) == sum(b),u_fe,bcs,solver_parameters={"linear_solver": "lu"})
c = fe.plot(u_fe)
ax=c.axes
ax.set_xlabel("x")
ax.set_ylabel("y")
plt.colorbar(c)
#%% Finite element assembly for direct problem
# W = H^1_0(Omega) x H^1(Omega_m) x L^2(Gamma_m)
#   = V_v x V_w x V_l
W = mpfe.BlockFunctionSpace([V_lag, V_lag,V_lag], restrict=[None,
                                        restriction_omega_m,restriction_gamma])
(v, w, l) = mpfe.block_split(mpfe.BlockTrialFunction(W))
(phi_v, phi_w, phi_l) = mpfe.block_split(mpfe.BlockTestFunction(W))
a, b = list(), list()
# -- Simple system for debug
# for tag in meshtag_2D['omega-d']:
#     a.append([fe.inner(fe.grad(v), fe.grad(phi_v))*dmesh.dx(tag)])
# for tag in meshtag_2D['omega-m']:
#     a.append([fe.inner(fe.grad(v), fe.grad(phi_v))*dmesh.dx(tag)])
#     a.append([fe.inner(w, phi_w)*dmesh.dx(tag)])
#     b.append([fe.inner(f, phi_w)*dmesh.dx(tag)])
# for tag in meshtag_1D['gamma-m']:
#     a.append([fe.inner(l("-"), phi_l("-"))*dmesh.dS(tag)])
#     b.append([fe.inner(f, phi_l("-"))*dmesh.dS(tag)])
# b.append([fe.inner(f, phi_v)*dmesh.dx])
# --
# -- PDE on v
for tag in meshtag_2D['omega-d']:
    a.append([sigma_1*fe.inner(fe.grad(v), fe.grad(phi_v))*dmesh.dx(tag)])
for tag in meshtag_2D['omega-m']:
    a.append([-sigma_2*fe.inner(fe.grad(v), fe.grad(phi_v))*dmesh.dx(tag)])
    a.append([sigma_2*fe.inner(2*fe.grad(w), fe.grad(phi_v))*dmesh.dx(tag)])
b.append([fe.inner(f, phi_v)*dmesh.dx])
# -- w=v and null boundary Lagrange multiplier
# for tag in meshtag_2D['omega-m']:
#     a.append([fe.inner(w, phi_w)*dmesh.dx(tag)])
#     a.append([fe.inner(v, phi_w)*dmesh.dx(tag)])
# for tag in meshtag_1D['gamma-m']:
#     a.append([fe.inner(l("-"), phi_l("-"))*dmesh.dS(tag)])
# -- w is lifting operator
for tag in meshtag_2D['omega-m']:
    a.append([fe.inner(fe.grad(w), fe.grad(phi_w))*dmesh.dx(tag)])
for tag in meshtag_1D['gamma-m']:
    a.append([fe.inner(l("-"), phi_w("-"))*dmesh.dS(tag)])
    a.append([fe.inner(w("-"), phi_l("-"))*dmesh.dS(tag)])
    a.append([-fe.inner(v("-"), phi_l("-"))*dmesh.dS(tag)])
A, B = mpfe.block_assemble(a), mpfe.block_assemble(b)
bcs = []
for tag in meshtag_1D['gamma-d']:
    bcs.append(mpfe.DirichletBC(W.sub(0), fe.Constant(0.), dmesh.boundaries,tag))
bcs_block = mpfe.BlockDirichletBC(bcs)
bcs_block.apply(A)

U = mpfe.BlockFunction(W)
mpfe.block_solve(A, U.block_vector(), B)


# TODO: plot T(v) to enable comparison
plt.figure()
c = fe.plot(U[0]);
#c = fe.plot(charfun*U[0]+2*U[1],vmin=-0.2,vmax=0.2); # T(v)
ax=c.axes
ax.set_xlabel("x")
ax.set_ylabel("y")
plt.colorbar(c)
print(f"Total DoF = {A.size(0)}")
#xdmfile = 'poisson-equation.xdmf'
#sol_idx=[0,1]; sol_name = ["u","w"]
#multiphenics_utils.export_xdmf(xdmfile,sol_idx,sol_name,[0],[U.block_vector().vec()],W)



#%% Finite element assembly for eigenvalue problem
# Standard H^1_0(Omega) space
V = fe.FunctionSpace(dmesh.mesh, 'P', 1) # 'CR' (best) or 'P'
(A_m,A_d) = PEP_utils.assemble_PEP_no_PML(V, dmesh, 
                        meshtag_2D['omega-m'], 
                        meshtag_2D['omega-d'], 
                        meshtag_1D['gamma-d'],
                        diag_Am=1e2, diag_Ad=1e-2)

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
f = plt.figure()
ax = f.add_subplot(1, 1, 1)
ax.plot(-np.real(eigval), -np.imag(eigval), label='FEM',
        linestyle='none', marker='x', mfc='none')
ax.plot(np.real(kappa_ex_ev),np.imag(kappa_ex_ev),linestyle='none',marker='o',markerfacecolor='none',label="Exact")
ax.plot(np.real(kappa_ex_odd),np.imag(kappa_ex_odd),linestyle='none',marker='o',markerfacecolor='none',label="Exact")
ax.legend()
ax.set_xlabel(r"$\Re(\kappa)$")
ax.set_ylabel(r"$\Im(\kappa)$")
solver="direct"
ax.set_title(f"{solver} (N={eigvec_r[0].size})")

#%% W = H^1_0(Omega) x H^1(Omega_m) x L^2(Gamma_m)
#   = V_v x V_w x V_l
W = mpfe.BlockFunctionSpace([V_lag, V_lag,V_lag], restrict=[None,
                                        restriction_omega_m,restriction_gamma])
(v, w, l) = mpfe.block_split(mpfe.BlockTrialFunction(W))
(phi_v, phi_w, phi_l) = mpfe.block_split(mpfe.BlockTestFunction(W))
a,b = list(), list()
# -- PDE on v
for tag in meshtag_2D['omega-d']:
    a.append([sigma_1*fe.inner(fe.grad(v), fe.grad(phi_v))*dmesh.dx(tag)])
for tag in meshtag_2D['omega-m']:
    b.append([fe.inner(fe.grad(v), fe.grad(phi_v))*dmesh.dx(tag)])
    b.append([-fe.inner(2*fe.grad(w), fe.grad(phi_v))*dmesh.dx(tag)])
# -- w is lifting operator
for tag in meshtag_2D['omega-m']:
    a.append([fe.inner(fe.grad(w), fe.grad(phi_w))*dmesh.dx(tag)])
for tag in meshtag_1D['gamma-m']:
    a.append([fe.inner(l("-"), phi_w("-"))*dmesh.dS(tag)])
    a.append([fe.inner(w("-"), phi_l("-"))*dmesh.dS(tag)])
    a.append([-fe.inner(v("-"), phi_l("-"))*dmesh.dS(tag)])
A, B = mpfe.block_assemble(a), mpfe.block_assemble(b)
bcs = []
for tag in meshtag_1D['gamma-d']:
    bcs.append(mpfe.DirichletBC(W.sub(0), fe.Constant(0.), dmesh.boundaries,tag))
bcs_block = mpfe.BlockDirichletBC(bcs)
vec = fe.PETScVector(); vec.init(A.mat().size[0])
for bc in bcs:
    bc.zero_columns(A,vec,1e2)
    bc.zero_columns(B,vec,1e-2)

A_petsc = A.mat(); B_petsc = B.mat()
OptDB = PETSc_utils.get_cleared_options_db()
OptDB = PETSc.Options()
solver="direct"
if solver=="iterative":
    OptDB["st_ksp_type"] = "gmres" # fgmres, gmres, gcr, ibcgs
    OptDB["st_pc_type"] = "ilu" # ilu, bjacobi, icc
    OptDB["st_ksp_rtol"] = 1e-10
else:     # Direct solver (MUMPS)
    OptDB["st_ksp_type"] = "preonly"
    OptDB["st_pc_type"] = "lu"
    OptDB["st_pc_factor_mat_solver_type"] = "mumps"
SLEPc_params = {
    'nev': 50,
    'target': -1/0.4,
    'shift': -1/0.4,
    'problem_type': SLEPc.EPS.ProblemType.GNHEP,
    'solver': SLEPc.EPS.Type.KRYLOVSCHUR,
    'tol': 1e-6,
    'max_it': 1000}


EPS = SLEPc_utils.solve_GEP_shiftinvert(A_petsc,B_petsc,**SLEPc_params)
(eigval,eigvec_r,eigvec_i) = SLEPc_utils.EPS_get_spectrum(EPS)
    # filter lambda=0
mask=[(np.abs(x)>1e-2) for x in eigval]
eigval = list(itertools.compress(eigval, mask))
eigvec_r = list(itertools.compress(eigvec_r, mask))
eigvec_i = list(itertools.compress(eigvec_i, mask))

    # kappa = 1/sigma_2
eigval = [1/v for v in eigval]

# Plot eigenvalues
f = plt.figure()
ax=f.add_subplot(1,1,1)
ax.plot(np.real(kappa_ex_ev),np.imag(kappa_ex_ev),linestyle='none',marker='o',markerfacecolor='none',label="Exact")
ax.plot(np.real(kappa_ex_odd),np.imag(kappa_ex_odd),linestyle='none',marker='o',markerfacecolor='none',label="Exact")
ax.plot(np.real(eigval),np.imag(eigval),linestyle='none',marker='x',label="FEM")
ax.set_xlabel(r"$\Re(\lambda)$")
ax.set_ylabel(r"$\Im(\lambda)$")
#ax.set_xlim([np.min(np.real(eigval)),0])
#ax.set_xlim([-1.5,-0.5])
ax.set_title(f"{solver} (N={eigvec_r[0].size})")
ax.legend()
    # Paraview export of all eigenmodes
xdmfile = "laplace-eigenvalues.xdmf"
output = fe.XDMFFile(xdmfile)
output.parameters["flush_output"] = True
output.parameters["rewrite_function_mesh"] = False
output.parameters["functions_share_mesh"] = True
print(f"Export to {xdmfile}...")
for i in range(len(eigvec_r)):
    idx_u, idx_l = 0, 1
    eigvec_p = eigvec_r[i] / abs(eigvec_r[i]).max()[1]
    u_fun = multiphenics_utils.get_function(eigvec_p, idx_u,f"u_{eigval[i]}",W)
    output.write(u_fun, 0)
    u_fun = multiphenics_utils.get_function(eigvec_p, idx_l,f"l_{eigval[i]}",W)
    output.write(u_fun, 0)
