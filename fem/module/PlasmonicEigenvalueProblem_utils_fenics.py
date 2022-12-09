#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions for solving the plasmonic eigenvalue problem (PEP) using fenics.

"""

import numpy as np
from petsc4py import PETSc
import SLEPc_utils
import warnings
import os
import fenics as fe
import fenics_utils
import PlasmonicEigenvalueProblem_utils_analytical as PEP_utils_ana

    # Module directory
DIR_MODULE = os.path.dirname(__file__)
    # Directory containing meshes
DIR_MESH = os.path.abspath(os.path.join(DIR_MODULE, os.pardir))
DIR_MESH = os.path.join(DIR_MESH,"mesh")


def assemble_PEP_no_PML(V, dmesh,
                        tags_omegaM, tags_omegaD, tags_Dirichlet,
                        diag_Am=1e2, diag_Ad=1e-2):
    """
    Assemble matrices for the PEP on a particle without corners:
        .. math:: \int_{\Omega_m}{\nabla u \cdot \nabla v}
                = lambda * \int_{\Omega_d}{\nabla u \cdot \nabla v},
    with Dirichlet boundary conditions.
    
    Parameters
    ----------
    V : dolfin.function.functionspace.FunctionSpace
        Function space.
    dmesh : fenics_utils.DolfinMesh
    tags_omegaM, tags_OmegaD : list of int
        Tags that identity the domains that constitutes Omega_m (Omega_d).
    tags_Dirichlet : list of int
        Tags that identify the Dirichlet boundary.
    diag_Am, diag_Ad : float, optional
        Diagonal penalty for enforcing Dirichlet boundary condition

    Returns
    -------
    A_m, A_d : dolfin.cpp.la.PETScMatrix
        Assembled matrices.

    """
        # surface measure
        # TODO: not needed, already in dmesh
    dx = fe.Measure("dx",domain=dmesh.mesh,subdomain_data=dmesh.domains)
        # Dirichlet boundary condition
    u_D = fe.Constant(0.0)
    bcs=[]
    for tag in tags_Dirichlet:
        bcs.append(fe.DirichletBC(V, u_D, dmesh.boundaries,tag))
        # Bilinear form over Omega M
    u = fe.TrialFunction(V)
    v = fe.TestFunction(V)
    for i in range(len(tags_omegaM)):
        if (i==0): # TODO: there must be a way to initialize a_m to avoid this
            a_m = fe.dot(fe.grad(u), fe.grad(v))*dx(tags_omegaM[i])
        else:
            a_m += fe.dot(fe.grad(u), fe.grad(v))*dx(tags_omegaM[i])
        # Bilinear form over Omega D
    for i in range(len(tags_omegaD)):
        if (i==0):
            a_d = fe.dot(fe.grad(u), fe.grad(v))*dx(tags_omegaD[i])
        else:
            a_d += fe.dot(fe.grad(u), fe.grad(v))*dx(tags_omegaD[i])
        # Assemble sparse eigenvalue problem
    A_m = fe.PETScMatrix(); A_d = fe.PETScMatrix()
    fenics_utils.assemble_GEP(A_m,a_m,A_d,a_d,bcs,diag_A=diag_Am,diag_B=diag_Ad)
    return A_m.mat(), A_d.mat()

def assemble_PEP_no_PML_mixed(W,dmesh,
                        tags_omegaM, tags_omegaD, tags_Dirichlet,
                        diag_Am=1e2, diag_Ad=1e-2):
    """
    Same as assemble_PEP_no_PML_fenics but using a mixed formulation
    H(div) x L^2.
    """
    dx = fe.Measure("dx",domain=dmesh.mesh,subdomain_data=dmesh.domains)
        # Dirichlet boundary condition
    bcs=[]
        # Define trial and test functions
    (phi, u) = fe.TrialFunctions(W)
    (psi, v) = fe.TestFunctions(W)
        # Right-hand side
    a_lhs, a_rhs = list(), list()      
    a_lhs.append( fe.dot(u,fe.div(psi))*dx)
    a_lhs.append( fe.dot(fe.div(phi),v)*dx)
    for i in range(len(tags_omegaM)):
        a_rhs.append(fe.dot(phi, psi)*dx(tags_omegaM[i]))
    for i in range(len(tags_omegaD)):
        a_lhs.append(fe.dot(phi, psi)*dx(tags_omegaD[i]))
        # Assemble sparse eigenvalue problem
    A_lhs = fe.PETScMatrix(); A_rhs = fe.PETScMatrix()
    fenics_utils.assemble_GEP(A_lhs,sum(a_lhs),A_rhs,sum(a_rhs),bcs,diag_A=diag_Am,diag_B=diag_Ad)
    return A_lhs.mat(), A_rhs.mat()

def assemble_PEP_no_PML_DG_BZ(p,dmesh,
                        tags_omegaM, tags_omegaD, tags_Dirichlet,
                        diag_Am=1e2, diag_Ad=1e-2):
    """
    Same as assemble_PEP_no_PML_fenics but using a DG formulation based
    on the Babuska-Zlamal method (primal). p is the polynomial degree.
    """
    
    dx = fe.Measure("dx",domain=dmesh.mesh,subdomain_data=dmesh.domains)
    dS = fe.Measure("dS",domain=dmesh.mesh,subdomain_data=dmesh.domains)
    ds = fe.Measure("ds",domain=dmesh.mesh,subdomain_data=dmesh.domains)
        # Define normal vector and mesh size
    n = fe.FacetNormal(dmesh.mesh)
    h = fe.CellDiameter(dmesh.mesh)
    h_avg = (h('+') + h('-'))/2
    eta, eta_b = 3, 1
    bcs = []
    spen = lambda h: h**(-2*p) # super-penalization
    V = fe.FunctionSpace(dmesh.mesh, 'DG', p)
    u, v = fe.TrialFunction(V), fe.TestFunction(V)
    from fenics import dot, grad, jump    
    a_lhs, a_rhs = list(), list()      
    if eta!=0:
        a_lhs.append(spen(h_avg)*eta/h_avg*dot(jump(u, n),jump(v, n))*dS)
    for tag in tags_Dirichlet:
        a_lhs.append(spen(h)*(eta_b/h)*u*v*ds(tag))
    for tag in tags_omegaD:
        a_lhs.append(dot(grad(u), grad(v))*dx(tag))
    for tag in tags_omegaM:
        a_rhs.append(dot(grad(u), grad(v))*dx(tag))
        # Assemble sparse eigenvalue problem
    A_lhs = fe.PETScMatrix(); A_rhs = fe.PETScMatrix()
    fenics_utils.assemble_GEP(A_lhs,sum(a_lhs),A_rhs,sum(a_rhs),bcs,diag_A=diag_Am,diag_B=diag_Ad)
    return V, A_lhs.mat(), A_rhs.mat()

def assemble_PEP_no_PML_DG_BZM(p,dmesh,
                        tags_omegaM, tags_omegaD, tags_Dirichlet,
                        diag_Am=1e2, diag_Ad=1e-2):
    """
    Same as assemble_PEP_no_PML_fenics but using a DG formulation based
    on the Babuska-Zlamal method (mixed). p is the polynomial degree.
    """
    
    dx = fe.Measure("dx",domain=dmesh.mesh,subdomain_data=dmesh.domains)
    dS = fe.Measure("dS",domain=dmesh.mesh,subdomain_data=dmesh.domains)
    ds = fe.Measure("ds",domain=dmesh.mesh,subdomain_data=dmesh.domains)
        # Define normal vector and mesh size
    n = fe.FacetNormal(dmesh.mesh)
    h = fe.CellDiameter(dmesh.mesh)
    h_avg = (h('+') + h('-'))/2
    eta, eta_b = 3, 1
    bcs = []
    spen = lambda h: h**(-2*p) # super-penalization
    DGv = fe.VectorElement("DG", dmesh.mesh.ufl_cell(), degree=p,dim=2)
    DG  = fe.FiniteElement("DG", dmesh.mesh.ufl_cell(), degree=p)
    V = fe.FunctionSpace(dmesh.mesh, DGv * DG)
    (phi, u) = fe.TrialFunctions(V)
    (psi, v) = fe.TestFunctions(V)    
    from fenics import dot, grad, jump
    a_lhs, a_rhs = list(), list()
    a_rhs.append(dot(phi,grad(v))*dx - dot(grad(u),psi) *dx)
    if eta!=0:
        a_rhs.append(spen(h_avg)*eta/h_avg*dot(jump(u, n),jump(v, n))*dS)
    for tag in tags_Dirichlet:
        a_rhs.append(spen(h)*(eta_b/h)*u*v*ds(tag))
    for tag in tags_omegaD:
        a_rhs.append(dot(phi, psi)*dx(tag))
    for tag in tags_omegaM:
        a_lhs.append(dot(phi, psi)*dx(tag))
        # Assemble sparse eigenvalue problem
    A_lhs = fe.PETScMatrix(); A_rhs = fe.PETScMatrix()
    fenics_utils.assemble_GEP(A_lhs,sum(a_lhs),A_rhs,sum(a_rhs),bcs,diag_A=diag_Am,diag_B=diag_Ad)
    return V, A_lhs.mat(), A_rhs.mat()


def assemble_PEP_no_PML_DG_LDG(p,dmesh,
                        tags_omegaM, tags_omegaD, tags_Dirichlet,
                        diag_Am=1e2, diag_Ad=1e-2):
    """
    Same as assemble_PEP_no_PML_fenics but using a DG formulation based
    on the Local DG method (mixed). p is the polynomial degree.
    """
    
    dx = fe.Measure("dx",domain=dmesh.mesh,subdomain_data=dmesh.domains)
    dS = fe.Measure("dS",domain=dmesh.mesh,subdomain_data=dmesh.domains)
    ds = fe.Measure("ds",domain=dmesh.mesh,subdomain_data=dmesh.domains)
        # Define normal vector and mesh size
    n = fe.FacetNormal(dmesh.mesh)
    h = fe.CellDiameter(dmesh.mesh)
    h_avg = (h('+') + h('-'))/2
    eta, eta_b = 1e3, 1e4
    beta = fe.as_vector([0,0])
    # beta = n('+')/2
    bcs = []
    DGv = fe.VectorElement("DG", dmesh.mesh.ufl_cell(), degree=p,dim=2)
    DG  = fe.FiniteElement("DG", dmesh.mesh.ufl_cell(), degree=p)
    V = fe.FunctionSpace(dmesh.mesh, DGv * DG)
    (phi, u) = fe.TrialFunctions(V)
    (psi, v) = fe.TestFunctions(V)    
    from fenics import dot, grad, jump , avg   
    a_lhs, a_rhs = list(), list()
    a_rhs.append(dot(phi,grad(v))*dx - dot(grad(u),psi) *dx)
    a_rhs.append(dot(jump(u, n),avg(psi))*dS - dot(avg(phi),jump(v, n))*dS)
    if eta!=0:
        a_rhs.append((eta/h_avg)*dot(jump(u, n),jump(v, n))*dS)
    if beta[0]!=0 or beta[1]!=0:
        a_rhs.append(dot(jump(u,n),beta*jump(psi,n))*dS - dot(beta*jump(phi,n),jump(v,n))*dS)
    for tag in tags_Dirichlet:
        a_rhs.append(dot(u*n,psi)*ds(tag) - dot(phi,v*n)*ds(tag) + (eta_b/h)*u*v*ds(tag))
    for tag in tags_omegaD:
        a_rhs.append(dot(phi, psi)*dx(tag))
    for tag in tags_omegaM:
        a_lhs.append(dot(phi, psi)*dx(tag))
        # Assemble sparse eigenvalue problem
    A_lhs = fe.PETScMatrix(); A_rhs = fe.PETScMatrix()
    fenics_utils.assemble_GEP(A_lhs,sum(a_lhs),A_rhs,sum(a_rhs),bcs,diag_A=diag_Am,diag_B=diag_Ad)
    return V, A_lhs.mat(), A_rhs.mat()

def assemble_PEP_PML(FE,dmesh,per,
                        tags_omegaM, tags_omegaD, tags_Dirichlet,
                        tags_omegaM_PML,tags_omegaD_PML,
                        diag_Am=1e2, diag_Ad=1e-2):
        # surface measure
    dx = fe.Measure("dx",domain=dmesh.mesh,subdomain_data=dmesh.domains)
        # Mixed FE space
    print(f"Number of periodic boundary conditions: {per.N}")
    if (per.N == 0):
        V = fe.FunctionSpace(dmesh.mesh,fe.MixedElement((FE,FE)))
    else:
        V = fe.FunctionSpace(dmesh.mesh,fe.MixedElement((FE,FE)),constrained_domain=per)
    Vr = V.sub(0); Vi = V.sub(1)
        # Dirichlet boundary condition
    u_D = fe.Constant(0.0)
    bcs=[]
    for tag in tags_Dirichlet:
        bcs.append(fe.DirichletBC(Vr, u_D, dmesh.boundaries,tag))
        bcs.append(fe.DirichletBC(Vi, u_D, dmesh.boundaries,tag))
        # Bilinear forms 
    def a_grad(u,v,tags):
        a_list=[]
        for i in range(len(tags)):
                a_list.append(fe.inner(fe.grad(u), fe.grad(v))*dx(tags[i]))
        return sum(a_list)
    
    def a_z(u,v,tags):
        a_list=[]
        for i in range(len(tags)):
                a_list.append(fe.inner(fe.Dx(u,0), fe.Dx(v,0))*dx(tags[i]))
        return sum(a_list)
        
    def a_theta(u,v,tags):
       a_list=[]
       for i in range(len(tags)):
               a_list.append(fe.inner(fe.Dx(u,1), fe.Dx(v,1))*dx(tags[i]))
       return sum(a_list)
       
    (ur,ui) = fe.TrialFunction(V)
    (vr,vi) = fe.TestFunction(V)
    a_m = [None]*5; a_d = [None]*5;
        # x 1
    a_m[0] = a_grad(ur,vr,tags_omegaM) + a_grad(ui,vi,tags_omegaM)
    a_d[0] = a_grad(ur,vr,tags_omegaD) + a_grad(ui,vi,tags_omegaD)
        # x (real(alpha))
    a_m[1] = a_z(ur,vr,tags_omegaM_PML) + a_z(ui,vi,tags_omegaM_PML)
    a_d[1] = a_z(ur,vr,tags_omegaD_PML) + a_z(ui,vi,tags_omegaD_PML)
        # x (imag(alpha))
    a_m[2] = -a_z(ui,vr,tags_omegaM_PML) + a_z(ur,vi,tags_omegaM_PML)
    a_d[2] = -a_z(ui,vr,tags_omegaD_PML) + a_z(ur,vi,tags_omegaD_PML)
        # x (real(1/alpha))
    a_m[3] = a_theta(ur,vr,tags_omegaM_PML) + a_theta(ui,vi,tags_omegaM_PML)
    a_d[3] = a_theta(ur,vr,tags_omegaD_PML) + a_theta(ui,vi,tags_omegaD_PML)        
        # x (imag(1/alpha))
    a_m[4] = -a_theta(ui,vr,tags_omegaM_PML) + a_theta(ur,vi,tags_omegaM_PML)
    a_d[4] = -a_theta(ui,vr,tags_omegaD_PML) + a_theta(ur,vi,tags_omegaD_PML)
        
        # Assemble sparse eigenvalue problem
    A_m = [None]*5; A_d = [None]*5;
    bcs_list = [ [] ]*5;
    bcs_list[0] = bcs; # Dirichlet b.c. only needed for 'x 1' term
    for i in range(len(A_m)):
        A_m[i] = fe.PETScMatrix(); A_d[i] = fe.PETScMatrix();
        fenics_utils.assemble_GEP(A_m[i],a_m[i],A_d[i],a_d[i],bcs_list[i],diag_A=diag_Am,diag_B=diag_Ad)
        # Indices of real and imaginary parts
    idx_ur = Vr.dofmap().dofs()
    idx_ui = Vi.dofmap().dofs()
    return A_m, A_d, idx_ur, idx_ui, V.sub(0), V.sub(1)

# Definition of periodic boundary conditions

def get_pbc_corner(y_offset):
    """
    Get periodic boundary conditions for a corner in Euler coordinates (z,theta).
    
    Parameters
    ----------
    y_offset : scalar
        y-offset of the computational domain.
        
    Returns
    -------
    per : fenics_utils.PeriodicBoundary
        Periodic boundary conditions
        
    """
    def EulerBot_map(x): # Bottom boundary {y-y_offset=-pi}
        return fe.near(x[1],y_offset-fe.DOLFIN_PI)
    def EulerTop_map(x): # Top boundary {y-y_offset=pi}
        return fe.near(x[1],y_offset+fe.DOLFIN_PI)
    def map_EulerBot_to_EulerTop(x,y):
        y[0] = x[0]; y[1] = x[1]+2*fe.DOLFIN_PI
    per = fenics_utils.PeriodicBoundary(); per.init()
    per.append(EulerBot_map,EulerTop_map,map_EulerBot_to_EulerTop)
    
    return per
      
def get_pbc_ellipse1corner(a_m,b_m,phi,R_TR,x_offset,y_offset,tol=1e-5):
    """
    Compute periodic boundary conditions for an ellipse perturbed by a corner, 
    assuming the ellipse-corner junction is C^1.
    
    Parameters
    ----------
    a_m, b_m : scalar
        ellipse m semi-axes
    phi : scalar
        corner angle
    R_TR: scalar
        truncation radius
    x_offset, y_offset : scalar
        x and y-offset of the corner region.
    tol: scalar
        tolerance
        
    Returns
    -------
    per : fenics_utils.PeriodicBoundary
        Periodic boundary conditions

    """
    per = fenics_utils.PeriodicBoundary(); per.init()
    (x_c,y_c,x_m,y_m,R) = PEP_utils_ana.get_C1corner_ellipse(a_m,b_m,phi)

    # bottom of Euler region (src) = top of Euler region (dst)
    def EulerBot_map(x): # Bottom boundary {Log(R_TR) < x-x_offset < Log(R), y-y_offset=-pi}
        return fe.near(x[1],y_offset-fe.DOLFIN_PI) and ((x[0]-x_offset-np.log(R))<tol) and ((x[0]-x_offset-np.log(R_TR))>-tol)
    def EulerTop_map(x): # Top boundary {Log(R_TR) < x-x_offset < Log(R), y-y_offset=pi}
        return fe.near(x[1],y_offset+fe.DOLFIN_PI) and ((x[0]-x_offset-np.log(R))<tol) and ((x[0]-x_offset-np.log(R_TR))>-tol)    
    def map_EulerBot_to_EulerTop(x,y):
        y[0] = x[0]; y[1] = x[1]+2*fe.DOLFIN_PI
    per.append(EulerBot_map,EulerTop_map,map_EulerBot_to_EulerTop)
    
    # right of Euler region (dst) = corner circle (src)
    def Circle_map(x): # Corner circle {|(x,y)-(x_c,0)|=R}
        if (fe.near((x[0]-x_c)**2+x[1]**2,R**2)):
            print(f"[is_circle] Point ({x[0]},{x[1]}) Matched")
        return fe.near((x[0]-x_c)**2+x[1]**2,R**2)
    def EulerRight_map(x): # rectangle right {x-x_offset=log(R), -pi<y-y_offset<pi}
        return fe.near(x[0],np.log(R)+x_offset) and ((x[1]-y_offset-np.pi)<tol) and ((x[1]-y_offset+np.pi)>-tol)
    def map_Circle_to_EulerRight(x,y):
        y[0] = 0.5*np.log((x[0]-x_c)**2+x[1]**2)+x_offset; # ln(r(x,y))+x_offset
        y[1] = np.arctan2(x[1],x[0]-x_c)+y_offset; # theta(x,y)+y_offset
    per.append(Circle_map,EulerRight_map,map_Circle_to_EulerRight)

    return per