#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions for solving the plasmonic eigenvalue problem (PEP) using fenicsx.

"""

import numpy as np
from petsc4py import PETSc
import SLEPc_utils
import warnings
import os
import PlasmonicEigenvalueProblem_utils_analytical as PEP_utils_ana

    # Module directory
DIR_MODULE = os.path.dirname(__file__)
    # Directory containing meshes
DIR_MESH = os.path.abspath(os.path.join(DIR_MODULE, os.pardir))
DIR_MESH = os.path.join(DIR_MESH,"mesh")

    
try:
    import dolfinx
    import fenicsx_utils
    import ufl
except ImportError:
    warnings.warn('dolfinx, fenicsx_utils or ufl failed to import')

try:
    import dolfinx_mpc
    import dolfinx_mpc_utils
except ImportError:
    warnings.warn('dolfinx_mpc failed to import')

def assemble_PEP_PML_fenicsx(V,dmesh,tags,per,diag_Am=1e2,diag_Ad=1e-2,
                             split_real_imag=False):
    """
    Assemble matrices for the PEP. Handles: corner perturbation regions (Euler
    coordinates only), complex scaling regions (Euler coordinates only), and 
    Dirichlet boundary conditions.
        
    The generalized eigenvalue problem is assembled under one of the following
    two general forms:
        
        (i) If split_real_imag=False,

        (A_m0 + a*A_m1 + (1/a)*A_m2)*U
            = -kappa * (same with '_d' instead of '_m')*U,        
        
        (ii) If split_real_imag=True,

        (A_m0 + Re(a)*A_m1 + Im(a)*A_m2 + Re(1/a)*A_m3 + Im(1/a)*A_m4)*U
            = -kappa * (same with '_d' instead of '_m')*U,
            
    where all the matrices are real-valued, and a is the scaling parameter.
    Depending on the case considered some of these matrices can be null.
    
    Parameters
    ----------
    FE : dolfinx.fem.function.FunctionSpace
        Finite element approximation space for (ur,ui)
    dmesh : fenics_utils.DolfinMesh
    tags: dict(list(int))
        Tags that identifies the various domains.
        'gamma-d': Dirichlet boundary.
        'omega-m', 'omega-d':  Omega_m (Omega_d), outside of the PML region.
        -> Contributes to Am0 and Ad0.
        'omega-m-pml','omega-d-pml': Intersection of Omega_m (Omega_d) and 
        the PML region. This region must be meshed in Euler coordinates.
        -> Contributes to Am1, Am2, Am3, Am4 and Ad1, Ad2, Ad3, Ad4.
    per : dolfinx_mpc_utils.PeriodicBoundary
        Periodic boundary conditions (0 periodic boundary conditions possible)    
    diag_Am, diag_Ad : float, optional
        Diagonal penalty for enforcing Dirichlet boundary condition
    split_real_imag (bool)
        Split real and imaginary parts.

    Returns (If split_real_imag=False)
    -------
    A_m, A_d : dictionary of dolfin.cpp.la.PETScMatrix
        Dictionary of assembled matrices with entries:
            '1' = Am0, Ad0
            'a' = Am1, Ad1
            'ainv' = Am2, Ad2
        
    Returns (If split_real_imag=True)
    -------
    A_m, A_d : dictionary of dolfin.cpp.la.PETScMatrix
        Dictionary of assembled matrices with entries:
            '1' = Am0, Ad0
            'real-a' = Am1, Ad1
            'imag-a' = Am2, Ad2
            'real-ainv' = Am3, Ad3
            'imag-ainv' = Am4, Ad4
    idx_ur, idx_ui : list(int)
        Indices of real and imaginary parts.
    V_r, V_i: dolfin.function.functionspace.FunctionSpace
        Function spaces for real and imaginary components
        
    Remark
    ------
    The mesh should not contain overlapping physical entities. Otherwise,
    periodic boundary conditions break down in fenics.

    """
        # surface measure
    dx = dmesh.dx
        # -- Retrieve tags
    tags_Dirichlet = tags['gamma-d']
    tags_omegaM, tags_omegaD, tags_omegaM_PML, tags_omegaD_PML = [],[],[],[]
    if 'omega-m' in tags:
        tags_omegaM = tags['omega-m']
    if 'omega-d' in tags:
        tags_omegaD = tags['omega-d']
    if 'omega-m-pml' in tags:
        tags_omegaM_PML = tags['omega-m-pml']
    if 'omega-d-pml' in tags:
        tags_omegaD_PML = tags['omega-d-pml']
        # -- Dirichlet boundary condition
    uD = dolfinx.fem.Function(V)
    uD.vector.setArray(0)
    uD.vector.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)
    bcs = fenicsx_utils.create_DirichletBC(dmesh.mesh,dmesh.facet_tags,V,[uD],tags_Dirichlet)
        # -- Periodic boundary condition
    per.create_finalized_MultiPointConstraint(V, dmesh.facet_tags, bcs)
    mpc = per.get_MultiPointConstraint()
        # -- Definition of bilinear forms
    def a_grad(u,v,tags):
        a_list=[]
        for i in range(len(tags)):
                a_list.append(ufl.inner(ufl.grad(u), ufl.grad(v))*dx(tags[i]))
        return sum(a_list)
    
    def a_z(u,v,tags):
        a_list=[]
        for i in range(len(tags)):
                a_list.append(ufl.inner(ufl.Dx(u,0), ufl.Dx(v,0))*dx(tags[i]))
        return sum(a_list)        
    
    def a_theta(u,v,tags):
       a_list=[]
       for i in range(len(tags)):
               a_list.append(ufl.inner(ufl.Dx(u,1), ufl.Dx(v,1))*dx(tags[i]))
       return sum(a_list)

    if split_real_imag:
        (ur,ui) = ufl.TrialFunction(V)
        (vr,vi) = ufl.TestFunction(V)
    else:
        u = ufl.TrialFunction(V)
        v = ufl.TestFunction(V)
    a_m = dict(); a_d = dict();
    if (len(tags_omegaM)>0) and (len(tags_omegaD)>0):
        if split_real_imag:
                # x 1
            a_m['1'] = a_grad(ur,vr,tags_omegaM) + a_grad(ui,vi,tags_omegaM)
            a_d['1'] = a_grad(ur,vr,tags_omegaD) + a_grad(ui,vi,tags_omegaD)
        else:
            a_m['1'] = a_grad(u,v,tags_omegaM)
            a_d['1'] = a_grad(u,v,tags_omegaD)          
        
    if (len(tags_omegaM_PML)>0) and (len(tags_omegaD_PML)>0):
        if split_real_imag:
                # x (real(alpha))
            a_m['real-a'] = a_z(ur,vr,tags_omegaM_PML) + a_z(ui,vi,tags_omegaM_PML)
            a_d['real-a'] = a_z(ur,vr,tags_omegaD_PML) + a_z(ui,vi,tags_omegaD_PML)
                # x (imag(alpha))
            a_m['imag-a'] = -a_z(ui,vr,tags_omegaM_PML) + a_z(ur,vi,tags_omegaM_PML)
            a_d['imag-a'] = -a_z(ui,vr,tags_omegaD_PML) + a_z(ur,vi,tags_omegaD_PML)
                # x (real(1/alpha))
            a_m['real-ainv'] = a_theta(ur,vr,tags_omegaM_PML) + a_theta(ui,vi,tags_omegaM_PML)
            a_d['real-ainv'] = a_theta(ur,vr,tags_omegaD_PML) + a_theta(ui,vi,tags_omegaD_PML)        
                # x (imag(1/alpha))
            a_m['imag-ainv'] = -a_theta(ui,vr,tags_omegaM_PML) + a_theta(ur,vi,tags_omegaM_PML)
            a_d['imag-ainv'] = -a_theta(ui,vr,tags_omegaD_PML) + a_theta(ur,vi,tags_omegaD_PML)
        else:
                # x (alpha)
            a_m['a'] = a_z(u,v,tags_omegaM_PML)
            a_d['a'] = a_z(u,v,tags_omegaD_PML)
                # x (1/alpha)
            a_m['ainv'] = a_theta(u,v,tags_omegaM_PML)
            a_d['ainv'] = a_theta(u,v,tags_omegaD_PML)
        
        # -- Assemble sparse eigenvalue problem
    A_m = dict(); A_d =dict();
    for key in a_m:
        print(f"Assembly a_m ({key})")
        A_m[key] = dolfinx_mpc.assemble_matrix(dolfinx.fem.form(a_m[key]), mpc, bcs=bcs,diagval=diag_Am)
        print(f"Assembly a_d ({key})")
        A_d[key] = dolfinx_mpc.assemble_matrix(dolfinx.fem.form(a_d[key]), mpc, bcs=bcs,diagval=diag_Ad)

        # -- Indices of real and imaginary parts
        # Shortcut: N = A_m[0].size[0]; idx_ur = np.r_[0:N:2]; idx_ui = np.r_[1:N:2]        
    if split_real_imag:
        idx_ur = np.unique(V.sub(0).dofmap.list.array)
        idx_ui = np.unique(V.sub(1).dofmap.list.array)
        out = (A_m, A_d, idx_ur, idx_ui, V.sub(0), V.sub(1))
    else:
        out = (A_m,A_d)
    return out

class PEP_with_PML_fenicsx:
    """ Plasmonic Eigenvalue Problem (PEP) for a particle perturbed by N
    corners. Each corner is discretized in Euler coordinates (z,theta) and a
    PML region is defined using a scaling parameter alpha. """
    def __init__(self,dmesh,tags,pbc=dolfinx_mpc_utils.PeriodicBoundary(),
                                 split_real_imag=False,case_name="",
                                 mesh_name=""):
        self.dmesh = dmesh  # Dolfin mesh
        self.pbc = pbc      # Periodic boundary condition (dolfinx_mpc_utils.PeriodicBoundary)
        self.tags = tags    # Tags of physical zones (dict)
        self.case_name = case_name
        self.mesh_name = mesh_name
            # Split real and imaginary part to assemble and solve
        self.split_real_imag = split_real_imag
            # Exact eigenvalues
        self.k_even = None
        self.k_odd = None
            # Exact deformed critical interval
        self.k_Ic_even = None
        self.k_Ic_odd = None
            # Unperturbed eigenvalues
        self.k_unp_even = None
        self.k_unp_odd = None
            # Assembly
        self.A_m_l = None   # lhs of PEP
        self.A_d_l = None   # rhs of PEP
        self.V = None       # function space
        self.idx_ur = None  # indices of real part (if split_real_imag==True)
        self.idx_ui = None  # indices of imaginary part (if split_real_imag==True)
        self.V_r = None     # subspace for real part (if split_real_imag==True)
        self.V_i = None     # subspace for imaginary part (if split_real_imag==True)
            # Eigenvalue problem
        self.alpha = None   # scaling parameter
        self.eigval = None
        self.eigvec_r = None # list of numpy array
        self.eigvec_i = None # list of numpy array
            
    def assemble(self,V_fem):
        """Assemble PEP.
        
        V_fem: (string, int)
            Name of finite element, degree 
        """
            # -- Function space
        if self.split_real_imag: # Solution on (ur,ui)
            self.V = dolfinx.fem.VectorFunctionSpace(self.dmesh.mesh, V_fem,2)
            print("Assembly on (ur,ui) (real/imag split)....")
        else: # Solution on u
            self.V = dolfinx.fem.FunctionSpace(self.dmesh.mesh, V_fem)
            print("Assembly on (u) (no real/imag split)....")
            
            # -- Assembly
        out = assemble_PEP_PML_fenicsx(self.V, self.dmesh, self.tags, self.pbc,
                                       diag_Am=1e2, diag_Ad=1e-2,
                                       split_real_imag=self.split_real_imag)
        self.A_m_l = out[0]; self.A_d_l = out[1]
        self.V_r = self.V
            # -- (If splitting) Extract Re/Im indices
        if self.split_real_imag:
            self.idx_ur = out[2]; self.idx_ui = out[3]
            self.V_r = out[4].collapse()[0]; self.V_i = out[5].collapse()[0]            
        print("Done.")
        
    def solve(self,SLEPc_params,alpha=1.0,tolReIm=1e-5):
        self.alpha = alpha; A_m_l, A_d_l = self.A_m_l, self.A_d_l
        coeff = dict()
        coeff = {'1':1.0, 'a': alpha, 'ainv': 1/alpha,
                 'real-a': np.real(alpha), 'imag-a': np.imag(alpha), 
                 'real-ainv': np.real(1/alpha), 'imag-ainv': np.imag(1/alpha)}
        A_m, A_d = 0, 0
        
        if ('a' in A_m_l.keys()) and (isinstance(alpha,complex)) and (PETSc.ComplexType!=PETSc.ScalarType):
            raise ValueError("PETSc.ScalarType is real. Set split_real_imag=True.")
            
        for key in A_m_l:
            A_m += coeff[key]*A_m_l[key]
            A_d += coeff[key]*A_d_l[key]
              # Solve PEP
        EPS = SLEPc_utils.solve_GEP_shiftinvert(A_m,A_d,**SLEPc_params)
            # Extract/Sort eigenvalues and eigenvectors
        if self.split_real_imag:
            idx_ur, idx_ui = self.idx_ur, self.idx_ui
            (eigval,eigvec)=SLEPc_utils.EPS_get_spectrum_ReImFormulation(EPS,
                                            idx_ur.astype(np.int32),
                                            idx_ui.astype(np.int32),
                                            tol=tolReIm)
            self.pbc.set_slave_dofs(eigvec)
            self.eigval = eigval
            self.eigvec_r = list(); self.eigvec_i = list()
            for vec in eigvec:
                self.eigvec_r.append(vec[idx_ur]); self.eigvec_i.append(vec[idx_ui])
        else:
            (eigval,eigvec_r,eigvec_i) = SLEPc_utils.EPS_get_spectrum(EPS)
            self.pbc.set_slave_dofs(eigvec_r)
            self.pbc.set_slave_dofs(eigvec_i)
            self.eigval = eigval
            self.eigvec_r = list(); self.eigvec_i = list()
            for i in range(len(eigval)):
                self.eigvec_r.append(eigvec_r[i].array)
                self.eigvec_i.append(eigvec_i[i].array)
    
    def set_deformed_critical_interval(self,phi_list):
        """ Define deformed critical interval. """
        self.k_Ic_even, self.k_Ic_odd = list(), list()
        for phi in phi_list:
            psi = lambda eta : np.sinh(eta*np.pi)/np.sinh(eta*(np.pi-phi));
            self.k_Ic_even.append(lambda eta,alpha : (1-psi(eta*alpha))/(1+psi(eta*alpha)))
            self.k_Ic_odd.append(lambda eta,alpha : (1+psi(eta*alpha))/(1-psi(eta*alpha)))

    def set_exact_eigenvalues(self,k_even,k_odd):
        """ Define maps k_even and k_odd: one input argument N."""
        self.k_even = k_even
        self.k_odd = k_odd

    def set_unperturbed_eigenvalues(self,k_even,k_odd):
        """ Define maps k_even and k_odd: one input argument N."""
        self.k_unp_even = k_even
        self.k_unp_odd = k_odd

    def plot_spectrum(self,f,eta=0,N=10,plot_mesh_name=False):
        alpha = self.alpha
        ax=f.add_subplot(1,1,1)
        Ic_linestyle = ['-','--','-.',':']
        if self.k_Ic_even is not None:
            for i in range(len(self.k_Ic_even)):
                if len(self.k_Ic_even)==1: # do not number interval if there is only one corner
                    label = r'$I_{c,\alpha}^{\rm{even}}$'
                else:
                    label = r'$I_{c,\alpha}^{\rm{even}}$'+f'(#{i+1})'
                ax.plot(np.real(self.k_Ic_even[i](eta,alpha)),np.imag(self.k_Ic_even[i](eta,alpha)),label=label,linestyle=Ic_linestyle[np.mod(i,len(Ic_linestyle))],marker='None',mfc='none',color='b')
        if self.k_Ic_odd is not None:
            for i in range(len(self.k_Ic_odd)):
                if len(self.k_Ic_odd)==1: # do not number interval if there is only one corner
                    label = r'$I_{c,\alpha}^{\rm{odd}}$'
                else:
                    label = r'$I_{c,\alpha}^{\rm{odd}}$'+f'(#{i+1})'
                ax.plot(np.real(self.k_Ic_odd[i](eta,alpha)),np.imag(self.k_Ic_odd[i](eta,alpha)),label=label,linestyle=Ic_linestyle[np.mod(i,len(Ic_linestyle))],marker='None',mfc='none',color='r')
        if (self.k_even is not None) and (len(self.k_even(N))!=0):
            ax.plot(np.real(self.k_even(N)),np.imag(self.k_even(N)),label=r'$\kappa_{\rm{exact}}^{\rm{even}}$',linestyle='none',marker='o',mfc='none',color='b')
        if (self.k_odd is not None) and (len(self.k_odd(N))!=0):
            ax.plot(np.real(self.k_odd(N)),np.imag(self.k_odd(N)),label=r'$\kappa_{\rm{exact}}^{\rm{odd}}$',linestyle='none',marker='o',mfc='none',color='r')
        if (self.k_unp_even is not None) and (len(self.k_unp_even(N))!=0):
            ax.plot(np.real(self.k_unp_even(N)),np.imag(self.k_unp_even(N)),label=r'$\kappa_{\rm{unp}}^{\rm{even}}$',linestyle='none',marker='s',mfc='none',color='b')
        if (self.k_unp_odd is not None) and (len(self.k_unp_odd(N))!=0):
            ax.plot(np.real(self.k_unp_odd(N)),np.imag(self.k_unp_odd(N)),label=r'$\kappa_{\rm{unp}}^{\rm{odd}}$',linestyle='none',marker='s',mfc='none',color='r')
        ax.plot(-np.real(self.eigval),-np.imag(self.eigval),label=r'FEM ('+f'{self.eigvec_r[0].size:d} DoF)',linestyle='none',marker='x',mfc='none',color='k')        
        ax.set_xlabel(r"$\Re(\kappa)$")
        ax.set_ylabel(r"$\Im(\kappa)$")
        title = r"Spectrum with $\arg(\alpha)=$"+f"{np.rad2deg(np.angle(alpha)):3.4g}°"+f"\n{self.case_name}"
        if plot_mesh_name:
            title += f" <{self.mesh_name}>"            
        ax.set_title(title)
        ax.set_xlim([-3,0])
        ax.set_ylim([-1,1])
        ax.legend(bbox_to_anchor=[0,-0.25],loc="upper left",ncol=4)
        ax.set_aspect('equal')
        return ax

    def plot_eigenfunction_pyvista(self,ip):
        import pyvista as pv
        # pv.set_plot_theme("document")
        colormap = 'viridis'
        colorbar_args = dict(vertical=True,height=0.8,position_x=0.9, position_y=0.1,
                             fmt="%.1g",title_font_size=1,label_font_size=20)
        i_plot = np.min([ip,len(self.eigval)-1])               
        U_r = np.real(self.eigvec_r[i_plot])
        U_i = np.real(self.eigvec_i[i_plot])
        try: # Approach w/o interpolation (fast but works for isoparametric element only)
            U_n = np.abs(U_r+1j*U_i); 
            U_r, U_i, U_n = U_r/np.max(np.abs(U_n)), U_i/np.max(np.abs(U_n)), U_n/np.max(np.abs(U_n))
            grid = fenicsx_utils.create_pyvista_UnstructuredGrid_from_mesh(self.dmesh.mesh)
            grid.clear_data()
            grid["ur"] = U_r
            grid["ui"] = U_i
            try_interpolation = False
        except ValueError:
            import warnings
            warnings.warn("No-interpolation approach failed. Trying interpolation...")
            grid = None
            try_interpolation = True
        if try_interpolation: # Interpolation to mesh coordinates
            u_topology, u_cell_types, u_geometry = dolfinx.plot.create_vtk_mesh(self.V_r)
            grid = pv.UnstructuredGrid(u_topology, u_cell_types, u_geometry)
            grid.point_data["ur"] = np.real(U_r)
            grid.point_data["ui"] = np.real(U_i)
        
        plotter = pv.Plotter(shape=(2, 1))
        plotter.subplot(0,0)
        if try_interpolation:
            plotter.add_mesh(grid,scalars="ur",show_edges=False, show_scalar_bar=True,scalar_bar_args=colorbar_args,cmap=colormap,clim=[-1,1])
        else:
            plotter.add_mesh(grid, scalars="ur",show_edges=False, show_scalar_bar=True,scalar_bar_args=colorbar_args,cmap=colormap,clim=[-1,1])
        plotter.view_xy()
        plotter.show_bounds(xlabel="x",ylabel="y",all_edges=True,minor_ticks=True,font_size=10)
        plotter.add_title(f"PML angle {np.rad2deg(np.angle(self.alpha)):3.4g} deg\n kappa={-self.eigval[i_plot]:1.2g} [{i_plot}/{len(self.eigval)}]\n Re(u)/|u|_inf",font_size=10)
        plotter.subplot(1,0)
        if try_interpolation:
            plotter.add_mesh(grid,scalars="ui",show_edges=False, show_scalar_bar=True,scalar_bar_args=colorbar_args,cmap=colormap,clim=[-1,1])
        else:
            plotter.add_mesh(grid.copy(), scalars="ui",show_edges=False, show_scalar_bar=True,scalar_bar_args=colorbar_args,cmap=colormap,clim=[-1,1])
        plotter.view_xy()
        plotter.show_bounds(xlabel="x",ylabel="y",all_edges=True,minor_ticks=True,font_size=10)
        plotter.add_title("\n Im(u)/|u|_inf",font_size=10)
        plotter.show()

    def plot_eigenfunction_pyvista_matplotlib(self,fig,ip):
        import pyvista as pv
        import pyvista_utils
        # pv.set_plot_theme("document")
        colormap = 'viridis'
        i_plot = np.min([ip,len(self.eigval)-1])
        U_r = np.real(self.eigvec_r[i_plot])
        U_i = np.real(self.eigvec_i[i_plot])
        U_n = np.abs(U_r+1j*U_i); 
        U_r, U_i = U_r/np.max(U_n), U_i/np.max(U_n)
        u_topology, u_cell_types, u_geometry = dolfinx.plot.create_vtk_mesh(self.V_r)
        grid = pv.UnstructuredGrid(u_topology, u_cell_types, u_geometry)
        grid.point_data["ur"] = np.real(U_r)
        plotter = pv.Plotter(shape=(1, 1),window_size=[1000, 1000])
        plotter.add_mesh(grid,scalars="ur",show_edges=False, show_scalar_bar=False,cmap=colormap)
        plotter.view_xy()
        img = pyvista_utils.get_trimmed_screenshot(plotter)
        x_max=np.max(self.dmesh.mesh.geometry.x,0); x_min=np.min(self.dmesh.mesh.geometry.x,0)
        ax = fig.add_subplot(1,1,1)
        coll = ax.imshow(img,vmin=-1,vmax=1,extent=(x_min[0],x_max[0],x_min[1],x_max[1]),cmap=colormap)
        fig.colorbar(coll,ax=ax)
        ax.set_title(f"$ Re(u)/|u|_inf, \kappa$={-self.eigval[i_plot]:4.4g} DoF={self.dmesh.mesh.geometry.x.shape[0]}")
        ax.set_xlabel('$x$')
        ax.set_ylabel('$y$')
        ax.grid(False)
        ax.set_aspect('equal')
        return ax
    
    def export_eigenfunction_pyvista(self,fname):
        grid = fenicsx_utils.create_pyvista_UnstructuredGrid_from_mesh(self.dmesh.mesh)
        grid.clear_arrays()
        for i in range(len(self.eigval)):
            name = f"k_{-self.eigval[i]:4.3g}"

            U_r = np.real(self.eigvec_r[i])
            U_i = np.real(self.eigvec_i[i])
            try: # Approach w/o interpolation (fast but works for isoparametric element only)
                U_n = np.abs(U_r+1j*U_i); 
                U_r, U_i = U_r/np.max(np.abs(U_n)), U_i/np.max(np.abs(U_n))
                grid[name+"_ur"] = U_r
                grid[name+"_ui"] = U_i
                try_interpolation = False
            except ValueError:
                import warnings
                warnings.warn("No-interpolation approach failed. Trying interpolation...")
                grid = None
                try_interpolation = True
            if try_interpolation: # Interpolation to mesh coordinates
                import pyvista as pv
                u_topology, u_cell_types, u_geometry = dolfinx.plot.create_vtk_mesh(self.V_r)
                grid = pv.UnstructuredGrid(u_topology, u_cell_types, u_geometry)
                grid.point_data[name+"_ur"] = np.real(U_r)
                grid.point_data[name+"_ui"] = np.real(U_i)
 
        fname = fname+".vtu"
        grid.save(fname)
        print(f"Eigenvectors saved to {fname}.")

    def export_eigenvalues_to_csv(self,fname):        
        ar = np.array(self.eigval)
        ar = np.hstack((np.real(ar[:,None]), np.imag(ar[:,None])))
        np.savetxt("FEM-"+fname,
                   ar,
                   header='Re,Im', comments='# Eigenvalues PEP\n',
                   fmt='%.6e', delimiter=',', newline='\n')
        print(f"Eigenvalues saved to {fname}.")        

# Definition of periodic boundary conditions

def get_pbc_corner_fenicsx(tag_bot,y_offset=0,atol=1e-14,rtol=1e-11):
    """ 
    Get periodic boundary conditions for a corner in Euler coordinates (z,theta).

    Parameters
    ----------
    y_offset : scalar
        y-offset of the computational domain.
        
    Returns
    -------
    per : dolfinx_mpc_utils.PeriodicBoundary()
        Periodic boundary conditions    
    """

    def EulerBot_map(x): # Bottom boundary {y-y_offset=-pi}
        return np.isclose(x[1],y_offset-np.pi,rtol=rtol,atol=atol)
    def EulerTop_map(x): # Top boundary {y-y_offset=pi}
        return np.isclose(x[1],y_offset+np.pi,rtol=rtol,atol=atol)
    def map_EulerBot_to_EulerTop(x):
        y = np.zeros(x.shape)
        y[0] = x[0]; y[1] = x[1]+2*np.pi; y[2] = x[2]
        return y
    pbc = dolfinx_mpc_utils.PeriodicBoundary()
    pbc.add_topological_condition(tag_bot, map_EulerBot_to_EulerTop,
                                  slave_map=EulerBot_map,
                                  master_map=EulerTop_map)
    return pbc

def get_pbc_ellipseNcorners_fenicsx(a_m,b_m,phi,pos,
                                    x_offset,y_offset,
                                    tag_EulerBot_bnd,tag_EulerRight_bnd,
                                    atol=1e-14,rtol=1e-11):
    """ Build dolfinx_mpc_utils.PeriodicBoundary for an ellipse perturbed by 
    N corners, each meshed in Euler coordinates.
    
    Inputs
    ------
        a_m, b_m: float
            Semi-axes of the elliptical particle.
        phi: list(float)
            Angles of the N corners.
        pos: list(string)
            Position of the N corners.
        x_offset, y_offset: list(float)
            Offsets for the N Euler regions.
        tag_EulerBot_bnd, tag_EulerRight_bnd: list(list(int))
            Tags of the bottom and right boundaries of the N Euler regions.
    """
    pbc = dolfinx_mpc_utils.PeriodicBoundary()
    
    # bottom of Euler region (slave) = top of Euler region (master)
    def EulerBot_map(x,y_offset):
        return np.isclose(x[1],y_offset-np.pi,rtol=rtol,atol=atol)
    def EulerTop_map(x,y_offset):
        return np.isclose(x[1],y_offset+np.pi,rtol=rtol,atol=atol)
    def map_EulerBot_to_EulerTop(x):
        y = np.zeros(x.shape)
        y[0] = x[0]; y[1] = x[1]+2*np.pi; y[2] = x[2]
        return y    
    for i, tags in enumerate(tag_EulerBot_bnd):
        pbc.add_topological_condition(tags, map_EulerBot_to_EulerTop,
            slave_map=lambda x, y_offset=y_offset[i]: EulerBot_map(x,y_offset),
            master_map=lambda x, y_offset=y_offset[i]: EulerTop_map(x,y_offset))
    
    # right of Euler region (slave) = corner circle (master)
    def map_EulerRight_to_LeftCircle(x,x_c,y_c,R,x_offset,y_offset):
        y = np.zeros(x.shape)
        # print(f"Radius on corner right boundary: R={np.exp(x[0]-x_offset)}")
        y[0] = x_c+R*np.cos(x[1]-y_offset); y[1] = y_c+R*np.sin(x[1]-y_offset); y[2]=x[2]
        return y
    def map_EulerRight_to_TopCircle(x,x_c,y_c,R,x_offset,y_offset):
        y = np.zeros(x.shape)
        # print(f"Radius on corner right boundary: R={np.exp(x[0]-x_offset)}")
        y[0] = x_c+R*np.sin(x[1]-y_offset); y[1] = y_c-R*np.cos(x[1]-y_offset); y[2]=x[2]
        return y
    for i, tags in enumerate(tag_EulerRight_bnd):
        (x_c,y_c,x_m,y_m,R) = PEP_utils_ana.get_C1corner_ellipse(a_m,b_m,phi[i],pos=pos[i])
        if (pos[i] == "left"):
            slave_to_master_map = lambda x,x_c=x_c,y_c=y_c,R=R,x_offset=x_offset[i],y_offset=y_offset[i]: \
                                    map_EulerRight_to_LeftCircle(x,x_c,y_c,R,x_offset,y_offset)
        elif (pos[i] == "top"):
            slave_to_master_map = lambda x,x_c=x_c,y_c=y_c,R=R,x_offset=x_offset[i],y_offset=y_offset[i]: \
                                    map_EulerRight_to_TopCircle(x,x_c,y_c,R,x_offset,y_offset)
        else:
            raise ValueError("Corner position not implemented.")
        pbc.add_topological_condition(tags, slave_to_master_map)        
    return pbc