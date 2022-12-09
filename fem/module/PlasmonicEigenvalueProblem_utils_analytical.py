#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analytical functions for solving the plasmonic eigenvalue problem (PEP).

"""

import numpy as np
import os

    # Module directory
DIR_MODULE = os.path.dirname(__file__)


def ellipse_eigenvalues(a_m,b_m,a_d,b_d,N):
    """
    Exact eigenvalues for an elliptical particle (Omega_m) embedded in a confocal
    elliptic domain (Omega_d).

    Parameters
    ----------
    a_m : float
        Major semi-axis of Omega_m.
    b_m : float
        Minor semi-axis of Omega_m.
    a_d : float
        Major semi-axis of Omega_d.
    b_d : float
        Minor semi-axis of Omega_d.
    N : int
        Number of eigenvalues.

    Returns
    -------
    kappa_ex_ev : numpy.ndarray
        Eigenvalues associated with eigenfunctions even w.r.t the major axis.
    kappa_ex_odd : numpy.ndarray
        Eigenvalues associated with eigenfunctions odd w.r.t the major axis.

    Remark
    -------
    
    The formula used for the eigenvalues assumes that the ellipses Omega_m
    and Omega_d are confocal, see Prop. 31 of 10.1016/j.jcp.2021.110433.
    
    """    
    c = np.sqrt(a_m**2-b_m**2) # a_m > b_m since it is the major axis
    mu_m = np.arccosh(a_m/c)
    mu_d = np.arccosh(a_d/c)    
    n = np.r_[1:N]
    n_ev = np.concatenate(( [0], n )) # add 0 for even eigenvalues
    kappa_ex_ev = -np.tanh(n_ev*(mu_d-mu_m))*np.tanh(n_ev*mu_m)
    kappa_ex_odd = -np.tanh(n*(mu_d-mu_m))/np.tanh(n*mu_m)
    return (kappa_ex_ev,kappa_ex_odd)

def ellipse_parameters(case_name):
    """
    Return value of ellipse parameters (a_m,b_m,a_d,b_d) corresponding to
    case_name.
    """
    # JCP2021 = 10.1016/j.jcp.2021.110433
    if case_name == "JCP2021_figure_9":
        a_d =3; b_d = 2
        c = np.sqrt(a_d**2-b_d**2)
        a_m =1.10*c; b_m = np.sqrt(a_m**2-c**2)
    elif case_name == "JCP2021_figure_10":
        a_m =2.5; b_m = 1; c = np.sqrt(a_m**2-b_m**2)
        a_d =3; b_d = np.sqrt(a_d**2-c**2)        
    elif case_name == "JCP2021_figure_12":
        mu_m = 0.5; mu_d = 1.25;
        a_m = 2.5; b_m = 1; c = np.sqrt(a_m**2-b_m**2) # arbitrary
        a_d = a_m * np.cosh(mu_d)/np.cosh(mu_m);
        b_d = np.sqrt(a_d**2-c**2)
    elif case_name == "JCP2021_case_A":
        (a_m,b_m,a_d,b_d) = perturbed_ellipse_parameters(case_name)[0:4]
    elif case_name == "JCP2021_case_B":
        (a_m,b_m,a_d,b_d) = perturbed_ellipse_parameters(case_name)[0:4]
    else:
        raise ValueError("Unknown case name.")
     
    return (a_m,b_m,a_d,b_d)

def perturbed_ellipse_parameters(case_name):
    """
    Return geometrical values describing an ellipse perturbed by corners.

    Parameters
    ----------
    case_name : str
        Reference name.

    Returns
    -------
    a_m, b_m : float
        Semi-axes of ellipse m.
    a_d, b_d : TYPE
        Semi-axes of ellipse d.
    x_c, y_c : list(float)
        Coordinates of each corner.
    x_m, y_m : list(float)
        Coordinates of junction point for each corner.
    phi : list(float)
        Angle for each corner.
    corner_pos : list(str)
        Position of each corner.
        
    """
        
    if case_name == "JCP2021_case_A":
        a_m =2.5; b_m = 1
        c = np.sqrt(a_m**2-b_m**2)
        phi = [0.75*np.pi]
        corner_pos = ["left"]
        (x_c,y_c,x_m,y_m,R) = get_C1corner_ellipse(a_m,b_m,phi[0],pos=corner_pos[0])
        a_d =np.abs(x_c) + 1.5*R; b_d = np.sqrt(a_d**2-c**2)
        x_c=[x_c]; y_c=[y_c]; x_m=[x_m]; y_m=[y_m]; R=[R]
    elif case_name == "JCP2021_case_B":
        a_m =2.5; b_m = 1
        c = np.sqrt(a_m**2-b_m**2)
        phi = [0.63*np.pi]
        corner_pos = ["left"]
        (x_c,y_c,x_m,y_m,R) = get_C1corner_ellipse(a_m,b_m,phi[0],pos=corner_pos[0])
        a_d =np.abs(x_c) + 1.5*R; b_d = np.sqrt(a_d**2-c**2)
        x_c=[x_c]; y_c=[y_c]; x_m=[x_m]; y_m=[y_m]; R=[R]

    return (a_m,b_m,a_d,b_d,x_c,y_c,x_m,y_m,phi,corner_pos)


def get_C1corner_ellipse(a,b,phi,pos="left"):
    """
    Compute geometrical parameters of the C^1 corner perturbation with angle phi
    in (0,pi), for an ellipse of semi-axes (a,b).

    Parameters
    ----------
    a, b : float
        Ellipse semi-axes.
    phi : float
        Corner angle in (0,pi).
    pos : string, optional
        Corner position ("left","top","right").

    """
     
    def ellipse_get_junction_point(a,b,phi):
       """ Return point (x,y) with slope tan(phi/2) with y>0. """
       x_j = np.sqrt(np.tan(phi/2)**2+(b/a)**2);
       x_j = -a * np.tan(phi/2) / x_j;
       y_j = b*np.sqrt(1-(x_j/a)**2);
       return (x_j,y_j)
           
    if (pos=="left"):
            # top junction point (y_j>0)
        (x_j,y_j) = ellipse_get_junction_point(a,b,phi)
            # abscissa of top point of corner
        x_c = x_j - y_j/np.tan(phi/2);
        y_c = 0;
    elif (pos=="top"):
        psi = np.pi - phi
        (x_j,y_j) = ellipse_get_junction_point(a,b,psi)
            # abscissa of top point of corner
        x_c = 0
        y_c = y_j - x_j*np.tan(psi/2);        
    elif (pos=="right"):
            # top junction point (y_j>0)
        (x_j,y_j) = ellipse_get_junction_point(a,b,-phi)
            # abscissa of top point of corner
        x_c = x_j + y_j/np.tan(phi/2);
        y_c = 0      
        # radius of corner circle
    R = np.sqrt((x_c-x_j)**2+(y_c-y_j)**2); 
    return (x_c,y_c,x_j,y_j,R)