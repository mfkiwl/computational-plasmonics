#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Design of two-corner perturbations.

This script plots various geometries perturbed by corner. It is useful to find
the geometrical parameters that will lead to embedded eigenvalues and 
complex resonances.

"""
    # Matrix
import numpy as np
import matplotlib.pyplot as plt
from numpy import pi
import PlasmonicEigenvalueProblem_utils as PEP_utils

def ellipse_exact_eigenvalues(N,mu_m,mu_d):
    """ Exact eigenvalues for PEP. mu_m, mu_d: elliptical coordinates. """
    n = np.r_[1:N]
    n_ev = np.concatenate(( [0], n )) # add 0 for even eigenvalues
    kappa_ex_ev = -np.tanh(n_ev*(mu_d-mu_m))*np.tanh(n_ev*mu_m)
    kappa_ex_odd = -np.tanh(n*(mu_d-mu_m))/np.tanh(n*mu_m)
    return (kappa_ex_ev,kappa_ex_odd)

#%% Elliptical particle
# Corner on both major and minor axis
plot_geometry = True
plot_spectrum = True

#a_d=1.1*a_m; b_d=np.sqrt(a_d**2-c**2)

    # JCP case A
a_m=2.5; b_m=1; c = np.sqrt(a_m**2-b_m**2)
a_d = 2.799469178605668; b_d = 1.608423974567369
phi1 = 0.75*pi # major axis, angle in (0,pi)
phi2 = 0.8*pi # minor axis, angle in (0,pi)
    # JCP case B
a_m=2.5; b_m=1; c = np.sqrt(a_m**2-b_m**2)
phi1 = 0.75*pi # major axis, angle in (0,pi)
phi2 = 0.88*pi # minor axis, angle in (0,pi)
(x_c,y_c,x_j,y_j,R) = PEP_utils.get_C1corner_ellipse(a_m,b_m,phi1)
a_d = a_m + 3.5*R; b_d = np.sqrt(a_d**2-c**2)
    # Quasi-circular case, with large perturbation
a_m =2.5; b_m = 2.45; c = np.sqrt(a_m**2-b_m**2)
phi1 = 0.65*np.pi
phi2 = 0.65*np.pi # minor axis, angle in (0,pi)
(x_c,y_c,x_j,y_j,R) = PEP_utils.get_C1corner_ellipse(a_m,b_m,phi1)
a_d =np.abs(x_c) + 1.5*R; b_d = np.sqrt(a_d**2-c**2)
    # Another quasi-circular case, with smaller perturbation
a_m =2.5; b_m = 2.45; c = np.sqrt(a_m**2-b_m**2)
phi1 = 0.75*np.pi
phi2 = 0.75*np.pi # minor axis, angle in (0,pi)
(x_c,y_c,x_j,y_j,R) = PEP_utils.get_C1corner_ellipse(a_m,b_m,phi1)
a_d =np.abs(x_c) + 4*R; b_d = np.sqrt(a_d**2-c**2)

    # Normalized radius of region mapped in Euler coordinates in (0,1)
R_PML1 = 1
R_PML2 = 1

if plot_geometry==True:
    f = plt.figure()
    ax=f.add_subplot(1,1,1)
    t = np.linspace(0,2*pi,int(1e3))
    
    ax.plot(a_m*np.cos(t), b_m*np.sin(t),color='k')
    ax.plot(a_d*np.cos(t), b_d*np.sin(t),color='k')
    ax.plot(c,0,marker='.',color='k')
    ax.plot(-c,0,marker='.',color='k')
    
    (x_c,y_c,x_j,y_j,R) = PEP_utils.get_C1corner_ellipse(a_m,b_m,phi1)
    ax.plot([x_j,x_c,x_j],[y_j,y_c,-y_j],marker='+',color='b')
    ax.plot(x_c+R_PML1*R*np.cos(t),y_c+R_PML1*R*np.sin(t),color='b',linestyle='--')
    R1 = R
    (x_c,y_c,x_j,y_j,R) = PEP_utils.get_C1corner_ellipse(a_m,b_m,phi2,pos="top")    
    ax.plot([x_j,x_c,-x_j],[y_j,y_c,y_j],marker='+',color='r')
    ax.plot(x_c+R_PML2*R*np.cos(t),y_c+R_PML2*R*np.sin(t),color='r',linestyle='--')
    ax.plot(x_c+R1*np.cos(t),y_c+R1*np.sin(t),color='b',linestyle='--')
    
    title = "Elliptical particle perturbed by two corners\n"+f"$(a_m,b_m)$=({a_m:2.2g},{b_m:2.2g}), $(a_d,b_d)$=({a_d:2.2g},{b_d:2.2g})\n"+r"$(\phi_{major}/2,\phi_{minor}/2)$"+f"=({np.rad2deg(phi1/2):3g}°,{np.rad2deg(phi2/2):3g}°)"
    ax.set_title(title)
    ax.set_ylim([-1.1*b_d,1.1*b_d])
    ax.set_xlim([-1.1*a_d,1.1*a_d])
    ax.set_aspect('equal') # orthonormal axis

if plot_spectrum == True:
         # Even major / even minor
    kappa_ex_ev_ev_fun = lambda N : ellipse_exact_eigenvalues(N, np.arccosh(a_m/c), np.arccosh(a_d/c))[0][0::2]
         # Even major / odd minor
    kappa_ex_ev_odd_fun = lambda N : ellipse_exact_eigenvalues(N, np.arccosh(a_m/c), np.arccosh(a_d/c))[0][1::2]
        # Odd major / even minor    
    kappa_ex_odd_ev_fun = lambda N : ellipse_exact_eigenvalues(N, np.arccosh(a_m/c), np.arccosh(a_d/c))[1][0::2]
        # Odd major / odd minor
    kappa_ex_odd_odd_fun = lambda N : ellipse_exact_eigenvalues(N, np.arccosh(a_m/c), np.arccosh(a_d/c))[1][1::2]
        
    f = plt.figure()
    ax=f.add_subplot(1,1,1)
    ax.plot(kappa_ex_ev_ev_fun(10),0*kappa_ex_ev_ev_fun(10),marker="o",linestyle='none',label='Even/Even',color='blue',fillstyle='none')
    ax.plot(kappa_ex_ev_odd_fun(10),0*kappa_ex_ev_odd_fun(10),marker="o",linestyle='none',label='Even/Odd',color='blue')
    ax.plot(kappa_ex_odd_odd_fun(10),0*kappa_ex_odd_odd_fun(10),marker="s",linestyle='none',label='Odd/Odd',color='red',fillstyle='none')
    ax.plot(kappa_ex_odd_ev_fun(10),0*kappa_ex_odd_ev_fun(10),marker="s",linestyle='none',label='Odd/Even',color='red')
    
    ax.set_title(title)
    
    PEP = PEP_utils.PEP_with_PML_fenicsx([],[],[])
    eta = np.linspace(1e-16,8,200); # in (0,inf)
        # Corner on major axis
    PEP.set_deformed_critical_interval([phi1])
    k = PEP.k_Ic_even[0](eta,1.0)
    ax.plot(np.real(k),-0.01+np.imag(k),label=r'$I_c^e$ (major)',color='b',linestyle='--')
    k = PEP.k_Ic_odd[0](eta,1.0)
    ax.plot(np.real(k),-0.01+np.imag(k),label=r'$I_c^o$ (major)',color='r',linestyle='--')
        # Corner on major axis
    PEP.set_deformed_critical_interval([phi2])
    k = PEP.k_Ic_even[0](eta,1.0)
    ax.plot(np.real(k),0.01+np.imag(k),label=r'$I_c^e$ (minor)',color='b')
    k = PEP.k_Ic_odd[0](eta,1.0)
    ax.plot(np.real(k),0.01+np.imag(k),label=r'$I_c^o$ (minor)',color='r')
    ax.legend()
    ax.set_xlim([-1,-0])

#%% Elliptical annulus
    # Outer sign-change interface
a_m=2.5; b_m=1; c = np.sqrt(a_m**2-b_m**2)

    # Inner sign-changing interface
a_mi=1.03*c;
# a_mi=0.9*a_m; 
b_mi=np.sqrt(a_mi**2-c**2); # a_mi > c
# a_mi=0.5*a_m; 
# b_mi=0.6*b_m;
    # Dirichlet ellipse (outer)
a_d=1.4*a_m; b_d=np.sqrt(a_d**2-c**2)
    # Angle on outer sign-changing interface (0,pi)
phi1 = pi/2
R_PML1 = 0.5 # factor in (0,1)
    # Angle on inner sign-changing interface (pi,2*pi)
phi2 = 2*pi-phi1
R_PML2 = 0.5 # factor in (0,1)

f = plt.figure()
ax=f.add_subplot(1,1,1)

t = np.linspace(0,2*pi,int(1e3))

ax.plot(a_m*np.cos(t), b_m*np.sin(t),color='k')
ax.plot(a_mi*np.cos(t), b_mi*np.sin(t),color='k')
ax.plot(c,0,marker='.',color='k')
ax.plot(-c,0,marker='.',color='k')

(x_c,y_c,x_j,y_j,R) = PEP_utils.get_C1corner_ellipse(a_m,b_m,phi1)
ax.plot([x_j,x_c,x_j],[y_j,y_c,-y_j],marker='+',color='b')
ax.plot(x_c+R_PML1*R*np.cos(t),y_c+R_PML1*R*np.sin(t),color='b')

(x_c,y_c,x_j,y_j,R) = PEP_utils.get_C1corner_ellipse(a_m,b_m,2*pi-phi2,pos="right")
ax.plot([x_j,x_c,x_j],[y_j,y_c,-y_j],marker='+',color='r')
ax.plot(x_c+R_PML2*R*np.cos(t),y_c+R_PML2*R*np.sin(t),color='r')

(x_c,y_c,x_j,y_j,R) = PEP_utils.get_C1corner_ellipse(a_mi,b_mi,2*pi-phi2)
ax.plot([x_j,x_c,x_j],[y_j,y_c,-y_j],marker='+',color='r')
ax.plot(x_c+R_PML1*R*np.cos(t),y_c+R_PML1*R*np.sin(t),color='r')

title = f"Elliptical annulus\n$(a_m,b_m)$=({a_m:2.2g},{b_m:2.2g}),"+r" $(a_{mi},b_{mi})$="+f"({a_mi:2.2g},{b_mi:2.2g})\n $(\phi_1/2,\phi_2/2)$=({np.rad2deg(phi1/2):3g}°,{np.rad2deg(phi2/2):3g}°)"
ax.set_title(title)
ax.set_ylim([-b_d,b_d])
ax.set_xlim([-a_d,a_d])
# ax.set_ylim([-2*y_j,2*y_j])
# ax.set_xlim([-1.1*a_m,-0.9*a_mi])
ax.set_aspect('equal') # orthonormal axis