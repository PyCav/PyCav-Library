import numpy as np
from scipy.integrate import quad
from copy import copy

import matplotlib.pyplot as plt

def numerov(x,dx,V,E,initial_values,params):
    steps = x.shape[0]-1
    psi = np.zeros(steps+1)
    psi[:2] = initial_values
    
    def g_n(n):
        return 2.*params[0]*(E-V(x[n]))
    
    def f_n(n):
        return (1.+g_n(n)*dx**2/12.)
    
    def inv_phi(phi,n):
        return phi/f_n(n)
    
    for i in range(1,steps):
        phi_2 = (12.-10.*f_n(i))*psi[i]-f_n(i-1)*psi[i-1]
        psi[i+1] = inv_phi(phi_2,i+1)
        
    return psi

def bisection_search(x,dx,V,params,bracket_E,tolerance = 0.5,max_evals = 1000):
    N = 0
    boundary = False
    def _boundary_check(val):
        if abs(val) < tolerance:
            return True
        else:
            return False
    
    order_bracket = np.sign(numerov(x,dx,V,bracket_E[0],[0.,dx],params)[-1])
    
    while not boundary:
        N += 1
        E_bisect = (bracket_E[0]+bracket_E[1])/2.0
        psi = numerov(x,dx,V,E_bisect,[0.,dx],params)
        boundary = _boundary_check(psi[-1])
        if not boundary:
            if order_bracket*psi[-1] < 0:
                bracket_E[1] = E_bisect
            else:
                bracket_E[0] = E_bisect
        if N > max_evals:
            print('Search reached max_evals, change bracketed values!')
            boundary = True
                
    return psi, E_bisect

def plot_levels(unperturb_wf,unperturb_erg,potential,params,x_lims = [-0.1,0.1], plot_zoom = [200]):

    x = np.linspace(x_lims[0],x_lims[1],500)

    fig_HO = plt.figure(figsize=(12,8))
    ax = plt.subplot(111)
    ax.set_xlim((np.min(x),np.max(x)))

    V = potential(params)
    ax.plot(x,V(x),'k')

    for i in range(4):
        psi_n = unperturb_wf(params,i)
        ax.plot(x,plot_zoom[0]*psi_n(x)+unperturb_erg(params,i))

    ax.set_ylim((0,unperturb_erg(params,4)))

    ax.set_yticks([])
    ax.set_xticks([])

    ax.set_ylabel('Energy')
    ax.set_xlabel('x')

def first_order_energy_sft(n,H,unperturb_wf,params,limits = [-np.inf,np.inf]):
    
    def combine_functions(m):
        def integrand(x):
            psi_m = unperturb_wf(params,m)
            return H(x)*psi_m(x)*psi_m(x)
        return integrand
    
    if hasattr(n, "__len__"):
        E = []
        for i in n:
                E_n, E_errn = quad(combine_functions(i),limits[0],limits[1])
                E.append(E_n)
        return E
    else:
        E, E_err = quad(combine_functions(n),limits[0],limits[1])
        return E

def first_order_wf(n,H,unperturb_wf,unperturb_erg,params,tolerance = 0.01, limits = [-np.inf,np.inf], return_list = False):
    
    def combine_functions(i,j):
        def integrand(x):
            psi_i = unperturb_wf(params,i)
            psi_j = unperturb_wf(params,j)
            return H(x)*psi_i(x)*psi_j(x)
        return integrand
      
    def E_diff(i,j):
        E_i = unperturb_erg(params,i)
        E_j = unperturb_erg(params,j)
        return E_i-E_j

    def find_k_values(even,k_list,I_list):
        truncate = False
        if even:
            k = 0
        else:
            k = 1
        while not truncate:
            if n != k:
                I, I_err = quad(combine_functions(n,k),limits[0],limits[1])
                E_nk = E_diff(n,k)
        
                if abs(I/E_nk) < tolerance:
                    truncate = True
                else:
                    I_list.append(I/E_nk)
                    k_list.append(k)
                    k += 2
            else:
                k += 2
        return k_list,I_list
    
    k_list = []
    I_list = []
    
    k_list,I_list = find_k_values(True,k_list,I_list)
    k_list,I_list = find_k_values(False,k_list,I_list)
       
    def sum_functions(n):
        def perturbed_wf(x):
            perturb_total = 0.0
            perturb_wf_n = unperturb_wf(params,n)
            perturb_total = perturb_wf_n(x)
            for k_sig in k_list:
                perturb_wf_k = unperturb_wf(params,k_sig)
                perturb_total = perturb_total + I_list[k_list.index(k_sig)]*perturb_wf_k(x)
            return perturb_total
        return perturbed_wf

    if not return_list:
        return sum_functions(n)
    else:
        return sum_functions(n),copy(k_list),copy(I_list)

def plot_perturb(H,unperturb_wf,unperturb_erg,potential,params,x_lims = [-0.1,0.1], tolerance = 0.01, limits = [-np.inf,np.inf], plot_zoom = [200,100]):
    
    colour_dict = {0 : 'r', 1 : 'b', 2 : 'g', 3 : 'c'}
    
    x = np.linspace(x_lims[0],x_lims[1],500)

    fig_HO = plt.figure(figsize=(12,8))
    ax = plt.subplot(111)
    ax.set_xlim((np.min(x),np.max(x)))
    
    perturb_h = [H(y) for y in x]
    
    V = potential(params)
    ax.plot(x,V(x)+perturb_h,'k')
    ax.set_ylim((0,unperturb_erg(params,4)))

    for i in range(4):
        psi_n = first_order_wf(i,H,unperturb_wf,unperturb_erg,params,tolerance = tolerance, limits = limits)
        ax.plot(x,plot_zoom[0]*psi_n(x)+unperturb_erg(params,i)+plot_zoom[1]*first_order_energy_sft(i,H,unperturb_wf,params,limits = limits), c = colour_dict[i])

    for i in range(4):
        psi_n = unperturb_wf(params,i)
        ax.plot(x,plot_zoom[0]*psi_n(x)+unperturb_erg(params,i), c = colour_dict[i], alpha = 0.5)

    ax.set_yticks([])
    ax.set_xticks([])

    ax.set_ylabel('Energy')
    ax.set_xlabel('x')
