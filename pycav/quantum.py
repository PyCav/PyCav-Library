import numpy as np
from scipy.integrate import quad
from copy import copy
import matplotlib.pyplot as plt

h_bar = 1.0

def create_Su(s):
    """
    Return the 'up' angular momentum operator

    Parameters
    ----------
    s : float or int
        Angular momentum quantum number (half-integer or integer)
    """

    dim = int(2*s+1)
    Su = np.zeros((dim,dim), dtype = np.complex128)
    for i in range(dim-1):
        m = s-1-i
        Su[i,i+1] = h_bar*np.sqrt(s*(s+1)-m*(m+1))
    return Su



def create_Sd(s):
    """
    Return the 'down' angular momentum operator

    Parameters
    ----------
    s : float or int
        Angular momentum quantum number (half-integer or integer)
    """

    dim = int(2*s+1)
    Sd = np.zeros((dim,dim), dtype = np.complex128)
    for i in range(1,dim):
        #magnetic quantum number
        m = s+1-i
        Sd[i,i-1] = h_bar*np.sqrt(s*(s+1)-m*(m-1))
    return Sd



def create_Sx(s):
    """
    Return the x angular momentum operator in the z basis

    Parameters
    ----------
    s : float or int
        Angular momentum quantum number (half-integer or integer)
    """

    return 0.5*(create_Su(s)+create_Sd(s))



def create_Sy(s):
    """
    Return the y angular momentum operator in the z basis

    Parameters
    ----------
    s : float or int
        Angular momentum quantum number (half-integer or integer)
    """

    return -0.5j*(create_Su(s)-create_Sd(s))



def create_Sz(s):
    """
    Return the z angular momentum operator in the z basis

    Parameters
    ----------
    s : float or int
        Angular momentum quantum number (half-integer or integer)
    """

    dim = int(2*s+1)
    Sz = np.zeros((dim,dim), dtype = np.complex128)
    for i in range(dim):
        Sz[i,i] = h_bar*(s-i)
    return Sz



class SpinSystem:
    """
    class representing a system of interacting spins in an inhomogeneous
    magnetic field. Interactions between spins are included in the 
    Hamiltonian as Kronecker products of the one-particle operators.
    """

    def __init__(self, spins, couplings, B_field = None):
        """
        Parameters
        ----------
        spins: list
            Spin quantum numbers of the particles (integer or half-integer)
        couplings: list
            NxN list for a system of N spins. The element J[i,j] gives the
            coupling strength between spins i and j. Positive coupling favours
            spin alignment.
        B_field: list, optional
            List of 3-dimensional vectors representing the magnetic flux 
            density at each spin. 
        """
        try:
            self.N = spins.shape[0]
        except AttributeError:
            self.N = len(spins)

        self.Sx = [create_Sx(s) for s in spins]
        self.Sy = [create_Sy(s) for s in spins]
        self.Sz = [create_Sz(s) for s in spins] 
        self.I = [np.identity(int(2*s+1), dtype = np.complex128) for s in spins]
        
        #interactions
        self.J = np.array(couplings, dtype = np.complex128)
        if B_field is not None:
            self.B = np.array(B_field, dtype = np.complex128)
        else:
            self.B = None
        
        #Hamiltonian
        self.H_size = 1
        for s in spins:
            self.H_size *= int(2*s+1)
        self.create_hamiltonian()

    def create_hamiltonian(self):
        """
        Create a Hamiltonian with spin-pairing. A linear coupling between
        spins and the magnetic field is included if a B_field was used to
        initialise the SpinSystem.
        The Hamiltonian acts on Kronecker products of the spins' states.
        """
        self.H = np.zeros((self.H_size, self.H_size), dtype = np.complex128)
        for i in range(self.N):
            #count from i+1 to N to avoid self-interaction and double counting
            for j in range(i+1, self.N):
                #create term of the form
                #I x S x I x I x S ...
                if i==0:
                    term_x = copy(self.Sx[0])
                    term_y = copy(self.Sy[0])
                    term_z = copy(self.Sz[0])
                else:
                    term_x = copy(self.I[0])
                    term_y = copy(self.I[0])
                    term_z = copy(self.I[0])
                for k in range(1, self.N):
                    if k == i or k == j:
                        term_x = np.kron(term_x, self.Sx[k])
                        term_y = np.kron(term_y, self.Sy[k])
                        term_z = np.kron(term_z, self.Sz[k])
                    else:
                        term_x = np.kron(term_x, self.I[k])
                        term_y = np.kron(term_y, self.I[k])
                        term_z = np.kron(term_z, self.I[k])
                #add this term to the Hamiltonian
                self.H -= self.J[i,j]*(term_x + term_y + term_z)
        #Include linear coupling to the magnetic field.
        if self.B is not None:
            for i in range(self.N):
                if i == 0:
                    B_term = self.B[0,0]*self.Sx[0]+self.B[0,1]*self.Sy[0]+self.B[0,2]*self.Sz[0]
                else:
                    B_term = self.I[0]
                for k in range(1, self.N):
                    if k == i:
                        B_term = np.kron(B_term, self.B[k,0]*self.Sx[k]+self.B[k,1]*self.Sy[k]+self.B[k,2]*self.Sz[k])
                    else:
                        B_term = np.kron(B_term, self.I[k])
                self.H -= B_term

    def get_energies(self):
        """
        Calculate the energy levels and eigenstates of the Hamiltonian
        """
        self.energies, self.states = np.linalg.eigh(self.H)

    def count_multiplicities(self, tolerance = 0.0001):
        """
        Create lists of multiplicities and unrepeated energy levels.
        """
        if hasattr(self, 'energies'):
            self.multiplicities = []
            #levels without degeneracy
            self.reduced_energies = []

            val = self.energies[0]
            count = 0
            #find the multiplicities
            for i in range(self.energies.shape[0]):
                if abs(self.energies[i] - val) < tolerance:
                    count += 1
                else:
                    self.multiplicities.append(count)
                    self.reduced_energies.append(self.energies[i-1])
                    count = 1
                    val = self.energies[i]
            self.multiplicities.append(count)
            self.reduced_energies.append(self.energies[i])
        #create the class attributes that
        #count_multiplicities depends upon
        else:
            self.get_energies()
            self.count_multiplicities()



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
        ascending = True
        if even:
            k = n+2
        else:
            k = n+1
        while not truncate:
            if n != k:
                if k < 2 and ascending == False:
                    truncate = True

                else:
                    I, I_err = quad(combine_functions(n,k),limits[0],limits[1])
                    E_nk = E_diff(n,k)
        
                    if abs(I/E_nk) < tolerance:
                        if n > 2 and ascending == True:
                            if even:
                                k = n-2
                            else:
                                k = n-1
                            ascending = False

                        else:
                            truncate = True
                    else:
                        I_list.append(I/E_nk)
                        k_list.append(k)

                        if ascending == True:
                            k += 2
                        elif ascending == False:
                            k -= 2

            else:
                if ascending == True:
                    k += 2
                elif ascending == False:
                    k -= 2

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
