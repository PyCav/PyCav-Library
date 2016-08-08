Quantum
====================================

Introduction to the Quantum module
--------------------------

This module contains functions for use in 1st order perturbation theory calculations and for solving 1d boundary value problems using the Shooting method. It also contains a class designed to represent systems of interacting spins. Documentation explaining the use of each function and the algorithms used within will be presented in individual sections. 

In-depth Documentation
-------

.. toctree::
   :maxdepth: 1

   numerov
   bisection
   perturbation
   spinsystem
   plotting

The list below summarises the functions, their input arguments and their outputs for quick reference for the informed user:

Functions
---------

numerov(x,dx,V,E,initial_values,params)
^^^^^^^^^^^^^^^^^^^

   This function performs a Numerov integration at the given energy within the given domain and returns the un-normalised wavefunction evaluated at the given spatial coords

   **Parameters:**

   *x: numpy array*

   An N element numpy array of equally spaced points in space (creating using numpy linspace is advised) at which the wavefunction will be evaluated

   *dx: float*

   Must give the spacing between points in the x array
   
   *V: function*
   
   Pass a function which takes x as an argument and returns the value of potential at that point, V(x)
   
   *E: float*
   
   The energy of the time independent wavefunction. This will give exponential growth if E does not correspond to a bound or free state energy eigenvalue.
   
   *inital_value: list*
   
   A list with 2 elements. The first of these is the boundary condition for the wavefunction at the lowest value spatial coordinate. The second of these is used to initialise the next point in from this, it must be non-zero but apart from this requirement it only effects normalisation (see in detail documentation)
   
   *params: list*
   
   List which can be used within your code to hold various physical parameter. The first element is required to be equal to the particle mass but apart from this size of the list and the other parameters are not called
   
   **Returns:**

   A *numpy array* containing the approximated wavefunction evaluated at x


bisection_search(x,dx,V,params,bracket_E,tolerance = 0.5,max_evals = 1000)
^^^^^^^^^^^^^^^^^^^^^^^^^^
   
   Uses the Numerov method to perform a bisection search in energy for a wavefunction which goes to zero at the boundaries. Shares similar arguments with the numerov function above see *x, dx, V, params*

   **Unique Parameters:**

   *bracket_E: list*

   A 2 element list which expresses a range in which the energy eigenvalue(s) lie. For the 2 energy values in this list, one of the values will have the Numerov approximated wavefunction above 0 and the other below 0. The ordering of these is handled by the function.

   *tolerance: float*

   The tolerance of the bisection search i.e. if the absolute value of the wavefunction at right hand boundary is less than the tolerance then the search is complete.

   *max_evals: int*

   The number of search evaluations taken before the search is given up. It is more likely your bracket_E list is not wide enough or your tolerance is too low than max_evals is too low
   
   **Returns:**

   The un-normalised wavefunction, as a numpy array, which has satisfied the bisection search (or max_evals has been reached) and the energy eigenvalue estimate as a float


first_order_energy_sft(n,H,unperturb_wf,params,limits = [-np.inf,np.inf])
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   
   Works out the 1st order energy shift from time independent perturbation theory for a given unperturbed system and the applied perturbation

   **Parameters:**

   *n: int or list*

   Principal quantum number (or list of these) which labels the unperturbed wavefunctions
   
   *H: function*

   The applied perturbation as a function of position 
   
   *unperturb_wf: function*

   A function which is passed params and n and returns a function of position e.g.
   
   .. code-block:: python
	
    def unperturb_wf(params,n):
      a = params[0]
      m = n+1
      def psi_n(x):
        return np.sqrt(2./a)*np.sin((m*np.pi/a)*(x+a/2.))
      return psi_n
      
   Here psi_n(x) is the returned function of position
  
   *params: list*
   
   List which can be used within your code to hold various physical parameter used by unperturb_wf and other functions (see later)
  
   *limits: list*
  
   List containing the integration limits for the inner product of the wavefunctions and the perturbation. Default is the whole space. For hard wall potentials adjust these limits to respect the boundaries

   **Returns:**

   A list or float depending on the input argument n containing the first energy shifts to the unperturbed wavefunctions for the given quantum numbers in n

first_order_wf(n,H,unperturb_wf,unperturb_erg,params,tolerance = 0.01, limits = [-np.inf,np.inf], return_list = False)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   Calculates the 1st order perturbed wavefunction for a given unperturbed system and the applied perturbation. The system is defined by its known unperturbed wavefunctions and energies. This function takes similar parameters as first_order_energy_sft (see above) so only new parameters will be defined. Further documentation can be found on the functions documentation page

   **Unique Parameters:**

   *unperturb_erg: function*

   A function which is passed params and n and returns the energy of the n-th unperturbed eigenstate e.g. for a harmonic oscillator

   .. code-block:: python
  
    def unperturb_erg(params,n):
      return (n+0.5)*params[1]

   where params[1] contains the angular frequency (for hbar = 1)

   *tolerance: float*

   The value below which terms in the 1st order wavefunction sum are ignored

   *return_list: float*

   Set to True if you require the perturbation sum prefactors and values of the principal quantum numbers of the unperturbed wavefunctions

   **Returns:**

   A function of position which corresponds to the 1st order perturbed wavefunction and if return_list = True, copies of the principal quantum number lists and the sum prefactors list which were used to calculate the resultant perturbed wavefunction

create_Su(s)
^^^^^^^^^^^
   Returns the 'up' angular momentum ladder operator in the z basis.

   **Parameters:**

   *s: int or float*

   The angular momentum quantum number, 1/2 for isolated electrons, 1 for a hydrogenic p-orbital etc. If s is a float, it must be a half-integer.

   **Returns:**

   *numpy array with complex elements* 

create_Sd(s)
^^^^^^^^^^^
   Returns the 'down' angular momentum ladder operator in the z basis. Usage is identical to create_Su.

create_Sx(s)
^^^^^^^^^^^
   Returns the 'x' angular momentum operator in the z basis. Usage is identical to create_Su.

create_Sy(s)
^^^^^^^^^^^
   Returns the 'y' angular momentum operator in the z basis. Usage is identical to create_Su.

create_Sz(s)
^^^^^^^^^^^
   Returns the 'z' angular momentum operator in the z basis. Usage is identical to create_Su.


